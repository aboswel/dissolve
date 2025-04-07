// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

#include "analyser/typeDefs.h"
#include "base/sysFunc.h"
#include "generator/box.h"
#include "math/regression.h"
#include "modules/clustering/clustering.h"

bool ClusteringModule::setUp(ModuleContext &moduleContext, Flags<KeywordBase::KeywordSignal> actionSignals)
{
    // Existing oldConfigurations in the restart file generate a new tab which is problematic. Only want this to run we've
    // freshly loaded a save
    if (!clusterConfig_)
    {
        auto *oldConfig = moduleContext.dissolve().coreData().findConfiguration("clusters");
        if (oldConfig)
        {
            moduleContext.dissolve().coreData().removeConfiguration(oldConfig);
            Messenger::print("Removed existing cluster configuration");
        }
    }
    // Clear all existing data
    selectedBonds_.clear();
    return true;
}

Module::ExecutionResult ClusteringModule::process(ModuleContext &moduleContext)
{
    clustersConfigGenerated_ = false;
    auto &moduleData = moduleContext.dissolve().processingModuleData();
    LineParser parser;
    std::string filename = std::format("{}.neighbourdata.txt", targetConfiguration_->niceName());

    if (!(userSites_.a_ && userSites_.b_ && (userSites_.cutOff > 0)))
    {
        Messenger::error("Cluster definition invalid!");
        return ExecutionResult::Failed;
    }

    selectedBonds_.emplace_back(userSites_);

    // Output file for diagnostics
    if (!parser.openOutput(filename))
    {
        Messenger::error("Failed to open output file '{}'.\n", filename);
        return ExecutionResult::Failed;
    }

    // Write header with site information
    parser.writeLineF("# Analysis for sites: {} - {}\n", selectedBonds_[0].a_->parent()->name(),
                      selectedBonds_[0].b_->parent()->name());

    // Process code

    // Produce NeighbourMap
    neighbourMap_.clear();
    for (auto bond : selectedBonds_)
    {
        Analyser::SiteMap neighbourMapA, neighbourMapB;

        // Transform SiteObjects to instances ready for the filter function
        SiteSelector selectionA(targetConfiguration_, std::vector<const SpeciesSite *>{bond.a_});
        const Analyser::SiteVector &siteVectorA = selectionA.sites();

        SiteSelector selectionB(targetConfiguration_, std::vector<const SpeciesSite *>{bond.b_});
        const Analyser::SiteVector &siteVectorB = selectionB.sites();

        SiteFilter filterA(targetConfiguration_, siteVectorA);
        std::tie(std::ignore, neighbourMapA) =
            filterA.filterBySiteProximity(siteVectorB, Range(0.001, bond.cutOff), 1, 100); // min of 0.01 to avoid self-bonding

        // If we're dealing with a bond between the same site, it's already symmetric. Running twice will add dupes
        if (bond.a_ != bond.b_)
        {
            SiteFilter filterB(targetConfiguration_, siteVectorB);
            std::tie(std::ignore, neighbourMapB) =
                filterB.filterBySiteProximity(siteVectorA, Range(0.001, bond.cutOff), 1, 100);
        }

        // Combining the neighbour maps into a single map. Because keys may already exist, need to check for them and add
        // neighbours if exists.
        for (auto neighbourMap : {neighbourMapA, neighbourMapB})
        {
            for (const auto &[site, neighbours] : neighbourMap)
            {
                if (neighbourMap_.contains(site))
                    neighbourMap_[site].insert(neighbourMap_[site].end(), neighbours.begin(), neighbours.end());
                else
                    neighbourMap_.insert({site, neighbours});
            }
        }
    }

    // ClusterMap generation 2.0
    {
        clusterMap_.clear();
        std::unordered_set<const Site *> processed;
        int clusterTrack{1};
        for (const auto &[clusterStart, _] : neighbourMap_)
        {
            if (!processed.contains(clusterStart))
            {
                std::unordered_set<const Site *> visited{clusterStart};
                buildCluster(clusterStart, visited);
                processed.insert(visited.begin(), visited.end());
                clusterMap_.insert({clusterTrack, std::vector<const Site *>(visited.begin(), visited.end())});
                clusterTrack++;
            }
        }
    }

    // Cluster size distribution
    sizeDistribution_.clear();
    for (const auto &[_, members] : clusterMap_)
        sizeDistribution_[members.size()]++;

    // Cluster mass calculation
    clusterMasses_.clear();
    for (const auto &[clusterID, memberVec] : clusterMap_)
    {
        float clusterMass{0};
        for (const auto &member : memberVec)
            clusterMass += member->parent()->parent()->mass();

        clusterMasses_[clusterID] = clusterMass;
    }

    // Cluster mass distribution
    massDistribution_.clear();
    for (const auto &[clusterID, clusterMass] : clusterMasses_)
        massDistribution_[clusterMass].emplace_back(clusterID);

    // Generation of radius of gyration
    {
        radiusOfGyration_.clear();
        const Box *box = targetConfiguration_->box();
        // Iterate through cluster map, skip clusters below size min
        for (const auto &[clusterID, clusterVec] : clusterMap_)
        {
            if (clusterVec.size() < gyrationMinSize_)
                continue;

            // Now calculate the centre of mass of the given cluster with regards to a reference site (first in clustermap)
            // Collect the coordinates of each member, multiply by mass of parent, accumlate a total, then divide by the mass of
            // the cluster
            Vec3<double> massWeightedTotalVec{0, 0, 0};
            const Site *refSite{clusterVec[0]}; // Define a reference site (the first member of the cluster)
            for (const auto &clusterMem : clusterVec)
                // Accumlate mass weighted total vector from reference site
                massWeightedTotalVec +=
                    (box->minimumVector(refSite->origin(), clusterMem->origin())) * clusterMem->parent()->parent()->mass();

            massWeightedTotalVec /= clusterMasses_[clusterID];
            clusterCoM_[clusterID] = massWeightedTotalVec;
            double massWeightedDistanceSqrd{0};
            // Now time to calculate the Radius of Gyration
            // Need to run through the clusterMap again, get the mim sqrd distance of site from CoM using ref as origin
            for (const auto &clusterMem : clusterVec)
                massWeightedDistanceSqrd +=
                    (box->minimumDistanceSquared(box->minimumVector(clusterMem->origin(), refSite->origin()),
                                                 clusterCoM_[clusterID])) *
                    clusterMem->parent()->parent()->mass();

            radiusOfGyration_[clusterID] = std::sqrt(massWeightedDistanceSqrd / clusterMasses_[clusterID]);
        }
    }

    // Calculate the fractal dimension from radius of gyration
    {
        // Create a Data1D object of the log log plot, perform linear regression, return gradient
        Data1D loglog;
        loglog.initialise(radiusOfGyration_.size(), false); // Add errors to this at some point?

        // Generate Data1D
        for (const auto &[clusterID, rg] : radiusOfGyration_)
            loglog.addPoint(std::log(radiusOfGyration_[clusterID]), std::log(clusterMasses_[clusterID]));

        // Perform linear regression
        fractalDimension_ = (Regression::linear(loglog, radiusOfGyration_.size()));
    }

    viewingReady = true;

    // Write clusters and other diagnostics...
    // Now compute the cluster size distribution and output it
    parser.writeLineF("\n=== Cluster Size Distribution ===\n");
    for (const auto &[clusterSize, count] : sizeDistribution_)
        parser.writeLineF("  Cluster Size {} : {} clusters\n", clusterSize, count);

    // Compute mass distribution diagnostic
    parser.writeLineF("\n=== Cluster Mass Distribution ===\n");
    for (const auto &[clusterMass, clusterVec] : massDistribution_)
    {
        for (const auto &cluster : clusterVec)
            parser.writeLineF("  Cluster Mass {:.3f} : cluster ID: {}\n", clusterMass, cluster);
    }

    // Write radius of gyration diagnostics to output file
    parser.writeLineF("\n=== Radius of Gyration ===\n");
    for (const auto &[clusterID, radius] : radiusOfGyration_)
        parser.writeLineF("Cluster {}: {:.3f}\n", clusterID, radius);

    // Write fractal dimension diagnostics to output file
    parser.writeLineF("\n=== Fractal Dimension ===\n");
    parser.writeLineF("Fractal dimension: {:.3f}\n", fractalDimension_);

    return ExecutionResult::Success;
}

void ClusteringModule::buildCluster(const Site *startSite, std::unordered_set<const Site *> &visited)
{
    for (const auto &[neighbour, _] : neighbourMap_[startSite])
    {
        if (!visited.contains(neighbour))
        {
            visited.emplace(neighbour);
            buildCluster(neighbour, visited);
        }
    }
}

Configuration *ClusteringModule::generateClustersConfig(Dissolve *dissolve, Configuration *source)
{
    // If the cluster config has already been generated this iteration, just return it - Is this still necessary?
    if (clustersConfigGenerated_)
        return clusterConfig_;

    // If the config does not exist at all, make it
    if (!clusterConfig_)
    {
        clusterConfig_ = dissolve->coreData().addConfiguration();
        clusterConfig_->setName("clusters");
        clusterConfig_->createBox(source->box()->axisLengths(), source->box()->axisAngles());
    }
    // If the configuration exists, but needs updating, empty it
    else
        clusterConfig_->empty();

    // Add the molecules we want.
    for (const auto &[clusterID, mems] : clusterMap_)
    {
        if (mems.size() >= minClusterForConfig_ && mems.size() <= maxClusterForConfig_)
        {
            for (const auto &site : mems)
            {
                auto copyableMolecule = std::make_shared<Molecule>(*site->molecule());
                clusterConfig_->copyMolecule(copyableMolecule);
            }
        }
    }
    clusterConfig_->updateObjectRelationships();
    clusterConfig_->incrementContentsVersion();

    if (clusterConfig_->nAtoms() == 0)
        Messenger::error("No displayable clusters were found! Tweak min/max values.");
    else
        Messenger::print("clusters config successfully generated");

    clustersConfigGenerated_ = true;
    return clusterConfig_;
}

Configuration *ClusteringModule::getSourceConfig() { return targetConfiguration_; }
