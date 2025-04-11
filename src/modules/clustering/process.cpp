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
    return true;
}

Module::ExecutionResult ClusteringModule::process(ModuleContext &moduleContext)
{
    auto &moduleData = moduleContext.dissolve().processingModuleData();

    if (!(a_ && b_ && (cutoff_ > 0)))
    {
        Messenger::error("Cluster definition invalid!");
        return ExecutionResult::Failed;
    }

    // Process code
    // Produce NeighbourMap
    neighbourMap_.clear();
    Analyser::SiteMap neighbourMapA, neighbourMapB;

    // Transform SiteObjects to instances ready for the filter function
    SiteSelector selectionA(targetConfiguration_, std::vector<const SpeciesSite *>{a_});
    const Analyser::SiteVector &siteVectorA = selectionA.sites();

    SiteSelector selectionB(targetConfiguration_, std::vector<const SpeciesSite *>{b_});
    const Analyser::SiteVector &siteVectorB = selectionB.sites();

    SiteFilter filterA(targetConfiguration_, siteVectorA);
    std::tie(std::ignore, neighbourMapA) = filterA.filterBySiteProximity(siteVectorB, Range(0.001, cutoff_), 1,
                                                                         100); // min of 0.01 to avoid self-selection

    // If we're dealing with a definition between the same site, it's already symmetric. Running twice will add dupes
    if (a_ != b_)
    {
        SiteFilter filterB(targetConfiguration_, siteVectorB);
        std::tie(std::ignore, neighbourMapB) = filterB.filterBySiteProximity(siteVectorA, Range(0.001, cutoff_), 1, 100);
    }

    // Combining the neighbour maps into a single map. Because keys may already exist, need to check for them and add
    // neighbours if exists.
    for (auto neighbourMap : {neighbourMapA, neighbourMapB})
        for (const auto &[site, neighbours] : neighbourMap)
        {
            if (neighbourMap_.contains(site))
                neighbourMap_[site].insert(neighbourMap_[site].end(), neighbours.begin(), neighbours.end());
            else
                neighbourMap_.insert({site, neighbours});
        }

    // ClusterMap generation 2.0
    clusterMap_.clear();
    std::unordered_set<const Site *> processed;
    auto clusterTrack = 1;
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

    // Cluster size distribution
    sizeDistribution_.clear();
    for (const auto &[clusterID, members] : clusterMap_)
        sizeDistribution_[members.size()].emplace_back(clusterID);

    auto &sizeData = moduleData.realise<Data1D>("SizeDist", name());
    sizeData.clear();
    for (const auto &[size, mems] : sizeDistribution_)
        sizeData.addPoint(std::log(size), std::log(mems.size()));

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

    auto &massData = moduleData.realise<Data1D>("MassDist", name());
    massData.clear();
    for (const auto &[mass, mems] : massDistribution_)
        massData.addPoint(std::log(mass), std::log(mems.size()));

    // Generation of radius of gyration
    radiusOfGyration_.clear();
    const auto *box = targetConfiguration_->box();
    for (const auto &[clusterID, clusterVec] : clusterMap_)
    {
        if (clusterVec.size() < gyrationMinSize_)
            continue;

        // CoM mass weighted calc from reference site
        Vec3<double> massWeightedTotalVec{0, 0, 0};
        const auto *refSite{clusterVec[0]}; // Reference as first member in cluster
        for (const auto &clusterMem : clusterVec)
            massWeightedTotalVec +=
                (box->minimumVector(refSite->origin(), clusterMem->origin())) * clusterMem->parent()->parent()->mass();

        massWeightedTotalVec /= clusterMasses_[clusterID];
        clusterCoM_[clusterID] = massWeightedTotalVec;
        auto massWeightedDistanceSqrd = 0.0;

        // Run through again for mass weighted distance squared
        for (const auto &clusterMem : clusterVec)
            massWeightedDistanceSqrd +=
                (box->minimumDistanceSquared(box->minimumVector(clusterMem->origin(), refSite->origin()),
                                             clusterCoM_[clusterID])) *
                clusterMem->parent()->parent()->mass();

        radiusOfGyration_[clusterID] = std::sqrt(massWeightedDistanceSqrd / clusterMasses_[clusterID]);
    }

    // Fractal Dimension: Create a Data1D object of the log log plot, perform linear regression, return gradient
    Data1D loglog;
    loglog.initialise(radiusOfGyration_.size(), false); // Add errors to this at some point?

    // Generate Data1D
    for (const auto &[clusterID, rg] : radiusOfGyration_)
        loglog.addPoint(std::log(radiusOfGyration_[clusterID]), std::log(clusterMasses_[clusterID]));

    // Perform linear regression
    fractalDimension_ = (Regression::linear(loglog, radiusOfGyration_.size()));

    viewingReady = true;

    // Output for debugging
    if (false)
    {
        LineParser parser;
        std::string filename = std::format("{}.neighbourdata.txt", targetConfiguration_->niceName());

        // Output file for diagnostics
        if (!parser.openOutput(filename))
        {
            Messenger::error("Failed to open output file '{}'.\n", filename);
            return ExecutionResult::Failed;
        }

        // Write header with site information
        parser.writeLineF("# Analysis for sites: {} - {}\n", a_->parent()->name(), b_->parent()->name());

        // Write clusters and other diagnostics...
        // Now compute the cluster size distribution and output it
        parser.writeLineF("\n=== Cluster Size Distribution ===\n");
        for (const auto &[clusterSize, mems] : sizeDistribution_)
            parser.writeLineF("  Cluster Size {} : {} clusters\n", clusterSize, mems.size());

        // Compute mass distribution diagnostic
        parser.writeLineF("\n=== Cluster Mass Distribution ===\n");
        for (const auto &[clusterMass, clusterVec] : massDistribution_)
            for (const auto &cluster : clusterVec)
                parser.writeLineF("  Cluster Mass {:.3f} : cluster ID: {}\n", clusterMass, cluster);

        // Write radius of gyration diagnostics to output file
        parser.writeLineF("\n=== Radius of Gyration ===\n");
        for (const auto &[clusterID, radius] : radiusOfGyration_)
            parser.writeLineF("Cluster {}: {:.3f}\n", clusterID, radius);

        // Write fractal dimension diagnostics to output file
        parser.writeLineF("\n=== Fractal Dimension ===\n");
        parser.writeLineF("Fractal dimension: {:.3f}\n", fractalDimension_);
    }

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

Configuration *ClusteringModule::generateClustersConfig(Dissolve *dissolve, Configuration *source, int displaySize,
                                                        int displayID)
{
    // If the config does not exist at all, make it
    if (!dissolve->coreData().findConfiguration("clusters"))
    {
        clusterConfig_ = dissolve->coreData().addConfiguration();
        clusterConfig_->setName("clusters");
        clusterConfig_->createBox(source->box()->axisLengths(), source->box()->axisAngles());
    }
    // If the configuration exists, empty it
    else
        clusterConfig_->empty();
    // Display all clusters
    if (displaySize == 0)
    {
        for (const auto &[clusterID, mems] : clusterMap_)
            for (const auto &site : mems)
            {
                auto copyableMolecule = std::make_shared<Molecule>(*site->molecule());
                clusterConfig_->copyMolecule(copyableMolecule);
            }
    }
    // Display all clusters of size displaySize
    else if (displaySize != 0 && displayID == 0)
    {
        for (const auto &[clusterID, mems] : clusterMap_)
            if (mems.size() == displaySize)
                for (const auto &site : mems)
                {
                    auto copyableMolecule = std::make_shared<Molecule>(*site->molecule());
                    clusterConfig_->copyMolecule(copyableMolecule);
                }
    }
    // Display cluster with ID displayID
    else if (displaySize != 0 && displayID != 0)
    {
        for (const auto &site : clusterMap_[displayID])
        {
            auto copyableMolecule = std::make_shared<Molecule>(*site->molecule());
            clusterConfig_->copyMolecule(copyableMolecule);
        }
    }

    clusterConfig_->updateObjectRelationships();
    clusterConfig_->incrementContentsVersion();

    if (clusterConfig_->nAtoms() == 0)
        Messenger::error("No displayable clusters were found! Tweak min/max values.");
    else
        Messenger::print("clusters config successfully generated");

    return clusterConfig_;
}

Configuration *ClusteringModule::getSourceConfig() { return targetConfiguration_; }
