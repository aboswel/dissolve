// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

#include "analyser/typeDefs.h"
#include "base/sysFunc.h"
#include "classes/localMolecule.h"
#include "data/elements.h"
#include "generator/box.h"
#include "generator/copy.h"
#include "math/regression.h"
#include "modules/clustering/clustering.h"

bool ClusteringModule::setUp(ModuleContext &moduleContext, Flags<KeywordBase::KeywordSignal> actionSignals)
{
    // Check user definitions
    if (!(a_ && b_ && (cutoff_ > 0)))
        Messenger::error("Cluster definition invalid!");

    return true;
}

Module::ExecutionResult ClusteringModule::process(ModuleContext &moduleContext)
{
    auto &moduleData = moduleContext.dissolve().processingModuleData();

    // Produce NeighbourMap
    neighbourMap_.clear();
    Analyser::SiteMap neighbourMapA, neighbourMapB;

    // Transform SiteObjects to instances ready for the filter function
    SiteSelector selectionA(targetConfiguration_, std::vector<const SpeciesSite *>{a_});
    const auto &siteVectorA = selectionA.sites();

    SiteSelector selectionB(targetConfiguration_, std::vector<const SpeciesSite *>{b_});
    const auto &siteVectorB = selectionB.sites();

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

    if (neighbourMap_.empty())
        Messenger::error("No neighbours found!");

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

    if (saveSizeDist_)
    {
        LineParser parser;
        parser.openOutput(std::format("{}.sizedist.txt", targetConfiguration_->niceName()));
        parser.writeLineF("# Analysis for sites: {} - {}\n", a_->parent()->name(), b_->parent()->name());
        parser.writeLineF("\n=== Cluster Size Distribution ===\nCluster Size : no. of clusters\n");
        for (const auto &[clusterSize, mems] : sizeDistribution_)
            parser.writeLineF("{} : {}\n", clusterSize, mems.size());
    }
    if (saveMassDist_)
    {
        LineParser parser;
        parser.openOutput(std::format("{}.massdist.txt", targetConfiguration_->niceName()));
        parser.writeLineF("# Analysis for sites: {} - {}\n", a_->parent()->name(), b_->parent()->name());
        parser.writeLineF("\n=== Cluster Mass Distribution ===\nCluster Mass : no. of clusters\n");
        for (const auto &[clusterMass, mems] : massDistribution_)
            parser.writeLineF("{:.3f} : {}\n", clusterMass, mems.size());
    }
    if (saveRgMass_)
    {
        LineParser parser;
        parser.openOutput(std::format("{}.massRg.txt", targetConfiguration_->niceName()));
        parser.writeLineF("# Analysis for sites: {} - {}\n", a_->parent()->name(), b_->parent()->name());
        parser.writeLineF("\n=== Mass - Radius of gyration ===\nCluster Mass : Radii\n");
        for (const auto &[clusterID, radius] : radiusOfGyration_)
            parser.writeLineF("{} : {}\n", clusterMasses_[clusterID], radius);
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

// This ends up running (at least) four times after each iteration?
void ClusteringModule::generateClustersConfig(Dissolve &dissolve, int displaySize, int displayID)
{
    if (clusterConfig_.generator().node("clusters"))
        clusterConfig_.clear();
    else
        clusterConfig_.setName("clusters");

    // Molecule transfer only works with a generator
    clusterConfig_.generator().createRootNode<CopyGeneratorNode>("clusters", targetConfiguration_);
    clusterConfig_.generate({dissolve.worldPool(), dissolve});
    clusterConfig_.removeMolecules(clusterConfig_.molecules());

    // Display all clusters
    if (displaySize == 0)
    {
        for (const auto &[clusterID, mems] : clusterMap_)
            for (const auto &site : mems)
            {
                auto mol = clusterConfig_.addMolecule(site->parent()->parent());
                for (auto &&[molAtom, sourceAtom] : zip(mol->atoms(), site->molecule()->atoms()))
                {
                    molAtom->setCoordinates(sourceAtom->r());
                    clusterConfig_.updateAtomLocation(molAtom);
                }
            }
    }

    // Display all clusters of size displaySize
    else if (displaySize != 0 && displayID == 0)
    {
        for (const auto &[clusterID, mems] : clusterMap_)
            if (mems.size() == displaySize)
                for (const auto &site : mems)
                {
                    auto mol = clusterConfig_.addMolecule(site->parent()->parent());
                    for (auto &&[molAtom, sourceAtom] : zip(mol->atoms(), site->molecule()->atoms()))
                    {
                        molAtom->setCoordinates(sourceAtom->r());
                        clusterConfig_.updateAtomLocation(molAtom);
                    }
                }
    }
    // Display cluster with ID displayID
    else if (displaySize != 0 && displayID != 0)
    {
        for (const auto &site : clusterMap_[displayID])
        {
            auto mol = clusterConfig_.addMolecule(site->parent()->parent());
            for (auto &&[molAtom, sourceAtom] : zip(mol->atoms(), site->molecule()->atoms()))
            {
                molAtom->setCoordinates(sourceAtom->r());
                clusterConfig_.updateAtomLocation(molAtom);
            }
        }
    }

    clusterConfig_.incrementContentsVersion();
    clusterConfig_.updateObjectRelationships();

    if (clusterConfig_.nAtoms() == 0)
    {
        Messenger::error("No clusters!");
    }
    else
        Messenger::print("Cluster visualisation generated");
}
