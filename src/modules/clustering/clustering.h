// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

#pragma once

#include "analyser/siteFilter.h"
#include "analyser/siteSelector.h"
#include "base/processPool.h"
#include "classes/configuration.h"
#include "classes/coreData.h"
#include "generator/context.h"
#include "main/dissolve.h"
#include "module/context.h"
#include "module/module.h"
#include <unordered_set>

class Atom;

// Clustering Module
class ClusteringModule : public Module
{
    public:
    ClusteringModule();
    ~ClusteringModule() override = default;

    /*
     * Definition
     */
    private:
    // Target configurations
    Configuration *targetConfiguration_{nullptr};

    // Simulation data...
    // Struct
    struct BondInfo
    {
        const SpeciesSite *a_, *b_;
        double cutOff;

        BondInfo(const SpeciesSite *site1, const SpeciesSite *site2, double cut) : a_(site1), b_(site2), cutOff(cut) {}
    };

    // GUI site and cutoff selection
    BondInfo userSites_{nullptr, nullptr, -1};
    // Vector of selected cluster definitions (when multiple can be defined)
    std::vector<BondInfo> selectedBonds_;
    // Map of neighbours for every site (symmetric)
    Analyser::SiteMap neighbourMap_{};
    // Map of every member in every cluster
    std::map<int, std::vector<const Site *>> clusterMap_{};
    // Map of cluster ID to size of cluster
    std::map<int, int> sizeDistribution_{};
    // Map of cluster ID to mass of cluster
    std::map<int, float> clusterMasses_{};
    // Map of mass to vector of cluster IDs (useful for processing)
    std::map<float, std::vector<int>> massDistribution_{};
    // {Site A : {Site B, C ... : Coordination Number for each pair}} where Site A - B, C ... are all the combinations defined
    std::map<const SpeciesSite *, std::map<const SpeciesSite *, float>> clusterSpeciesCoordNo_;
    // Minimum cluster size to count in coordination number calc
    int minCNSize_{0};
    // Maximum cluster size to count in coordination number calc
    int maxCNSize_{0};
    // Map of cluster ID to radius of gyration
    std::map<int, float> radiusOfGyration_;
    // Minimum size to include in radius of gyration calculation
    int gyrationMinSize_{10};
    // Map of cluster ID to CoM vector from reference site (first member in cluster map) - Required in Rg calc. Might be useful
    // down the road
    std::map<int, Vec3<double>> clusterCoM_;
    // Fractal dimension from linear regression of LogLog Radius of gyration - cluster mass
    double fractalDimension_{-1};
    // Minimum cluster size to display
    int minClusterForConfig_{5};
    // Maximum cluster size to display
    int maxClusterForConfig_{10000};
    // Tracks if the clusters have been generated yet
    bool clustersConfigGenerated_{false};
    // Pointer to the cluster display configuration
    Configuration *clusterConfig_{nullptr};

    /*
     * Processing
     */

    private:
    // Run main processing
    Module::ExecutionResult process(ModuleContext &moduleContext) override;
    // Recursion for cluster generation
    void buildCluster(const Site *startSite, std::unordered_set<const Site *> &visited);

    public:
    // Set up module for processing
    bool setUp(ModuleContext &moduleContext, Flags<KeywordBase::KeywordSignal> actionSignals) override;
    // Generation of the cluster visualisation configuration
    Configuration *generateClustersConfig(Dissolve *dissolve, Configuration *source);
    // Getter for the target configuration of the module
    Configuration *getSourceConfig();
    // Radius of Gyration: computes "compactness" of each cluster (used in fractal dimension as well)
    void generateRadiusOfGyration();
    // Fractal dimension: Performs linear regression on log(Radius of Gyration) - log(cluster mass), returns gradient
    void generateFractalDimension();
    // Ensures we don't visualise the config before it's ready (still necessary?)
    bool viewingReady{false};
};
