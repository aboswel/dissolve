// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2024 Team Dissolve and contributors

// Define the main clustering class to store user defined parameters and analysis settings, as well as function prototypes and helper methods
// clustering.h
#pragma once

#include "module/module.h"
#include "module/context.h"
#include "analyser/siteFilter.h"
#include "analyser/siteSelector.h"

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
    // Need to retrieve all the species in a configuration, then iterate through them to get the sites. Site combinations are then selected in the gui
    // Need a struct which contains all the user specified info for the intermolecular bonding
    
    // Struct 
    struct BondInfo 
    {
    const SpeciesSite *a_, *b_;
    double cutOff;
    
        BondInfo(const SpeciesSite* site1,
                 const SpeciesSite* site2,
                 double cut)
        : a_(site1), b_(site2), cutOff(cut) {}
    };

    // Vector of intermolecular bonds specified
    std::vector<BondInfo> selectedBonds;

    /*
     * Processing
     */
    public:
    // Set up module for processing
    // bool setUp();

    private:
    // Run main processing
    Module::ExecutionResult process(ModuleContext &moduleContext) override;
    // Basic analysis routine
    // Generates a symmetric adjacency list from specified intermolecular bonds
    std::pair<Analyser::SiteVector, Analyser::SiteMap> siteFiltering(Configuration *cfg_, std::vector<BondInfo> bonds);
    // Generates a cluster map from adjacency list
    std::map<int, std::vector<const Site*>> generateClusters(Analyser::SiteMap neighbourMap);
    // Basic metric computation
    // Calculates cluster size (No. of members) distribution from cluster map
    std::map<int, int> sizeDistribution(std::map<int, std::vector<const Site*>> clusterMap);
    // Calculates cluster mass distribution from cluster map
    std::map<float, int> massDistribution(std::map<int, std::vector<const Site*>> clusterMap);
};
