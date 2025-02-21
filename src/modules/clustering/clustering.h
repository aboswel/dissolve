// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

// Define the main clustering class to store user defined parameters and analysis settings, as well as function prototypes and helper methods
// clustering.h
#pragma once

#include "module/module.h"
#include "module/context.h"
#include "analyser/siteFilter.h"
#include "analyser/siteSelector.h"
#include <unordered_set>

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
    Configuration *targetConfiguration_{nullptr}; // Will set in GUI
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
    std::vector<BondInfo> selectedBonds_;
    Analyser::SiteMap neighbourMap_;
    std::map<int, std::vector<const Site*>> clusterMap_;
    std::map<int, int> sizeDistribution_;
    std::map<float, std::vector<int>> massDistribution_;
    std::map<const SpeciesSite*, std::map<const SpeciesSite*, float>> clusterSpeciesCoordNo_;
    int minCNSize_{0};
    int maxCNSize_{0};
    std::map<int, float> radiusOfGyration_;
    int gyrationMinSize_{3};

    // Getters
    public:
    std::vector<BondInfo>& getSelectedBonds();
    Analyser::SiteMap& getNeighbourMap();
    std::map<int, std::vector<const Site*>>& getClusterMap();
    std::map<int, int>& getSizeDistribution();
    std::map<float, std::vector<int>>& getMassDistribution(); // Use .size() to get the number of clusters of a mass.
    std::map<const SpeciesSite*, std::map<const SpeciesSite*, float>>& getClusterSpeciesCoordNo();
    std::map<int, float>& getRadiusOfGyration();

    /*
     * Processing
     */
    public:
    // Set up module for processing
    // bool setUp();

    private:
    // Run main processing
    Module::ExecutionResult process(ModuleContext &moduleContext) override;

    // Generates a symmetric adjacency list from specified intermolecular bonds
    void generateNeighbourMap();
    // Generates a cluster map from adjacency list
    void generateClusterMap();

    // Calculates cluster size (No. of members) distribution from cluster map
    void generateSizeDistribution();
    // Calculates cluster mass distribution from cluster map
    void generateMassDistribution();
    // Calculates clustering coordination numbers symmetrically for each intermolecular bond. Min and Max cluster sizes to include in calc. can be specified.
    void  generateClusterSpeciesCoordNo();
    
    // Radius of Gyration: computes "compactness" of each cluster (used in fractal dimension as well)
    void generateRadiusOfGyration();
};
