// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

// process.cpp
#include "modules/clustering/clustering.h"
#include "analyser/typeDefs.h"
#include "base/sysFunc.h"
#include "main/dissolve.h"

// Main processing:
//
//      - Neighbour search and distance calculations performed, valid intermolecular bonds stored defined by sites and cutoffs
//      - Separate handling of hydrogen bonds (loose and strict definitions i.e. with directionality)
//      - Perform graph based analysis to identify clusters
//      - Utilise graph networks to compute analysis metrics:
//          - Coordination numbers (for clustered moelcules)
//          - Cluster size distributions
//          - Fractal dimension (via radius of gyration??)
//          - Other potential graph network metrics: clustering coefficient, betweenness centrality?
//
//      - setUp function at beginning of each iteration:
//          - Clean up previous clustering datastructures (or save maybe - Clustering over a variety of iteration steps likely necessary to determine statistical significance of findings...)
//
//      - Output functions - option to save data?
//      - Figure out how to organise data for visualisation?


// Setup (runs when the simulation starts, not each iteration)
// bool ClusteringModule::setUp()
//{
    // Prepare module for simulation
    //  - Initialise selected sites and cut-off distances
    // std::vector<BondInfo> selectedBonds; // This needs to be linked to GUI - MUST be initialised
//}

Module::ExecutionResult ClusteringModule::process(ModuleContext &moduleContext)
{
    // Set configuration (TESTING - should be done in GUI)
    // Get the first configuration from coreData
    auto &configs = moduleContext.dissolve().coreData().configurations();
    if (configs.empty())
    {
        Messenger::error("No configurations found in simulation.\n");
        return ExecutionResult::Failed;
    }

    // Set the first configuration as our target
    targetConfiguration_ = configs.front().get();
    
    if (!targetConfiguration_)
    {
        Messenger::error("Failed to set target configuration.\n");
        return ExecutionResult::Failed;
    }

    Messenger::print("Successfully set target configuration to '{}'.\n", targetConfiguration_->name());

    // TESTING CODE - FIND TWO SITES IN THE CONFIGURATION AND SET THOSE TO SELECTED SITES
    // SelectedSites is a vector of vectors of const SpeciesSite pointers. This will be managed by GUI later, contructing BondInfo objects directly
    std::vector<const SpeciesSite *> selectedSites;
    const auto& molecules = targetConfiguration_->molecules();
    
    // Keep track of species we've already processed
    std::set<const Species*> processedSpecies;
    
    // Iterate through molecules to find unique species and their sites
    for (const auto& mol : molecules)
    {
        const auto* species = mol->species();
        
        // Skip if we've already processed this species
        if (processedSpecies.find(species) != processedSpecies.end())
            continue;
            
        // Add this species to processed set
        processedSpecies.insert(species);
        
        // Get sites from this species
        const auto& speciesSites = species->sites();  // Added const
        for (const auto& site : speciesSites)         // Changed to const reference
        {
            selectedSites.emplace_back(site.get()); 
        }
    }

    // Check if we have enough sites
    if (selectedSites.size() < 2)
    {
        Messenger::error("Not enough sites found. Expected at least 2, found {}\n", selectedSites.size());
        return ExecutionResult::Failed;
    }

    // At this point, all weve done is generate a vector of const SpeciesSite pointers. Only need two for testing. These will be selected in the GUI later.
    // Now we need to specify a bond between the two sites. 

    //Generate bonds between the first two sites - again, just for testing purposes. 
    BondInfo interBond(selectedSites[0], selectedSites[0], 2.5);
    // BondInfo interBond2(selectedSites[0], selectedSites[2], 2);
    // BondInfo interBond3(selectedSites[2], selectedSites[2], 2);
    selectedBonds_ = {interBond}; //, interBond2, interBond3}; // Will for loop this later to add each selected bond, just need one for now

    // MODULE CODE
    auto& moduleData = moduleContext.dissolve().processingModuleData();
    LineParser parser;
    std::string filename = std::format("{}.neighbourdata.txt", targetConfiguration_->niceName());
    
    // Output file for diagnostics
    if (!parser.openOutput(filename))
    {
        Messenger::error("Failed to open output file '{}'.\n", filename);
        return ExecutionResult::Failed;
    }

    // Write header with site information
    parser.writeLineF("# Analysis for sites: {} - {}\n", 
                     selectedBonds_[0].a_->parent()->name(),
                     selectedBonds_[0].b_->parent()->name()   );
                     //selectedBonds_[1].a_->parent()->name(),
                     //selectedBonds_[1].b_->parent()->name());

    
    for (const auto& [site, neighbours] : getNeighbourMap())
    {
        parser.writeLineF("\nSite '{}', at coordinates ({:.3f}, {:.3f}, {:.3f}) : {} neighbours\n\n", 
                    site->parent()->name(), 
                    site->origin().x, site->origin().y, site->origin().z,
                    neighbours.size());
        
        for (const auto& [neighbour, index] : neighbours)
        {
            parser.writeLineF("  Neighbour '{}', uniqueSiteIndex '{}' at coordinates ({:.3f}, {:.3f}, {:.3f}) : distance = {}\n", 
                            neighbour->parent()->name(), 
                            index,
                            neighbour->origin().x, neighbour->origin().y, neighbour->origin().z,
                            targetConfiguration_->box()->minimumDistance(site->origin(), neighbour->origin()));
        }
    }

    // After clustersOut is generated, print diagnostics:
    parser.writeLineF("\n=== Cluster Diagnostics ===\n");
    for (const auto& [clusterId, clusterSites] : getClusterMap())
    {
        parser.writeLineF("Cluster {}: {} site(s)\n", clusterId, clusterSites.size());
        for (const auto& site : clusterSites)
        {
            parser.writeLineF("  Site '{}', coords: ({:.3f}, {:.3f}, {:.3f})\n", 
                site->parent()->name(), site->origin().x, site->origin().y, site->origin().z);
        }
        parser.writeLineF("\n");
    }

    // Write clusters and other diagnostics...
    // Now compute the cluster size distribution and output it
    parser.writeLineF("\n=== Cluster Size Distribution ===\n");
    for (const auto& [clusterSize, count] : getSizeDistribution())
    {
        parser.writeLineF("  Cluster Size {} : {} clusters\n", clusterSize, count);
    }

    // Compute mass distribution diagnostic
    parser.writeLineF("\n=== Cluster Mass Distribution ===\n");
    for (const auto& [clusterMass, clusterVec] : getMassDistribution())
    {
        for (const auto& cluster : clusterVec)
        {
                parser.writeLineF("  Cluster Mass {:.3f} : cluster ID: {}\n", clusterMass, cluster);
        }

    }

    // Write coordination number diagnostics to output file
    parser.writeLineF("\n=== Coordination Numbers ===\n");
    for (const auto& [siteA, innerMap] : getClusterSpeciesCoordNo())
    {
        for (const auto& [siteB, coord] : innerMap)
        {
            parser.writeLineF("From site {} to site {}: {:.3f}\n",
                                siteA->parent()->name(),
                                siteB->parent()->name(),
                                coord);
        }
    }

    // Write radius of gyration diagnostics to output file
    parser.writeLineF("\n=== Radius of Gyration ===\n");
    for (const auto& [clusterID, radius] : getRadiusOfGyration())
    {
        parser.writeLineF("Cluster {}: {:.3f}\n", clusterID, radius);
    }

    return ExecutionResult::Success;
}

