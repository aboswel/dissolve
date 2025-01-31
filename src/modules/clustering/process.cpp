// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2024 Team Dissolve and contributors

// process.cpp
#include "base/sysFunc.h"
#include "main/dissolve.h"
#include "analyser/siteFilter.h"
#include "analyser/siteSelector.h"
#include "modules/clustering/clustering.h"

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

    // TESTING CODE
    std::vector<std::reference_wrapper<const std::unique_ptr<SpeciesSite>>> selectedSites;  // Changed to const
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
            selectedSites.emplace_back(std::cref(site));  // Changed to cref
        }
    }

    // Check if we have enough sites
    if (selectedSites.size() < 2)
    {
        Messenger::error("Not enough sites found. Expected at least 2, found {}\n", selectedSites.size());
        return ExecutionResult::Failed;
    }

    BondInfo interBond = {selectedSites[0], selectedSites[1], 4.5};
    std::vector<BondInfo> selectedBonds{interBond};

    // MODULE CODE
    auto& moduleData = moduleContext.dissolve().processingModuleData();
    LineParser parser;
    std::string filename = fmt::format("{}.neighbourdata.txt", targetConfiguration_->niceName());
    
    for (const auto& bonding : selectedBonds)
    {
        std::vector<const SpeciesSite*> sitesA;
        sitesA.push_back(bonding.site1_.get().get());

        std::vector<const SpeciesSite*> sitesB;
        sitesB.push_back(bonding.site2_.get().get());

        SiteSelector selectionA(targetConfiguration_, sitesA);
        const Analyser::SiteVector& siteVectorA = selectionA.sites();

        SiteSelector selectionB(targetConfiguration_, sitesB);
        const Analyser::SiteVector& siteVectorB = selectionB.sites();

        SiteFilter filter(targetConfiguration_, siteVectorA);
        auto [filteredASites, neighbourMap] = filter.filterBySiteProximity(siteVectorB, Range(0, bonding.cutOff), 1, 10);

        // Open output file
        if (!parser.openOutput(filename))
        {
            Messenger::error("Failed to open output file '{}'.\n", filename);
            return ExecutionResult::Failed;
        }

        // Write header with site information
        parser.writeLineF("# Analysis for sites: {} - {}\n", 
                         bonding.site1_.get().get()->parent()->name(),
                         bonding.site2_.get().get()->parent()->name());
        
        // Write the number of filtered sites
        parser.writeLineF("# Number of filtered sites: {}\n", filteredASites.size());

        // Write data for each filtered site and its neighbours
        for (const auto& [site, neighbours] : neighbourMap)
        {
            parser.writeLineF("Site '{}': {} neighbours\n", 
                        site->parent()->name(), 
                        neighbours.size());
            
            for (const auto& [neighbour, distance] : neighbours)
            {
                parser.writeLineF("  Neighbour '{}': distance = {}\n", 
                                neighbour->parent()->name(), 
                                distance);
            }
        }
    }

    return ExecutionResult::Success;
}