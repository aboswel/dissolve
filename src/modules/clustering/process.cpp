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
    // Now we need to generate a bond between the two sites. 

    //Generate a bond between the first two sites - again, just for testing purposes. 
    BondInfo interBond(selectedSites[0], selectedSites[1], 4.5);
    std::vector<BondInfo> selectedBonds{interBond}; // Will for loop this later to add each selected bond, just need one for now

    // MODULE CODE
    auto& moduleData = moduleContext.dissolve().processingModuleData();
    LineParser parser;
    std::string filename = std::format("{}.neighbourdata.txt", targetConfiguration_->niceName());
    
    for (const auto& bonding : selectedBonds)
    {
        // Open output file
        if (!parser.openOutput(filename))
        {
            Messenger::error("Failed to open output file '{}'.\n", filename);
            return ExecutionResult::Failed;
        }

        // With revised BondInfo struct, we can more easily generate site instances from species sites (the role of SiteSelector)
        SiteSelector selectionA(targetConfiguration_, std::vector<const SpeciesSite*>{bonding.target_});
        const Analyser::SiteVector& siteVectorA = selectionA.sites();
        for (const auto& [site, index] : siteVectorA)
        {
            parser.writeLineF("Filter index (SiteAVector): {}\n", index);
        }

        // Do the same for the neighbour sites
        SiteSelector selectionB(targetConfiguration_, std::vector<const SpeciesSite*>{bonding.nbr_});
        const Analyser::SiteVector& siteVectorB = selectionB.sites();
        for (const auto& [site, index] : siteVectorB)
        {
            parser.writeLineF("Filter index (SiteBVector): {}\n", index);
        }

        // Site filter contains the target configuration, and the sites to filter by as SiteVectors
        SiteFilter filter(targetConfiguration_, siteVectorA);
        auto [filteredASites, neighbourMap] = filter.filterBySiteProximity(siteVectorB, Range(0, bonding.cutOff), 1, 100);

        // Write header with site information
        parser.writeLineF("# Analysis for sites: {} - {}\n", 
                         bonding.target_->parent()->name(),
                         bonding.nbr_->parent()->name());
        
        // Write the number of filtered sites
        parser.writeLineF("# Number of filtered sites: {}\n", filteredASites.size());
        // Write data for each filtered site and its neighbours
        for (const auto& [site, index] : filteredASites)
        {
            parser.writeLineF("Filter index: {}, coords: {:.3f}\n", index, site->origin().x);
        }

        int i{0};
        for (const auto& [site, neighbours] : neighbourMap)
        {
            parser.writeLineF("Site '{}', uniqueSiteIndex '{}' at coordinates ({:.3f}, {:.3f}, {:.3f}) : {} neighbours\n\n", 
                        site->parent()->name(), 
                        std::get<1>(filteredASites[i]),
                        site->origin().x, site->origin().y, site->origin().z,
                        neighbours.size());
            
            for (const auto& [neighbour, index] : neighbours)
            {
                parser.writeLineF("  Neighbour '{}', uniqueSiteIndex '{}' at coordinates ({:.3f}, {:.3f}, {:.3f}) : distance = {}\n", 
                                neighbour->parent()->name(), 
                                neighbour->uniqueSiteIndex().value_or(-1),
                                neighbour->origin().x, neighbour->origin().y, neighbour->origin().z,
                                targetConfiguration_->box()->minimumDistance(site->origin(), neighbour->origin()));
            }
            i++;
        }
    }

    return ExecutionResult::Success;
}