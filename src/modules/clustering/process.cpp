// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2024 Team Dissolve and contributors

// process.cpp
#include "base/sysFunc.h"
#include "main/dissolve.h"
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
    // Now we need to specify a bond between the two sites. 

    //Generate bonds between the first two sites - again, just for testing purposes. 
    BondInfo interBond(selectedSites[0], selectedSites[1], 4.5);
    BondInfo interBond2(selectedSites[1], selectedSites[2], 4.5);
    std::vector<BondInfo> selectedBonds{interBond, interBond2}; // Will for loop this later to add each selected bond, just need one for now

    // MODULE CODE
    auto& moduleData = moduleContext.dissolve().processingModuleData();
    LineParser parser;
    std::string filename = std::format("{}.neighbourdata.txt", targetConfiguration_->niceName());
    
    // Generate symmetric adjacency list for the selected bond (index 0 for now)
    auto [filteredSites, neighbourMap] = SiteFiltering(targetConfiguration_, selectedBonds);

    // Output file for diagnostics
    if (!parser.openOutput(filename))
    {
        Messenger::error("Failed to open output file '{}'.\n", filename);
        return ExecutionResult::Failed;
    }

    // Write header with site information
    parser.writeLineF("# Analysis for sites: {} - {}, {} - {}\n", 
                     selectedBonds[0].a_->parent()->name(),
                     selectedBonds[0].b_->parent()->name(),
                     selectedBonds[1].a_->parent()->name(),
                     selectedBonds[1].b_->parent()->name());

    for (const auto& [site, index] : filteredSites)
    {
        parser.writeLineF("Filter name: {}, index: {}, coords: {:.3f}\n", site->parent()->name(), index, site->origin().x);
    }

    int i{0};
    for (const auto& [site, neighbours] : neighbourMap)
    {
        parser.writeLineF("Site '{}', uniqueSiteIndex '{}' at coordinates ({:.3f}, {:.3f}, {:.3f}) : {} neighbours\n\n", 
                    site->parent()->name(), 
                    std::get<1>(filteredSites[i]),
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
        i++;
    }

    /*  Pre seperated function for filtering
    for (const auto& bonding : selectedBonds)
    {
        // Open output file
        if (!parser.openOutput(filename))
        {
            Messenger::error("Failed to open output file '{}'.\n", filename);
            return ExecutionResult::Failed;
        }

        // With revised BondInfo struct, we can more easily generate site instances from species sites (the role of SiteSelector)
        SiteSelector selectionA(targetConfiguration_, std::vector<const SpeciesSite*>{bonding.a_});
        const Analyser::SiteVector& siteVectorA = selectionA.sites();
        for (const auto& [site, index] : siteVectorA)
        {
            parser.writeLineF("Filter index (SiteAVector): {}\n", index);
        }

        // Do the same for the neighbour sites
        SiteSelector selectionB(targetConfiguration_, std::vector<const SpeciesSite*>{bonding.b_});
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
                         bonding.a_->parent()->name(),
                         bonding.b_->parent()->name());
        
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
    */
    return ExecutionResult::Success;
}

// Lets move the site proximity filtering to a separate function, will to deal with all bond vectors at once. How to handle strict H bonding??
// This needs to generate a symmetric adjacency list for a given intermolecular bond - keep this in the format of filterBySiteProximity return?
std::pair<Analyser::SiteVector, Analyser::SiteMap> ClusteringModule::SiteFiltering(Configuration *cfg_, std::vector<BondInfo> bonds)
{
    Analyser::SiteVector combinedFilteredSites;
    Analyser::SiteMap combinedNeighbourMap;

    for (auto bond : bonds)
    {
        // Transform SiteObjects to instances ready for the filter function
        SiteSelector selectionA(cfg_, std::vector<const SpeciesSite*>{bond.a_});
        const Analyser::SiteVector& siteVectorA = selectionA.sites();

        SiteSelector selectionB(cfg_, std::vector<const SpeciesSite*>{bond.b_});
        const Analyser::SiteVector& siteVectorB = selectionB.sites();

        SiteFilter filterA(targetConfiguration_, siteVectorA);
        auto [filteredASites, neighbourMapA] = filterA.filterBySiteProximity(siteVectorB, Range(0.01, bond.cutOff), 1, 100); // min of 0.01 to avoid self-bonding

        // If we're dealing with a bond between the same site, it's already symmetric. Running twice will add dupes
        if (bond.a_ == bond.b_)
        {
            combinedFilteredSites.reserve(combinedFilteredSites.size() + filteredASites.size());
            combinedFilteredSites.insert(combinedFilteredSites.end(), filteredASites.begin(), filteredASites.end());
            combinedNeighbourMap.insert(neighbourMapA.begin(), neighbourMapA.end());
            continue;
        }

        // Doing this twice allows us to generate a symmetric adjacency list
        SiteFilter filterB(targetConfiguration_, siteVectorB);
        auto [filteredBSites, neighbourMapB] = filterB.filterBySiteProximity(siteVectorA, Range(0.01, bond.cutOff), 1, 100);

        // Combining the filtered sites into a single vector
        combinedFilteredSites.reserve(combinedFilteredSites.size() + filteredASites.size() + filteredBSites.size());
        combinedFilteredSites.insert(combinedFilteredSites.end(), filteredASites.begin(), filteredASites.end());
        combinedFilteredSites.insert(combinedFilteredSites.end(), filteredBSites.begin(), filteredBSites.end()); 

        // Combining the neighbour maps into a single map. Because keys may already exist, need to check for them and add neighbours if exists.
        for (auto neighbourMap : {neighbourMapA, neighbourMapB})
        {
            for (const auto& [site, neighbours] : neighbourMap)
            {
                if (combinedNeighbourMap.contains(site))
                {
                    combinedNeighbourMap[site].insert(combinedNeighbourMap[site].end(), neighbours.begin(), neighbours.end());
                }
                else
                {
                    combinedNeighbourMap.insert({site, neighbours});
                }
            }
        }
    }
    return {combinedFilteredSites, combinedNeighbourMap};
}
