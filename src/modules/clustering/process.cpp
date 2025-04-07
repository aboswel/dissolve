// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

// process.cpp
#include "modules/clustering/clustering.h"
#include "analyser/typeDefs.h"
#include "base/sysFunc.h"
#include "generator/box.h"
#include "math/regression.h"

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


bool ClusteringModule::setUp(ModuleContext &moduleContext, Flags<KeywordBase::KeywordSignal> actionSignals)
{
    // Existing oldConfigurations in the restart file generate a new tab which is problematic. Only want this to run we've freshly loaded a save
    if (!clusterConfig_)
    {
        for (const auto &config : moduleContext.dissolve().coreData().configurations())
        {
            Messenger::print("\nCurrent Configs : {}", config->name());
            functionTrack(std::format("Iteration{}{}.txt", moduleContext.dissolve().iteration(), config->name()), 1); 
        }
        auto *oldConfig = moduleContext.dissolve().coreData().findConfiguration("clusters");
        if (oldConfig)
        {
            moduleContext.dissolve().coreData().removeConfiguration(oldConfig);
            Messenger::print("Removed existing cluster configuration");
        }
    }
    // Clear all existing data
    selectedBonds_.clear();
    neighbourMap_.clear();
    clusterMap_.clear();
    sizeDistribution_.clear();
    massDistribution_.clear();
    clusterMasses_.clear();
    return true;
}

Module::ExecutionResult ClusteringModule::process(ModuleContext &moduleContext)
{
    /*
    // Set configuration (TESTING - should be done in GUI)
    // Get the first configuration from coreData
    auto &configs = moduleContext.dissolve().coreData().configurations();
    if (configs.empty())
    {
        Messenger::error("No configurations found in simulation.\n");
        return ExecutionResult::Failed;
    }

    if (!targetConfiguration_)
    {
        Messenger::error("Failed to set target configuration.\n");
        return ExecutionResult::Failed;
    }

    Messenger::print("Successfully set target configuration to '{}'.\n", targetConfiguration_->name());

    // Attempt to create and read an inputFile
    if (!inputCreated)
    {
        produceInput(selectedSites);
        inputCreated = true;
    }

    if (!readInput(selectedSites))
    {
        Messenger::error("Failed to retrieve inputs from clusteringSettings.txt - Please select sites!");
        return ExecutionResult::Failed;
    }
    */

    // Check if we have enough sites
    //if (selectedSites.size() < 2)
    //{
    //    Messenger::error("Not enough sites found. Expected at least 2, found {}\n", selectedSites.size());
    //    return ExecutionResult::Failed;
    //}

    // At this point, all weve done is generate a vector of const SpeciesSite pointers. Only need two for testing. These will be selected in the GUI later.
    // Now we need to specify a bond between the two sites. 

    //Generate bonds between the first two sites - again, just for testing purposes. 
    //BondInfo interBond(selectedSites[0], selectedSites[0], 4.4);
    // BondInfo interBond2(selectedSites[0], selectedSites[2], 2);
    // BondInfo interBond3(selectedSites[2], selectedSites[2], 2);
    //selectedBonds_ = {interBond}; //, interBond2, interBond3}; // Will for loop this later to add each selected bond, just need one for now

    // MODULE CODE
    clustersConfigGenerated_ = false;
    auto& moduleData = moduleContext.dissolve().processingModuleData();
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
    parser.writeLineF("# Analysis for sites: {} - {}\n", 
                     selectedBonds_[0].a_->parent()->name(),
                     selectedBonds_[0].b_->parent()->name()   );
                     //selectedBonds_[1].a_->parent()->name(),
                     //selectedBonds_[1].b_->parent()->name());


    // Process code

    // Produce NeighbourMap
    neighbourMap_.clear();
    for (auto bond : selectedBonds_) // Change to getter function when this is actually selected in GUI
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

    // ClusterMap generation 1.0
    /*
    // We need to traverse the cluster:
    //      We need a starting site (sites) - This is the first site in the cluster (iterate down neighbourmap for this)
    //      From this site, we list all the neighbours
    //      We purge neighbours that exist in the visited set
    //      Those that are left, are put in the visited set and added to the cluster
    //      Once added, these neighbours now become sites
    //      All neighbours of these sites are now found, purged, added and become next sites
    //      Continue until there are no more valid neighbours. Then we iterate down the siteMap until we find the next unvisited site, repeat process.
    //      Ought to see if theres a better way to do this... Recursive implementation?

    std::unordered_set<const Site*> visited;
    std::vector<const Site*> sites;
    int clusterTrack{1};

    for (const auto &[clusterStart, _] : neighbourMap_)
    {
        // check if the selected site has been visited.
        if (visited.contains(clusterStart))
        {
            continue;
        }

        // cluster start becomes sites
        sites.emplace_back(clusterStart);

        // While loop continues for this cluster until sites ends up empty, effectively until there are no more unvisited
        // neighbours
        while (!sites.empty())
        {
            // Add sites to visited and cluster
            for (const auto &s : sites)
            {
                visited.insert(s);
            }
            clusterMap_[clusterTrack].insert(clusterMap_[clusterTrack].end(), sites.begin(), sites.end());

            // Use an unordered_set to collect new neighbours (avoids duplicates)
            std::unordered_set<const Site *> newNeighbours;
            for (const auto &si : sites)
            {
                for (const auto &nbr : neighbourMap_[si])
                {
                    if (visited.contains(std::get<0>(nbr)))
                    {
                        continue;
                    }
                    newNeighbours.insert(std::get<0>(nbr));
                }
            }
            // Convert the set of newNeighbours to a vector and use as next sites
            sites.assign(newNeighbours.begin(), newNeighbours.end());
        }
        clusterTrack++;
    }
    */

    // ClusterMap generation 2.0
    {
        clusterMap_.clear();
        std::unordered_set<const Site*> processed;
        int clusterTrack{1};
        for (const auto &[clusterStart, _] : neighbourMap_)
        {
            if (!processed.contains(clusterStart)) 
            {
                std::unordered_set<const Site*> visited{clusterStart};
                buildCluster(clusterStart, visited);
                processed.insert(visited.begin(), visited.end());
                clusterMap_.insert({clusterTrack, std::vector<const Site*>(visited.begin(), visited.end())});
                clusterTrack++;
            }
        }
    }

    // Cluster size distribution
    sizeDistribution_.clear();
    for (const auto& [_, members] : clusterMap_)
    {
        sizeDistribution_[members.size()]++;
    }

    // Cluster mass calculation
    clusterMasses_.clear();
    for (const auto& [clusterID, memberVec] : clusterMap_)
    {
        float clusterMass{0};
        for (const auto& member : memberVec)
        {
            clusterMass += member->parent()->parent()->mass();
        }
        clusterMasses_[clusterID] = clusterMass;
    }
    
    // Cluster mass distribution
    massDistribution_.clear();
    for (const auto& [clusterID, clusterMass] : clusterMasses_)
    {
        massDistribution_[clusterMass].emplace_back(clusterID);
    }

    radiusOfGyration_.clear();
    generateRadiusOfGyration();
    generateFractalDimension();
    viewingReady = true;


    // Diagnostics
    /*
    for (const auto& [site, neighbours] : neighbourMap_)
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
    for (const auto& [clusterId, clusterSites] : clusterMap_)
    {
        parser.writeLineF("Cluster {}: {} site(s)\n", clusterId, clusterSites.size());
        for (const auto& site : clusterSites)
        {
            parser.writeLineF("  Site '{}', coords: ({:.3f}, {:.3f}, {:.3f})\n", 
                site->parent()->name(), site->origin().x, site->origin().y, site->origin().z);
        }
        parser.writeLineF("\n");
    }
    */

    // Write clusters and other diagnostics...
    // Now compute the cluster size distribution and output it
    parser.writeLineF("\n=== Cluster Size Distribution ===\n");
    for (const auto& [clusterSize, count] : sizeDistribution_)
    {
        parser.writeLineF("  Cluster Size {} : {} clusters\n", clusterSize, count);
    }

    
    // Compute mass distribution diagnostic
    parser.writeLineF("\n=== Cluster Mass Distribution ===\n");
    for (const auto& [clusterMass, clusterVec] : massDistribution_)
    {
        for (const auto& cluster : clusterVec)
        {
                parser.writeLineF("  Cluster Mass {:.3f} : cluster ID: {}\n", clusterMass, cluster);
        }

    }
    
    /*
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
    */

    // Write radius of gyration diagnostics to output file
    parser.writeLineF("\n=== Radius of Gyration ===\n");
    for (const auto& [clusterID, radius] : radiusOfGyration_)
    {
        parser.writeLineF("Cluster {}: {:.3f}\n", clusterID, radius);
    }

    // Write fractal dimension diagnostics to output file
    parser.writeLineF("\n=== Fractal Dimension ===\n");
    parser.writeLineF("Fractal dimension: {:.3f}\n", fractalDimension_);

    // Write a formatted output of cluster mass : radius of gyration
    parser.writeLineF("\n=== Rg : Mass ===\n");
    for (const auto &[clusterID, Rg] : radiusOfGyration_)
    {
        parser.writeLineF("{}\n", clusterMasses_[clusterID]);
    }
    parser.writeLineF("\n");
    for (const auto &[clusterID, Rg] : radiusOfGyration_)
    {
        parser.writeLineF("{}\n", Rg);
    }

    return ExecutionResult::Success;
}

void ClusteringModule::buildCluster(const Site* startSite, std::unordered_set<const Site*>& visited)
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

// Everything from here on down is temporary code allowing runtime setting adjustment and 3D viewing by hijacking a configuration to display clusters
// Shouldn't be necessary from here on out!
/*
void ClusteringModule::produceInput(std::vector<const SpeciesSite*> sites)
{
    LineParser parser;

    // Check if file exists. Skip if so.
    if (DissolveSys::fileExists("clusteringSettings.txt"))
        return;
    
    // Open file for writing
    if (!parser.openOutput("clusteringSettings.txt", true))
        Messenger::error("Failed to open clustering settings.\n");
    
    // Write header/banner
    parser.writeBannerComment("# Clustering Settings");
    parser.writeLineF("\n# Avaliable sites : ");

    for (const auto& site : sites)
        parser.writeLineF("{} ", site->parent()->name());

    // Write settings sections with values
    // Section 1: General Settings
    parser.writeLineF("\n\n# Intermolecular Bond\n");
    parser.writeLineF("Site1: \n");
    parser.writeLineF("Site2: \n");
    parser.writeLineF("Cutoff: \n\n");
    parser.writeLineF("Min Size For Fractal Dimension: \n\n");
}

bool ClusteringModule::readInput(std::vector<const SpeciesSite*> sites)
{
    LineParser parser;
    
    // Open file for reading
    if (!parser.openInput("clusteringSettings.txt"))
    {
        Messenger::error("Failed to open clustering settings file.\n");
        return false;
    }
    
    // Clear existing bonds
    selectedBonds_.clear();
    
    std::string site1Name, site2Name;
    int minFrac;
    double cutoff = 0.0;
    bool inBondSection = false;
    
    Messenger::print("Reading clustering settings file...\n");
    
    // Read through file line by line
    while (!parser.eofOrBlank())
    {
        if (parser.getArgsDelim() != LineParser::Success)
            continue;
        
        // Skip empty lines or comment-only lines
        if (parser.nArgs() == 0 || parser.argsv(0)[0] == '#')
            continue;
        
        // Check if we're starting a new bond with Site1
        if (parser.argsv(0) == std::string("Site1:"))
        {
            // If we already have a complete previous bond definition, store it
            if (inBondSection && !site1Name.empty() && !site2Name.empty() && cutoff > 0.0)
            {
                if (addBondFromNames(sites, site1Name, site2Name, cutoff))
                {
                    Messenger::print("Added previous bond before starting new one\n");
                }
            }
            
            // Start a new bond section
            inBondSection = true;
            site1Name.clear();
            site2Name.clear();
            cutoff = 0.0;
            
            // Store the new Site1 value
            if (parser.nArgs() > 1) {
                site1Name = parser.argsv(1);
                Messenger::print("  Found new bond starting with Site1: {}\n", site1Name);
            }
            else
                Messenger::print("  Warning: Site1 has no value\n");
        }
        else if (parser.argsv(0) == std::string("Site2:"))
        {
            if (!inBondSection) {
                Messenger::print("  Warning: Found Site2 without preceding Site1\n");
                inBondSection = true;
            }
            
            if (parser.nArgs() > 1) {
                site2Name = parser.argsv(1);
                Messenger::print("  Found Site2: {}\n", site2Name);
            }
            else
                Messenger::print("  Warning: Site2 has no value\n");
        }
        else if (parser.argsv(0) == std::string("Cutoff:"))
        {
            if (!inBondSection) {
                Messenger::print("  Warning: Found Cutoff without preceding Site1/Site2\n");
                inBondSection = true;
            }
            
            if (parser.nArgs() > 1) {
                cutoff = parser.argd(1);
                Messenger::print("  Found Cutoff: {}\n", cutoff);
            }
            else
                Messenger::print("  Warning: Cutoff has no value\n");
        }
        else if (parser.argsv(0) == std::string("MF:"))
        {   
            if (parser.nArgs() > 1) {
                minFrac = parser.argd(1);
                gyrationMinSize_ = minFrac;
                Messenger::print("  Found Min size for fractal calc: {}\n", minFrac);
            }
            else
                Messenger::print("  Warning: Min Fractal cluster size not set. Defaulting to 3.\n");
        }
    }
    
    // Don't forget to process the last bond
    if (inBondSection && !site1Name.empty() && !site2Name.empty() && cutoff > 0.0)
    {
        addBondFromNames(sites, site1Name, site2Name, cutoff);
    }
    
    // Check if we have at least one bond
    if (selectedBonds_.empty())
    {
        Messenger::error("No valid bond definitions found in clustering settings file.\n");
        
        // Show available site names
        Messenger::print("Available site names are: ");
        for (const auto& site : sites)
            Messenger::print("{} ", site->parent()->name());
        Messenger::print("\n");
        
        return false;
    }
    
    Messenger::print("Successfully loaded {} clustering bond definitions.\n", selectedBonds_.size());
    
    return true;
}

// Helper function to find sites and create bond
bool ClusteringModule::addBondFromNames(const std::vector<const SpeciesSite*>& sites, 
                                      const std::string& site1Name, 
                                      const std::string& site2Name,
                                      double cutoff)
{
    // Find the corresponding SpeciesSite pointers
    const SpeciesSite* site1Ptr = nullptr;
    const SpeciesSite* site2Ptr = nullptr;
    
    for (const auto& site : sites)
    {
        if (site->parent()->name() == site1Name)
            site1Ptr = site;
        if (site->parent()->name() == site2Name)
            site2Ptr = site;
    }
    
    // Check if both sites were found
    if (!site1Ptr || !site2Ptr)
    {
        Messenger::error("One or both specified sites not found: '{}', '{}'\n", site1Name, site2Name);
        return false;
    }
    
    // Create BondInfo object and add to selectedBonds_
    selectedBonds_.emplace_back(site1Ptr, site2Ptr, cutoff);
    
    Messenger::print("Added clustering bond: {} - {} (cutoff: {:.3f})\n", 
                   site1Name, site2Name, cutoff);
    
    return true;
}
*/
/*
Configuration* ClusteringModule::generateClustersConfig(Dissolve *dissolve, Configuration *source)
{
    Configuration* clusterConfig{nullptr};
    // As this gets called when the GUI is set-up (does it need to be??), we need to make sure we don't try to access anything that hasn't been created yet.
    if (!viewingReady)
        return nullptr;
    
    for (const auto &configs : dissolve->coreData().configurations())
        Messenger::print(configs->name());
    
    // Check if this configuration already exists, if so delete and regenerate
    auto *oldConfig = dissolve->coreData().findConfiguration("clusters");
    if (oldConfig)
    {
        dissolve->coreData().removeConfiguration(oldConfig);
        oldConfig = nullptr;
        Messenger::print("Removed old configuration, regenerating...");
    }
    
    // Set up the basic parameters for a new configuration
    auto *clusterConfig = dissolve->coreData().addConfiguration();
    clusterConfig->setName("clusters");
    clusterConfig->createBox(source->box()->axisLengths(), source->box()->axisAngles());
    
    // Add the molecules we want.
    // Get the molecules from the sites listed in the clustermap and add them to the configuration
    for (const auto &[clusterID, mems] : clusterMap_)
    {
        if (mems.size() >= minClusterForConfig && mems.size() <= maxClusterForConfig)
        {
            for (const auto &site : mems)
            {
                clusterConfig->copyMolecule(std::const_pointer_cast<Molecule>(site->molecule()));
            }
        }
    }

    if (clusterConfig)
    {
        if (clusterConfig->nAtoms() == 0)
        {
            Messenger::error("No displayable clusters were found! Tweak min/max values.");
            dissolve->coreData().removeConfiguration(clusterConfig);
            clusterConfig = nullptr;
        }
        else
            Messenger::print("clusters config successfully generated");
    }
        
    return clusterConfig;
    
}
*/
Configuration *ClusteringModule::generateClustersConfig(Dissolve *dissolve, Configuration *source)
{
    Messenger::print("Using cutoff {}", selectedBonds_[0].cutOff);
    functionTrack("NewConfiguration1.neighbourdata.txt", 31);
    // If the cluster config has already been generated this iteration, just return it
    if (clustersConfigGenerated_)
    {
        Messenger::print("Returned existing cluster config");
        functionTrack("NewConfiguration1.neighbourdata.txt", 311);
        return clusterConfig_;
    }

    // If the config does not exist at all, make it
    functionTrack("NewConfiguration1.neighbourdata.txt", 32);
    if (!clusterConfig_)
    {
        functionTrack("NewConfiguration1.neighbourdata.txt", 321);
        clusterConfig_ = dissolve->coreData().addConfiguration();
        functionTrack("NewConfiguration1.neighbourdata.txt", 322);
        clusterConfig_->setName("clusters");
        functionTrack("NewConfiguration1.neighbourdata.txt", 323);
        clusterConfig_->createBox(source->box()->axisLengths(), source->box()->axisAngles());
        functionTrack("NewConfiguration1.neighbourdata.txt", 324);
        Messenger::print("generated new cluster config");

    }
    else
    {
        // If the configuration exists, but needs updating, empty it
        functionTrack("NewConfiguration1.neighbourdata.txt", 325);
        clusterConfig_->empty();
        functionTrack("NewConfiguration1.neighbourdata.txt", 326);
        Messenger::print("emptied cluster config");
    }

    // Add the molecules we want.
    // Get the molecules from the sites listed in the clustermap and add them to the configuration
    
    for (const auto &[clusterID, mems] : clusterMap_)
    {
        if (mems.size() >= minClusterForConfig && mems.size() <= maxClusterForConfig)
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

    functionTrack("NewConfiguration1.neighbourdata.txt", 33);
    if (clusterConfig_->nAtoms() == 0)
    {
        Messenger::error("No displayable clusters were found! Tweak min/max values.");
        functionTrack("NewConfiguration1.neighbourdata.txt", 331);
    }
    else
    {
        Messenger::print("clusters config successfully generated");
        functionTrack("NewConfiguration1.neighbourdata.txt", 332);
    }
    clustersConfigGenerated_ = true;
    functionTrack("NewConfiguration1.neighbourdata.txt", 333);
    return clusterConfig_;
}

Configuration* ClusteringModule::getSourceConfig()
{
    return targetConfiguration_;
}


// Calculate the radius of gyration of cluster sizes above the min value
void ClusteringModule::generateRadiusOfGyration()
{
    const Box *box = targetConfiguration_->box();
    // Iterate through cluster map, skip clusters below size min
    for (const auto& [clusterID, clusterVec] : clusterMap_)
    {
        if (clusterVec.size() < gyrationMinSize_)
        {
            continue;
        }
        // Now calculate the centre of mass of the given cluster with regards to a reference site (first in clustermap)
        // Collect the coordinates of each member, multiply by mass of parent, accumlate a total, then divide by the mass of the cluster
        Vec3<double> massWeightedTotalVec{0,0,0};
        const Site* refSite{clusterVec[0]}; // Define a reference site (the first member of the cluster)
        for (const auto& clusterMem : clusterVec)
        {
            // Accumlate mass weighted total vector from reference site
            massWeightedTotalVec += (box->minimumVector(refSite->origin(), clusterMem->origin())) * clusterMem->parent()->parent()->mass();
        }
        massWeightedTotalVec /= clusterMasses_[clusterID];
        clusterCoM_[clusterID] = massWeightedTotalVec;
        double massWeightedDistanceSqrd{0};
        // Now time to calculate the Radius of Gyration
        // Need to run through the clusterMap again, get the mim sqrd distance of site from CoM using ref as origin
        for (const auto& clusterMem : clusterVec)
        {
            massWeightedDistanceSqrd += (box->minimumDistanceSquared(box->minimumVector(clusterMem->origin(), refSite->origin()), clusterCoM_[clusterID]))
                                            * clusterMem->parent()->parent()->mass();
        }
        radiusOfGyration_[clusterID] = std::sqrt(massWeightedDistanceSqrd / clusterMasses_[clusterID]);
    }
}

void ClusteringModule::generateFractalDimension()
{
    // Create a Data1D object of the log log plot, perform linear regression, return gradient
    Data1D loglog;
    loglog.initialise(radiusOfGyration_.size(), false); // Add errors to this at some point?

    // Generate Data1D
    for (const auto& [clusterID, rg] : radiusOfGyration_)
    {
        loglog.addPoint(std::log(radiusOfGyration_[clusterID]), std::log(clusterMasses_[clusterID]));
    }
    // Perform linear regression
    fractalDimension_ = (Regression::linear(loglog, radiusOfGyration_.size()));
}

void ClusteringModule::functionTrack(std::string_view output, int track)
{
    LineParser parser;
    parser.openOutput(output);
    parser.writeLineF("{}", track);
}