// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

#include "modules/clustering/clustering.h"
#include "analyser/typeDefs.h"
#include "base/sysFunc.h"
#include "main/dissolve.h"
#include "classes/box.h"
#include <unordered_set>
#include <tuple>
#include <cmath>

// Getter functions
Analyser::SiteMap& ClusteringModule::getNeighbourMap()
{
    if (neighbourMap_.empty())
    {
        generateNeighbourMap();
    }
    return neighbourMap_;
}

std::map<int, std::vector<const Site*>>& ClusteringModule::getClusterMap()
{
    if (clusterMap_.empty())
    {
        generateClusterMap();
    }
    return clusterMap_;
}

std::map<int, int>& ClusteringModule::getSizeDistribution()
{
    if (sizeDistribution_.empty())
    {
        generateSizeDistribution();
    }
    return sizeDistribution_;
}

std::map<float, std::vector<int>>& ClusteringModule::getMassDistribution()
{
    if (massDistribution_.empty())
    {
        generateMassDistribution();
    }
    return massDistribution_;
}

std::map<const SpeciesSite*, std::map<const SpeciesSite*, float>>&  ClusteringModule::getClusterSpeciesCoordNo()
{
    if (clusterSpeciesCoordNo_.empty())
    {
        generateClusterSpeciesCoordNo();
    }
    return clusterSpeciesCoordNo_;
}


std::map<int, float>& ClusteringModule::getRadiusOfGyration()
{
    if (radiusOfGyration_.empty())
    {
        generateRadiusOfGyration();
    }
    return radiusOfGyration_;
}

// Process functions
// Lets move the site proximity filtering to a separate function, will to deal with all bond vectors at once. How to handle strict H bonding??
// This needs to generate a symmetric adjacency list for a given intermolecular bond - keep this in the format of filterBySiteProximity return?
void ClusteringModule::generateNeighbourMap()
{
    Analyser::SiteMap combinedNeighbourMap;

    for (auto bond : selectedBonds_) // Change to getter function when this is actually selected in GUI
    {
        Analyser::SiteMap neighbourMapA, neighbourMapB;

        // Transform SiteObjects to instances ready for the filter function
        SiteSelector selectionA(targetConfiguration_, std::vector<const SpeciesSite*>{bond.a_});
        const Analyser::SiteVector& siteVectorA = selectionA.sites();

        SiteSelector selectionB(targetConfiguration_, std::vector<const SpeciesSite*>{bond.b_});
        const Analyser::SiteVector& siteVectorB = selectionB.sites();

        SiteFilter filterA(targetConfiguration_, siteVectorA);
        std::tie(std::ignore, neighbourMapA) = filterA.filterBySiteProximity(siteVectorB, Range(0.001, bond.cutOff), 1, 100); // min of 0.01 to avoid self-bonding

        // If we're dealing with a bond between the same site, it's already symmetric. Running twice will add dupes
        if (bond.a_ != bond.b_)
        {
            SiteFilter filterB(targetConfiguration_, siteVectorB);
            std::tie(std::ignore, neighbourMapB) = filterB.filterBySiteProximity(siteVectorA, Range(0.001, bond.cutOff), 1, 100);
        }

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
    neighbourMap_ = combinedNeighbourMap;
}

// Organise combined neighbour map into clusters
// {Cluster No (for ID) : (sites* contributing to cluster)}
void ClusteringModule::generateClusterMap()
{
    // We need to traverse the cluster:
    //      We need a starting site (sites) - This is the first site in the cluster (iterate down neighbourmap for this)
    //      From this site, we list all the neighbours
    //      We purge neighbours that exist in the visited set
    //      Those that are left, are put in the visited set and added to the cluster
    //      Once added, these neighbours now become sites
    //      All neighbours of these sites are now found, purged, added and become next sites
    //      Continue until there are no more valid neighbours. Then we iterate down the siteMap until we find the next unvisited site, repeat process.
    //      Ought to see if theres a better way to do this...

    std::unordered_set<const Site*> visited;
    std::map<int, std::vector<const Site*>> clusters; // {ClusterID : vector of sites belonging to cluster}
    std::vector<const Site*> sites;
    int clusterTrack{1};

    for (const auto& [clusterStart, _] : getNeighbourMap())
    {
        // check if the selected site has been visited.
        if (visited.contains(clusterStart))
        {
            continue;
        }

        // cluster start becomes sites
        sites.emplace_back(clusterStart);

        // While loop continues for this cluster until sites ends up empty, effectively until there are no more unvisited neighbours
        while (!sites.empty())
        {
            // Add sites to visited and cluster
            for (const auto& s : sites)
            {
                visited.insert(s);
            }
            clusters[clusterTrack].insert(clusters[clusterTrack].end(), sites.begin(), sites.end());

            // Use an unordered_set to collect new neighbours (avoids duplicates)
            std::unordered_set<const Site*> newNeighbours;
            for (const auto& si : sites)
            {
                for (const auto& nbr : neighbourMap_[si])
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
    clusterMap_ = clusters;
}

// Basic metric computation
// Cluster size distribution
void ClusteringModule::generateSizeDistribution()
{
    std::map<int, int> distribution; // {cluster size : no of clusters}

    // Iterate through the cluster map
    for (const auto& [_, members] : getClusterMap())
    {
        distribution[members.size()]++;
    }
    sizeDistribution_ = distribution;
}
// Cluster mass distribution
void ClusteringModule::generateMassDistribution()
{
    std::map<float, std::vector<int>> distribution; // {cluster mass : vector of clusterIDs with that mass}

    // Iterate through cluster map
    for (const auto& [clusterID, memberVec] : getClusterMap())
    {
        float clusterMass{0};
        for (const auto& member : memberVec)
        {
            clusterMass += member->parent()->parent()->mass();
        }
        distribution[clusterMass].emplace_back(clusterID);
    }
    massDistribution_ = distribution;
}
// Coordination numbers for each intermolecular bond (symmetric) in cluster sizes ranging from min to max
// Min and max values optional. 0 values indicate no restriction
// {(subject site, object site) : coord no.}
// I'm not sure how useful this actually is...
void ClusteringModule::generateClusterSpeciesCoordNo()
{
    std::map<const SpeciesSite*, std::map<const SpeciesSite*, float>> bonds; // {Subject Site : {Object Site : coord No.}}
    std::map<const SpeciesSite*, int> instances;
    bool fullRange = (minCNSize_ == 0 && maxCNSize_ == 0);

    // Start iterating through the cluster map.
    for (auto const& [_, clusterVec] : getClusterMap())
    {
        // Check if size restrictions apply
        if (fullRange || (clusterVec.size() >= minCNSize_ && clusterVec.size() <= maxCNSize_))
        {
            // Iterate through each member of the cluster
            for (auto const& clusterMem : clusterVec)
            {
                // For each species found, increment the total number of that species by one.
                instances[clusterMem->parent()]++;

                // Find the member in the neighbour map, 
                for (auto const& [clusterMemNbr, index] : getNeighbourMap()[clusterMem])
                {
                    bonds[clusterMem->parent()][clusterMemNbr->parent()]++;
                }
            }
        }
    }
    // Average the coordination numbers
    for (const auto& [siteA, num] : instances)
    {
        for (const auto& [siteB, coordNo] : bonds[siteA])
        {
            bonds[siteA][siteB] /= num;
        }
    }
    clusterSpeciesCoordNo_ = bonds;
}

// Calculate the radius of gyration of cluster sizes above the min value
void ClusteringModule::generateRadiusOfGyration()
{
    const Box *box = targetConfiguration_->box();
    // Iterate through cluster map, skip clusters below size min
    for (const auto& [clusterID, clusterVec] : getClusterMap())
    {
        if (clusterVec.size() < gyrationMinSize_)
        {
            continue;
        }
        // Now calculate the centre of mass with regards to the origin of the configuration
        // Collect the coordinates of each member, multiply by mass of parent, accumlate a total, then divide by the mass of the cluster
        Vec3<double> massWeightedTotalVec{0,0,0};
        double clusterMass{0}; // Easier the calculate this here rather than messing around with the mass distribution
        const Site* refSite{clusterVec[0]}; // Define a reference site (the first member of the cluster)
        for (const auto& clusterMem : clusterVec)
        {
            // Accumlate mass weighted total vector from reference site
            massWeightedTotalVec += (box->minimumVector(refSite->origin(), clusterMem->origin())) * clusterMem->parent()->parent()->mass();
            clusterMass += clusterMem->parent()->parent()->mass();
        }
        massWeightedTotalVec /= clusterMass;
        clusterCoM_[clusterID] = massWeightedTotalVec;
        double massWeightedDistanceSqrd{0};
        // Now time to calculate the Radius of Gyration
        // Need to run through the clusterMap again, get the mim sqrd distance of site from CoM using ref as origin
        for (const auto& clusterMem : clusterVec)
        {
            massWeightedDistanceSqrd += (box->minimumDistanceSquared((clusterMem->origin() - refSite->origin()), clusterCoM_[clusterID]))
                                            * clusterMem->parent()->parent()->mass();
        }
        radiusOfGyration_[clusterID] = std::sqrt(massWeightedDistanceSqrd / clusterMass);
    }
}