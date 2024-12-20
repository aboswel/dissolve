// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2024 Team Dissolve and contributors

// Constructor to initialise member variables with default values, and allocate data structures for storing results
// Destructor to clean up once analysis is done

// Move user settings from GUI into member variables, validate inputs
// Set up initial references to system objects i.e. pointers to simulation box

// Load system data into relevant data structures

// Basically ensure everything is ready for process.cpp to actually perfrom the analysis
// clustering.cpp
#include "modules/clustering/clustering.h"

ClusteringModule::ClusteringModule() : Module(ModuleTypes::Clustering) 
{

}
