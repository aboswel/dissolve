// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

#include "modules/clustering/clustering.h"
#include "keywords/double.h"
#include "keywords/integer.h"
#include "keywords/configuration.h"
#include "keywords/speciesSite.h"

ClusteringModule::ClusteringModule() : Module(ModuleTypes::Clustering) 
{
    // Run clustering despite no configuration changes (to investigate different cluster definitions on the same config)
    // Perhaps add an option to disable evolution from here for convenience?
    executeIfTargetsUnchanged_ = true;

    keywords_.addTarget<ConfigurationKeyword>("Configuration", "Set target configuration for the module", targetConfiguration_)
        ->setEditSignals({KeywordBase::ClearModuleData, KeywordBase::RecreateRenderables});

    keywords_.setOrganisation("Options", "Sites");
    keywords_.add<SpeciesSiteKeyword>("SiteA", "Choose the first site for cluster definition", userSites_.a_);
    keywords_.add<SpeciesSiteKeyword>("SiteB", "Choose the second site for cluster definition", userSites_.b_);
    keywords_.add<DoubleKeyword>("Cut-off", "The maximum distance between sites for them to be considered clustered", userSites_.cutOff);

    keywords_.setOrganisation("Options", "Display");
    keywords_.add<IntegerKeyword>("Minimum cluster size to display", "No clusters below this threshold will show in the output visualiser", minClusterForConfig_);
    keywords_.add<IntegerKeyword>("Maximum cluster size to display", "No clusters above this threshold will show in the output visualiser", maxClusterForConfig_);

}
