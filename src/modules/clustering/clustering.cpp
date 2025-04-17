// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

#include "keywords/bool.h"
#include "keywords/configuration.h"
#include "keywords/double.h"
#include "keywords/expression.h"
#include "keywords/integer.h"
#include "keywords/speciesSite.h"
#include "modules/clustering/clustering.h"

ClusteringModule::ClusteringModule() : Module(ModuleTypes::Clustering)
{
    // Run clustering despite no configuration changes (to investigate different cluster definitions on the same config)
    // Perhaps add an option to disable evolution from here for convenience?
    executeIfTargetsUnchanged_ = true;

    keywords_.addTarget<ConfigurationKeyword>("Configuration", "Set target configuration for the module", targetConfiguration_)
        ->setEditSignals({KeywordBase::ClearModuleData, KeywordBase::RecreateRenderables});

    keywords_.setOrganisation("Options", "Sites");
    keywords_.add<SpeciesSiteKeyword>("SiteA", "Choose the first site for cluster definition", a_);
    keywords_.add<SpeciesSiteKeyword>("SiteB", "Choose the second site for cluster definition", b_);
    keywords_.add<DoubleKeyword>(
        "Cut-off", "The maximum distance between sites for them to be considered part of the same cluster", cutoff_);

    keywords_.setOrganisation("Options", "Hydrogen Bonding");
    keywords_.add<BoolKeyword>("Directional hydrogen bonding (hover for info)?",
                               "Must be defined between two oxygens, with at least one fragment site specifing a hydroxyl "
                               "group, use -O(-H) NETA description for general hydroxyl groups when generating an empty "
                               "fragment string on your species",
                               strict_);
    keywords_.add<DoubleKeyword>("Maximum angle deviation from 180", "+/- from 180 degrees", angleDev_);

    keywords_.setOrganisation("Export", "Options");
    keywords_.add<BoolKeyword>("Export Size Distribution", "SizeDist.txt", saveSizeDist_);
    keywords_.add<BoolKeyword>("Export Mass Distribution", "MassDist.txt", saveMassDist_);
    keywords_.add<BoolKeyword>("Export Radius of Gyration - Cluster Mass", "RgMass.txt", saveRgMass_);
}
