// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2024 Team Dissolve and contributors

#include "modules/clustering/gui/clusteringWidget.h"
#include "modules/clustering/clustering.h"

ClusteringModuleWidget::ClusteringModuleWidget(QWidget *parent, ClusteringModule *module, Dissolve &dissolve)
    : ModuleWidget(parent, dissolve), module_(module)
{
    // Set up user interface
    ui_.setupUi(this);

    refreshing_ = false;
}