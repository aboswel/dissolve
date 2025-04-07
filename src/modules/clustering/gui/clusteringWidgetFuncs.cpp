// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

#include "modules/clustering/clustering.h"
#include "modules/clustering/gui/clusteringWidget.h"

ClusteringModuleWidget::ClusteringModuleWidget(QWidget *parent, ClusteringModule *module, Dissolve &dissolve)
    : ModuleWidget(parent, dissolve), module_(module)
{
    // Set up user interface
    ui_.setupUi(this);
    ui_.ViewerWidget->setConfiguration(clusterConfiguration_);
}

void ClusteringModuleWidget::updateControls()
{
    // Check eveything is ready
    if (!module_->viewingReady)
        return;

    if (refreshing_)
        return;

    refreshing_ = true;

    // Set the new configuration and redraw the visualisation
    clusterConfiguration_ = module_->generateClustersConfig(&dissolve_, module_->getSourceConfig());
    ui_.ViewerWidget->setConfiguration(clusterConfiguration_);
    ui_.ViewerWidget->dataModified(); // still needed?
    ui_.ViewerWidget->postRedisplay();

    refreshing_ = false;
}

void ClusteringModuleWidget::on_refreshButton_clicked() { updateControls(); }
