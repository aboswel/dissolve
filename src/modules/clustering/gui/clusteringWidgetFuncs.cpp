// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2024 Team Dissolve and contributors

#include "modules/clustering/gui/clusteringWidget.h"
#include "modules/clustering/clustering.h"

ClusteringModuleWidget::ClusteringModuleWidget(QWidget *parent, ClusteringModule *module, Dissolve &dissolve)
    : ModuleWidget(parent, dissolve), module_(module)
{
    // Set up user interface
    ui_.setupUi(this);
    ui_.ViewerWidget->setConfiguration(clusterConfiguration_);
}

void ClusteringModuleWidget::updateControls()
{
    module_->functionTrack("NewConfiguration1.neighbourdata.txt", 1);
    if (!module_->viewingReady)
        return;

    if (refreshing_)
        return;
    
    refreshing_ = true;
    
    clusterConfiguration_ = module_->generateClustersConfig(&dissolve_, module_->getSourceConfig());
    ui_.ViewerWidget->setConfiguration(clusterConfiguration_);
    ui_.ViewerWidget->dataModified();
    ui_.ViewerWidget->postRedisplay();

    refreshing_ = false;
}

void ClusteringModuleWidget::on_refreshButton_clicked()
{
    updateControls();
}
