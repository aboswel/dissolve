// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2024 Team Dissolve and contributors

#pragma once

#include "gui/configurationViewer.h"
#include "modules/clustering/gui/ui_clusteringWidget.h"
#include "modules/widget.h"
#include "gui/configurationViewerWidget.h"

// Forward Declarations
class ClusteringModule;

// Module Widget
class ClusteringModuleWidget : public ModuleWidget
{
    // All Qt declarations derived from QObject must include this macro
    Q_OBJECT

    private:
    // Associated Module
    ClusteringModule *module_;
    ConfigurationViewerWidget *viewerWidget_{nullptr};
    Configuration *clusterConfiguration_{nullptr}; 
    bool refreshing_{false};

    public:
    ClusteringModuleWidget(QWidget *parent, ClusteringModule *module, Dissolve &dissolve);
    
    private:
    void updateControls();

    /*
     * UI
     */
    private:
    // Main form declaration
    Ui::ClusteringModuleWidget ui_;

    private Q_SLOTS:
    void on_refreshButton_clicked();

    
};
