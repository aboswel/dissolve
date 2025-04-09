// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

#pragma once

#include "gui/configurationViewer.h"
#include "gui/configurationViewerWidget.h"
#include "modules/clustering/gui/ui_clusteringWidget.h"
#include "modules/widget.h"
#include <QListWidget>

// Forward Declarations
class ClusteringModule;

// Module Widget
class ClusteringModuleWidget : public ModuleWidget
{
    Q_OBJECT

    private:
    // Associated Module
    ClusteringModule *module_;
    // Pointer to cluster config for visualisation (maybe just use the module config?)
    Configuration *clusterConfiguration_{nullptr};
    // Integers for which clusters to display
    int displaySize_{0}, displayID_{0};

    public:
    ClusteringModuleWidget(QWidget *parent, ClusteringModule *module, Dissolve &dissolve);

    private:
    void updateControls();
    void buildSizeList();
    void buildIDList(QListWidgetItem *item);

    /*
     * UI
     */
    private:
    // Main form declaration
    Ui::ClusteringModuleWidget ui_;

    private Q_SLOTS:
    void on_refreshButton_clicked();
    void on_listWidget_itemClicked(QListWidgetItem *item);
    void on_listWidget2_itemClicked(QListWidgetItem *item);
};
