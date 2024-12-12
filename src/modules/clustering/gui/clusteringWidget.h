// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2024 Team Dissolve and contributors

#pragma once

#include "modules/clustering/gui/ui_clusteringWidget.h"
#include "modules/widget.h"

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

    public:
    ClusteringModuleWidget(QWidget *parent, ClusteringModule *module, Dissolve &dissolve);

    /*
     * UI
     */
    private:
    // Main form declaration
    Ui::ClusteringModuleWidget ui_;
};
