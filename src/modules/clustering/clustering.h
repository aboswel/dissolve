// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2024 Team Dissolve and contributors

#pragma once

#include "module/module.h"

// Clustering Module
class ClusteringModule : public Module
{
    public:
    ClusteringModule();
    ~ClusteringModule() override = default;

    /*
     * Definition
     */
    private:
    // Target configurations
    Configuration *targetConfiguration_{nullptr};

    /*
     * Processing
     */
    private:
    // Run main processing
    Module::ExecutionResult process(ModuleContext &moduleContext) override;
};
