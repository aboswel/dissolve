// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

#include "gui/render/renderableData1D.h"
#include "modules/clustering/clustering.h"
#include "modules/clustering/gui/clusteringWidget.h"

ClusteringModuleWidget::ClusteringModuleWidget(QWidget *parent, ClusteringModule *module, Dissolve &dissolve)
    : ModuleWidget(parent, dissolve), module_(module)
{
    ui_.setupUi(this);
    ui_.ViewerWidget->setConfiguration(clusterConfiguration_);

    // Set up the distribution plots
    sizeDist_ = ui_.SizePlot->dataViewer();
    auto &viewS = sizeDist_->view();
    viewS.setViewType(View::FlatXYView);
    viewS.axes().setTitle(0, "Log(Cluster Size)");
    viewS.axes().setMax(0, 10.0);
    viewS.axes().setTitle(1, "Log(No. of Clusters)");
    viewS.axes().setMin(1, 0.0);
    viewS.axes().setMax(1, 1.0);
    viewS.setAutoFollowType(View::AllAutoFollow);

    massDist_ = ui_.MassPlot->dataViewer();
    auto &viewM = massDist_->view();
    viewM.setViewType(View::FlatXYView);
    viewM.axes().setTitle(0, "Log(Cluster Mass)");
    viewM.axes().setMax(0, 10.0);
    viewM.axes().setTitle(1, "Log(No. of Clusters)");
    viewM.axes().setMin(1, 0.0);
    viewM.axes().setMax(1, 1.0);
    viewM.setAutoFollowType(View::AllAutoFollow);
}

void ClusteringModuleWidget::updateControls(const Flags<ModuleWidget::UpdateFlags> &updateFlags)
{
    if (refreshing_)
        return;

    refreshing_ = true;

    // Check eveything is ready
    if (!module_->viewingReady)
    {
        refreshing_ = false;
        return;
    }

    module_->generateClustersConfig(dissolve_, displaySize_, displayID_);
    clusterConfiguration_ = module_->getClusterConfig();
    ui_.ViewerWidget->setConfiguration(clusterConfiguration_);
    ui_.ViewerWidget->dataModified();
    ui_.ViewerWidget->postRedisplay();

    // Configure the size/mass histograms
    if (updateFlags.isSet(ModuleWidget::RecreateRenderablesFlag) || sizeDist_->renderables().empty() ||
        massDist_->renderables().empty())
    {        
        Messenger::error("BUILD LIST SHOULD BE UPDATED??????????");
        sizeDist_->clearRenderables();
        massDist_->clearRenderables();

        if (sizeDist_->renderables().empty())
            sizeDist_->createRenderable<RenderableData1D>(
                std::format("{}//SizeDist", module_->name()),
                std::format("SizeDist//{}", module_->keywords().getConfiguration("Configuration")->niceName()),
                module_->keywords().getConfiguration("Configuration")->niceName());

        if (massDist_->renderables().empty())
            massDist_->createRenderable<RenderableData1D>(
                std::format("{}//MassDist", module_->name()),
                std::format("MassDist//{}", module_->keywords().getConfiguration("Configuration")->niceName()),
                module_->keywords().getConfiguration("Configuration")->niceName());
    }

    if (!fromBuilder)
        buildSizeList();

    sizeDist_->validateRenderables(dissolve_.processingModuleData());
    massDist_->validateRenderables(dissolve_.processingModuleData());
    ui_.SizePlot->updateToolbar();
    ui_.MassPlot->updateToolbar();
    sizeDist_->postRedisplay();
    massDist_->postRedisplay();

    fromBuilder = false;
    refreshing_ = false;
}

void ClusteringModuleWidget::buildSizeList()
{
    // Produce a list of all the different cluster sizes
    auto &sizeDistribution = module_->getSizeDistribution();
    ui_.clusterSizeList->clear();
    ui_.clusterSizeList->addItem("View All");
    for (const auto &[size, _] : sizeDistribution)
    {
        auto *item = new QListWidgetItem(QString::number(size));
        item->setFlags(item->flags() | Qt::ItemIsEnabled | Qt::ItemIsSelectable);
        ui_.clusterSizeList->addItem(item);
    }
}

void ClusteringModuleWidget::buildIDList(QListWidgetItem *item)
{
    // Produce a list of all the cluster IDs of the selected cluster size
    auto sizeDistribution = module_->getSizeDistribution();
    ui_.clusterIDList->clear();
    if (item->text() == "View All")
    {
        displaySize_ = 0;
        updateControls();
        return;
    }

    ui_.clusterIDList->addItem("View All");
    displaySize_ = item->text().toInt();
    for (const auto &ID : sizeDistribution[displaySize_])
    {
        auto *cell = new QListWidgetItem(QString::number(ID));
        cell->setFlags(cell->flags() | Qt::ItemIsEnabled | Qt::ItemIsSelectable);
        ui_.clusterIDList->addItem(cell);
    }
}

void ClusteringModuleWidget::on_clusterSizeList_itemClicked(QListWidgetItem *item) { buildIDList(item); }

void ClusteringModuleWidget::on_clusterIDList_itemClicked(QListWidgetItem *item)
{
    item->text() == "View All" ? displayID_ = 0 : displayID_ = item->text().toInt();
    fromBuilder = true;
    updateControls();
}