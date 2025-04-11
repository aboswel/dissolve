// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

#include "modules/clustering/gui/clusteringWidget.h"
#include "gui/render/renderableData1D.h"
#include "modules/clustering/clustering.h"

ClusteringModuleWidget::ClusteringModuleWidget(QWidget *parent, ClusteringModule *module, Dissolve &dissolve)
    : ModuleWidget(parent, dissolve), module_(module)
{
    // Set up cluster configuration viewing
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
    // Set the new configuration and redraw the visualisation. Is there a way to set the style to spheres view from here (lines
    // arent showing)?
    clusterConfiguration_ = module_->generateClustersConfig(&dissolve_, module_->getSourceConfig(), displaySize_, displayID_);
    ui_.ViewerWidget->setConfiguration(clusterConfiguration_);
    ui_.ViewerWidget->dataModified();
    ui_.ViewerWidget->postRedisplay();

    // Configure the size/mass histograms
    if (updateFlags.isSet(ModuleWidget::RecreateRenderablesFlag) || sizeDist_->renderables().empty() ||
        massDist_->renderables().empty())
    {
        buildSizeList();
        sizeDist_->clearRenderables();
        massDist_->clearRenderables();

        if (sizeDist_->renderables().empty())
            sizeDist_->createRenderable<RenderableData1D>(std::format("{}//SizeDist", module_->name()),
                                                          std::format("SizeDist//{}", module_->getSourceConfig()->niceName()),
                                                          module_->getSourceConfig()->niceName());

        if (massDist_->renderables().empty())
            massDist_->createRenderable<RenderableData1D>(std::format("{}//MassDist", module_->name()),
                                                          std::format("MassDist//{}", module_->getSourceConfig()->niceName()),
                                                          module_->getSourceConfig()->niceName());
    }

    sizeDist_->validateRenderables(dissolve_.processingModuleData());
    massDist_->validateRenderables(dissolve_.processingModuleData());
    ui_.SizePlot->updateToolbar();
    ui_.MassPlot->updateToolbar();
    sizeDist_->postRedisplay();
    massDist_->postRedisplay();

    refreshing_ = false;
}

void ClusteringModuleWidget::buildSizeList()
{
    // Produce a list of all the different cluster sizes
    auto sizeDistribution = module_->getSizeDistribution();
    ui_.listWidget->clear();
    ui_.listWidget->addItem("View All");
    for (const auto &[size, _] : sizeDistribution)
    {
        QListWidgetItem *item = new QListWidgetItem(QString::number(size));
        item->setFlags(item->flags() | Qt::ItemIsEnabled | Qt::ItemIsSelectable);
        ui_.listWidget->addItem(item);
    }
}

void ClusteringModuleWidget::buildIDList(QListWidgetItem *item)
{
    // Produce a list of all the cluster IDs of the selected cluster size
    auto sizeDistribution = module_->getSizeDistribution();
    ui_.listWidget2->clear();
    if (item->text() == "View All")
    {
        displaySize_ = 0;
        updateControls();
        return;
    }

    ui_.listWidget2->addItem("View All");
    displaySize_ = item->text().toInt();
    for (const auto &ID : sizeDistribution[displaySize_])
    {
        QListWidgetItem *cell = new QListWidgetItem(QString::number(ID));
        cell->setFlags(cell->flags() | Qt::ItemIsEnabled | Qt::ItemIsSelectable);
        ui_.listWidget2->addItem(cell);
    }
}

void ClusteringModuleWidget::on_listWidget_itemClicked(QListWidgetItem *item) { buildIDList(item); }

void ClusteringModuleWidget::on_listWidget2_itemClicked(QListWidgetItem *item)
{
    if (item->text() == "View All")
    {
        displayID_ = 0;
        updateControls();
        return;
    }

    displayID_ = item->text().toInt();
    updateControls();
}