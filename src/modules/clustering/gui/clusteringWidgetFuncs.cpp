// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (c) 2025 Team Dissolve and contributors

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
    // Check eveything is ready
    if (!module_->viewingReady)
        return;

    if (refreshing_)
        return;

    // Set the new configuration and redraw the visualisation
    clusterConfiguration_ = module_->generateClustersConfig(&dissolve_, module_->getSourceConfig(), displaySize_, displayID_);
    ui_.ViewerWidget->setConfiguration(clusterConfiguration_);
    ui_.ViewerWidget->dataModified(); // still needed?
    ui_.ViewerWidget->postRedisplay();
}

void ClusteringModuleWidget::buildSizeList()
{
    // Produce a list of all the different cluster sizes
    auto sizeDistribution = module_->getSizeDistribution();
    ui_.listWidget->clear();
    int row{0};
    for (const auto &[size, _] : sizeDistribution)
    {
        QListWidgetItem *item = new QListWidgetItem(QString::number(size));
        item->setFlags(item->flags() | Qt::ItemIsEnabled | Qt::ItemIsSelectable);
        ui_.listWidget->addItem(item);
        row++;
    }
    ui_.listWidget->addItem("View All");
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

    displaySize_ = item->text().toInt();
    int r{0};
    for (const auto &ID : sizeDistribution[displaySize_])
    {
        QListWidgetItem *cell = new QListWidgetItem(QString::number(ID));
        cell->setFlags(cell->flags() | Qt::ItemIsEnabled | Qt::ItemIsSelectable);
        ui_.listWidget2->addItem(cell);
        r++;
    }
    ui_.listWidget2->addItem("View All");
}

void ClusteringModuleWidget::on_refreshButton_clicked()
{
    // Build the size table
    buildSizeList();
    updateControls();
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