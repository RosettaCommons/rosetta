// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/task/project_view.cpp
/// @brief  Project view class for ui library.
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu).

#include "project_view.h"
#include "ui_project_view.h"

#include <ui/task/project_model.h>


namespace ui {
namespace task {


ProjectView::ProjectView(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::ProjectView)
{
    ui->setupUi(this);

	ProjectModel *model = new ProjectModel();
        ui->project->setModel(model);

}

ProjectView::~ProjectView()
{
    delete ui;
}

} // namespace task
} // namespace ui
