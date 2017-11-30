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

#include <ui/task/task.h>
#include <ui/task/task_view.h>
#include <ui/config/config_dialog.h>

#include <QFileDialog>
#include <QDebug>

namespace ui {
namespace task {


ProjectView::ProjectView(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::ProjectView)
{
    ui->setupUi(this);

	project_model_ = new ProjectModel(this);
	ui->project->setModel(project_model_);

	project_ = std::make_shared<Project>();
	project_model_->set(project_);

	// QItemSelectionModel *selectionModel = new QItemSelectionModel(model);
	// ui->project->setSelectionModel(selectionModel);
	// ui->project->selectionModel()->setCurrentIndex(ui->project->currentIndex(), QItemSelectionModel::SelectCurrent);
}

ProjectView::~ProjectView()
{
    delete ui;
}


void ProjectView::on_new_task_clicked()
{
    qDebug() << "ProjectView::on_new_task_clicked()";

	static int task_counter = 0;
	//QString label = QString("Task-%1").arg( project_->size() );
	QString label = QString("Task-%1").arg( task_counter++ );

	project_->add(label, std::make_shared<Task>(label) );
	project_model_->set(project_);
}

/*
void ProjectView::on_submit_clicked()
{
    qDebug() << "ProjectView::on_submit_clicked()";

	if( ProjectModel * model = static_cast<ProjectModel*>( ui->project->model() ) ) {
		bool valid = false;

        QModelIndex index = ui->project->currentIndex();

		// check if selection is actually highlighted to avoid accidents
		QModelIndexList list = ui->project->selectionModel()->selectedIndexes();
		for(QModelIndex const &i : list) {
			//qDebug() << "row:" << i.row() << " column:" << i.column();
			if( index.row() == i.row() ) valid = true;

			//if( index.row() == i.row() and index.column() == i.column()  and  i.column() == 0 ) valid = true;
		}

		//qDebug() << "currentIndex: " << index;

		if( valid and index.isValid() ) {
			if( auto task = model->task(index) ) task->submit();
		}
	}
}
*/
void ProjectView::on_delete_task_clicked()
{
    qDebug() << "ProjectView::on_delete_task_clicked()";

	if( ProjectModel * model = static_cast<ProjectModel*>( ui->project->model() ) ) {
		bool valid = false;

        QModelIndex index = ui->project->currentIndex();

		// check if selection is actually highlighted to avoid accidents
		QModelIndexList list = ui->project->selectionModel()->selectedIndexes();
		for(QModelIndex const &i : list) {
			//qDebug() << "row:" << i.row() << " column:" << i.column();
			if( index.row() == i.row() ) valid = true;

			//if( index.row() == i.row() and index.column() == i.column()  and  i.column() == 0 ) valid = true;
		}

		//qDebug() << "currentIndex: " << index;

		if( valid and index.isValid() ) {
			if( auto task = model->task(index) ) {
				project_->erase(task);
				project_model_->set(project_);
			}
		}
	}
}


void ProjectView::on_action_new_project_triggered()
{
	project_ = std::make_shared<Project>();

	project_model_->set(project_);

	file_name_ = "";
}


void ProjectView::on_action_open_project_triggered()
{
    qDebug() << "ProjectView::on_project_open()";

	QString file_name = QFileDialog::getOpenFileName(this, tr("Open Project"), "", tr("RosettaUI Projects (*.rosetta)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
	if( not file_name.isEmpty() ) {
		QFile file(file_name);

		if (!file.open(QIODevice::ReadOnly) ) return;

		QDataStream in( &file );

		project_ = std::make_shared<Project>();
		in >> *project_;
		file.close();

		project_model_->set(project_);

		file_name_ = file_name;
	}
}

void ProjectView::save_project(bool avoid_asking_file_name_if_possible)
{

	if( file_name_.isEmpty()  or  (not avoid_asking_file_name_if_possible) ) {
		file_name_ = QFileDialog::getSaveFileName(this, tr("Save File"), "", tr("RosettaUI Projects (*.rosetta)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
	}

	if( not file_name_.isEmpty() ) {
		QFile file(file_name_);

		if (!file.open(QIODevice::WriteOnly) ) return;

		QDataStream out( &file );

		out << *project_;
		file.close();
	}
}

void ProjectView::on_action_save_project_triggered()
{
    qDebug() << "ProjectView::on_project_save()";
	save_project(true);
}

void ProjectView::on_action_save_project_as_triggered()
{
    qDebug() << "ProjectView::on_project_save_as()";
	save_project(false);
}


void ProjectView::on_action_preferences_triggered()
{
    qDebug() << "ProjectView::on_action_preferences_triggered()";

        config::ConfigDialog * preferences = new config::ConfigDialog(this);
    preferences->show();
}



} // namespace task
} // namespace ui

void ui::task::ProjectView::on_project_doubleClicked(const QModelIndex &index)
{
    qDebug() << "ProjectView::on_project_doubleClicked";

	//QModelIndex index = ui->project_view->currentIndex();
	TaskSP task = project_model_->task(index);

	if(task) {
		auto tv = new TaskView(task);
		tv->show();
	}
}
