// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/task/project_model.cpp
/// @brief  Project model class for ui library.
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu).

#include <ui/task/project_model.h>

#include <ui/task/task.h>

#include <QDebug>


namespace ui {
namespace task {

QVariant ProjectTasksModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (role == Qt::DisplayRole) {
		if (orientation == Qt::Horizontal) {
			switch (section) {
				case 0: return QString("name");
				case 1: return QString("state");
			}
		} else if (orientation == Qt::Vertical) {
			return QString::number(section+1);
		}
	}
    return QVariant();
}


int ProjectTasksModel::rowCount(const QModelIndex &/*parent*/) const
{
	return rows_.size();
}


int	ProjectTasksModel::columnCount(const QModelIndex &/*parent*/) const
{
	return 2;
}


Qt::ItemFlags ProjectTasksModel::flags(const QModelIndex &index) const
{
    if (!index.isValid()) return Qt::NoItemFlags;
	//if( index.column() == 0  and  editable_ ) return Qt::ItemIsEditable | QAbstractItemModel::flags(index);
	return QAbstractItemModel::flags(index);
}


QVariant ProjectTasksModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid()) return QVariant();
    //if (role != Qt::DisplayRole && role != Qt::EditRole) return QVariant();
	if (role != Qt::DisplayRole ) return QVariant();

	if( index.row() < static_cast<int>( rows_.size() ) ) {

		if      ( index.column() == 0 ) return rows_[index.row()].name;
		else if ( index.column() == 1 ) return rows_[index.row()].state;
	}
	return QVariant();
}


void ProjectTasksModel::update_from_tasks(std::vector<TaskSP> const &tasks)
{
	beginResetModel();

	rows_.clear();
	rows_.reserve( tasks.size() );

	for (auto const & t : tasks) rows_.push_back( Row{t->name(), Task::to_string( t->state() ) } );

	endResetModel();
}


} // namespace task
} // namespace ui
