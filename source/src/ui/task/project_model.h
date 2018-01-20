// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/task/project_model.h
/// @brief  Project model class for ui library.
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu).

#ifndef UI_TASK_PROJECTMODEL_H
#define UI_TASK_PROJECTMODEL_H

#include <ui/task/task.fwd.h>

#include <QAbstractTableModel>

namespace ui {
namespace task {

class ProjectTasksModel : public QAbstractTableModel
{
    Q_OBJECT

public:
	using QAbstractTableModel::QAbstractTableModel;

	QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;

	int rowCount(const QModelIndex &parent = QModelIndex()) const override;
	int	columnCount(const QModelIndex &parent = QModelIndex()) const override;
	Qt::ItemFlags flags(const QModelIndex &index) const override;

	QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;

	void update_from_tasks(std::vector<TaskSP> const &tasks);
private:
	struct Row {
		QString name, state;
	};

	std::vector<Row> rows_;
};

} // namespace task
} // namespace ui

#endif // UI_TASK_PROJECTMODEL_H
