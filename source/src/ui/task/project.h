// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/task/project.h
/// @brief  Project class for ui library.
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu).

#ifndef UI_TASK_PROJECT_H
#define UI_TASK_PROJECT_H

#include <ui/task/project.fwd.h>

#include <ui/task/project_model.h>

//#include <ui/task/node.h>
#include <ui/task/task.fwd.h>

#include <ui/util/exception.h>


#include <QObject>

namespace ui {
namespace task {

///
/// Project is intended to hold all data relevant for current user session.
/// In general saving and later restoring project should allow user to return to exactly the same configuration of tools (and possibly windows).
///
/// Note that this is *not* current configuration of windows (but it can store that as well)
/// instead it should be seen as a `bag` for all relevant user data that was used/created during current session.
///
class Project final : public QObject
{
    Q_OBJECT

public:
	//using Map = std::map< QString, TaskSP >;
	//using Key = Map::key_type;

	//explicit Project(QUuid _project_id);
	explicit Project();

	//std::string type() const override;

	/// Add new Task
    //void add(Key const &, TaskSP const &);
	void add_task(TaskSP const &task);

	// erase given Task from this project
	bool erase(TaskSP const &task);

	// return number of tasks in Project
	//int size() const { return tasks_.size(); }

	// GUI helper function
	// return key of task or nullptr if leaf could not be found
	//Key const * find(Task *leaf) const;

	// Assign 'project_' of given Task to this (assuming it was already inserted into tasks_)
	void assign_ownership(TaskSP const &t);

	// return i'th task from tasks list, return empty SP if index is invalid
	TaskSP task(int index);

	bool operator ==(Project const &r) const;
	bool operator !=(Project const &r) const { return not (*this == r); }


	QString file_name() const { return file_name_; }
	void file_name(QString const &file_name) { file_name_ = file_name; }


	std::vector<TaskSP> const & tasks() const { return tasks_; }

	ProjectTasksModel * model() { return &task_model_; }

	// serialization
	friend QDataStream &operator<<(QDataStream &, Project const&);
	friend QDataStream &operator>>(QDataStream &, Project &);

Q_SIGNALS:


private Q_SLOTS:
	void changed();

private:
	//void listen_to_updates();

	/// assign ownership of all tasks (needed for serialization)
	void assign_ownership();


private:
	std::vector<TaskSP> tasks_;
	ProjectTasksModel task_model_;
	QString file_name_;
};


// class ProjectBadFileFormatException : public ui::util::BadFileFormatException
// {
// public:
//     void raise() const { throw *this; }
//     ProjectBadFileFormatException *clone() const { return new ProjectBadFileFormatException(*this); }
// };


} // namespace task
} // namespace ui

#endif // UI_TASK_PROJECT_H
