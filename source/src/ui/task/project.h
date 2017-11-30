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


//#include <ui/task/node.h>
#include <ui/task/task.fwd.h>

#include <ui/util/exception.h>


#include <QObject>

namespace ui {
namespace task {


class Project;
using ProjectSP  = std::shared_ptr< Project >;
using ProjectCSP = std::shared_ptr< Project const >;


class Project final : public QObject
{
    Q_OBJECT

public:
	using Map = std::map< QString, TaskSP >;
	using Key = Map::key_type;

	//explicit Project(QUuid _project_id);
	explicit Project();

	//std::string type() const override;


	/// Add new Task
    void add(Key const &, TaskSP const &);

	// erase given Task from this project
	bool erase(TaskSP const &task);

	// return number of tasks in Project
	int size() const { return tasks_.size(); }


	// GUI helper function
	// return key of task or nullptr if leaf could not be found
	Key const * find(Task *leaf) const;

	// Assign 'project_' of given Task to this (assuming it was already inserted into tasks_)
	void assign_ownership(TaskSP const &t);

	bool operator ==(Project const &r) const;
	bool operator !=(Project const &r) const { return not (*this == r); }

	// serialization
	friend QDataStream &operator<<(QDataStream &, Project const&);
	friend QDataStream &operator>>(QDataStream &, Project &);

private:
	//void listen_to_updates();

private:

	Map tasks_;

	friend struct PNode;
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
