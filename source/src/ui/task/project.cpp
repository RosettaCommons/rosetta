// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/task/project.cpp
/// @brief  Project class for ui library.
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu).

#include <ui/task/project.h>

#include <ui/task/task.h>

//#include <ui/task/util.h>
#include <ui/util/serialization.h>

#include <cassert>

namespace ui {
namespace task {


//Project::Project(QUuid _project_id) : Node(QUuid()), project_id(_project_id)
Project::Project()
{
	//get_user_credentials();
}


// std::string Project::type() const
// {
// 	return "project";
// }


/// Add new Task
// void Project::add(Key const &key, TaskSP const &task)
// {
// 	if( task->project_ ) task->project_->erase(task);

// 	tasks_[key] = task;
//     task->project_ = this;
// }

/// Add new Task
void Project::add_task(TaskSP const &task)
{
	if( task->project_ ) task->project_->erase(task);

	tasks_.push_back(task);
    task->project_ = this;

	task_model_.update_from_tasks(tasks_);

	connect(task.get(), &Task::changed, this, &Project::changed);
}


// erase given Task from this project
bool Project::erase(TaskSP const &task)
{
	if( task->project_ == this ) {
		auto sz = tasks_.size();

		tasks_.erase( std::remove(tasks_.begin(), tasks_.end(), task), tasks_.end() );

		task->project_ = nullptr;

		disconnect(task.get(), &Task::changed, this, &Project::changed);
		bool res = sz != tasks_.size();

		if(res) Q_EMIT changed();

		return res;
	}
	else return false;
}

void Project::changed()
{
	qDebug() << "Project::changed()";
	task_model_.update_from_tasks(tasks_);
}

// GUI helper function
// return key of task or nullptr if leaf could not be found
// Project::Key const * Project::find(Task *leaf) const
// {
//     auto it = find_if(tasks_.begin(), tasks_.end(), [&](Map::value_type const &p) { return p.second.get() == leaf; } );

// 	if( it == tasks_.end() ) return nullptr;
// 	else return &it->first;
// }


void Project::assign_ownership(TaskSP const &task)
{
	assert( task->project_ == nullptr );
	task->project_ = this;
	connect(task.get(), &Task::changed, this, &Project::changed);
}

TaskSP Project::task(int index)
{
	if(index < 0  or  index >= static_cast<int>( tasks_.size() ) ) return TaskSP();

	return tasks_[index];
}


bool Project::operator ==(Project const &rhs) const
{
	return tasks_ == rhs.tasks_;

	// if( tasks_.size() != rhs.tasks_.size() ) return false;

	// if( not std::equal( tasks_.begin(), tasks_.end(), rhs.tasks_.begin(),
	// 				[](Map::value_type const &l, Map::value_type const &r) {
	// 					if( l.first  != r.first  ) return false;
	// 					if( *l.second != *r.second ) return false;
	// 					return true;
	// 				})
	// 	) {
	// 	return false;
	// }

	// return true;
}


// void Project::listen_to_updates()
// {
// 	for(auto & p : tasks_) p.listen_to_updates();
// }

quint32 const _Project_stream_version_ = 0x00000001;

QDataStream &operator<<(QDataStream &out, Project const &p)
{
	out.setVersion(QDataStream::Qt_5_6);

	out << ui::util::_Project_magic_number_;
	out << _Project_stream_version_;

	using namespace ui::util;
	out << p.tasks_;

	return out;
}


QDataStream &operator>>(QDataStream &in, Project &p)
{
	in.setVersion(QDataStream::Qt_5_6);

	quint64 magic;
	in >> magic;
	if( magic != ui::util::_Project_magic_number_ ) throw ui::util::BadFileFormatException( QString("Invalid _Project_magic_number_: read %1, was expecting %2...").arg(magic).arg(ui::util::_Project_magic_number_) );

	quint32 version;
	in >> version;
	if( version != _Project_stream_version_ ) throw ui::util::BadFileFormatException( QString("Invalid _Project_stream_version_: read %1, was expecting %2...").arg(magic).arg(_Project_stream_version_) );

	using namespace ui::util;
	in >> p.tasks_;

	for(auto & t : p.tasks_) {
		p.assign_ownership(t);
		//listen_to_updates();  // moved to Task operator>>
	}

	p.task_model_.update_from_tasks(p.tasks_);

	return in;
}



} // namespace task
} // namespace ui
