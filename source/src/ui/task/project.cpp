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

#include <ui/task/util.h>

#include <cassert>

namespace ui {
namespace task {


//Project::Project(QUuid _project_id) : Node(QUuid()), project_id(_project_id)
Project::Project()
{
}


// std::string Project::type() const
// {
// 	return "project";
// }


/// Add new Task
void Project::add(Key const &key, TaskSP const &task)
{
	if( task->project_ ) task->project_->erase( task.get() );

	tasks_[key] = task;
    task->project_ = this;
}


// erase given Task from this project
bool Project::erase(Task *task)
{
    auto it = find_if(tasks_.begin(), tasks_.end(), [&](Map::value_type const &p) { return p.second.get() == task; } );

	if( it == tasks_.end() ) return false;

    tasks_.erase(it);
	task->project_ = nullptr;
   	return true;
}

// GUI helper function
// return key of task or nullptr if leaf could not be found
Project::Key const * Project::find(Task *leaf) const
{
    auto it = find_if(tasks_.begin(), tasks_.end(), [&](Map::value_type const &p) { return p.second.get() == leaf; } );

	if( it == tasks_.end() ) return nullptr;
	else return &it->first;
}


void Project::assign_ownership(TaskSP const &t)
{
	assert( t->project_ == nullptr );

	t->project_ = this;
}


bool Project::operator ==(Project const &rhs) const
{


	if( tasks_.size() != rhs.tasks_.size() ) return false;

	if( not std::equal( tasks_.begin(), tasks_.end(), rhs.tasks_.begin(),
					[](Map::value_type const &l, Map::value_type const &r) {
						if( l.first  != r.first  ) return false;
						if( *l.second != *r.second ) return false;
						return true;
					})
		) {
		return false;
	}

	return true;
}




quint64 const _Project_magic_number_   = 0xFFF1D6BF94D5E57F;
quint32 const _Project_stream_version_ = 0x00000001;


QDataStream &operator<<(QDataStream &out, Project const &p)
{
	out.setVersion(QDataStream::Qt_5_6);

	out << _Project_magic_number_;
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
	if( magic != _Project_magic_number_ ) throw ui::util::BadFileFormatException();  // ProjectBadFileFormatException();

	quint32 version;
	in >> version;
	if( version != _Project_stream_version_ ) throw ui::util::BadFileFormatException();  // ProjectBadFileFormatException();

	in >> p.tasks_;

	for(auto &t : p.tasks_) p.assign_ownership(t.second);

	return in;
}



} // namespace task
} // namespace ui
