// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/task/task.h
/// @brief  Task class for ui library.
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu).

#ifndef TASK_H
#define TASK_H

#include <ui/task/task.fwd.h>

#include <ui/task/node.h>
#include <ui/task/project.h>

#include <ui/util/exception.h>

#include <QObject>

namespace ui {
namespace task {


/// Generic container that represent abstract chunk of data that will be stored as file on back-end.
/// Derive new classes from it if you need to create new types of Node while re-use serialization routines
class FileMixin
{
public:

	virtual QString type() const = 0;

private:
	virtual QByteArray data() const = 0;
	virtual void data(QByteArray const &) = 0;

	//virtual QByteArray file_data() const = 0;
	//virtual void file_data(QByteArray const &) = 0;

	virtual QString file_name() const = 0;


	// return true if underlying object contain no data
	virtual bool empty() const = 0;
};



/// File node that hold concreate data with concreate file-name
class File : public FileMixin
{
public:
	explicit File() {}
	explicit File(QString const & file_name);
	explicit File(QString const & file_name, QByteArray const & file_data) : file_name_(file_name), file_data_(file_data) {}

	QString type() const override { return "file"; };

	void init_from_file(QString const & file_name);

	// we pulling FileMixin methods into public because they now have double function as accessor to file data
	QByteArray data() const override { return file_data_; };
	void data(QByteArray const &_file_data) override { file_data_ = _file_data; }

    //QByteArray file_data() const override { return file_data_; };
	//void file_data(QByteArray const &_file_data) override { file_data_ = _file_data; }

	QString file_name() const override { return file_name_; }
	void file_name(QString const &_file_name) { file_name_ = _file_name; }

	bool empty() const override { return file_data_.isEmpty(); }
	//bool null() const override { return file_data_.isNull(); }

	File& operator=(File&& other) noexcept;


	bool operator ==(File const &r) const;
	bool operator !=(File const &r) const { return not (*this == r); }

	// serialization
	friend QDataStream &operator<<(QDataStream &, File const&);
	friend QDataStream &operator>>(QDataStream &, File &);

private:

	QString file_name_;
	QByteArray file_data_;
};



class Task : public QObject
{
	Q_OBJECT

public:
	enum class State {_draft_, _queued_, _running_, _finished_};

public:
	explicit Task() {}
	explicit Task(QString const &description);

	//std::string type() const override;

	/// return reference to 'input' node which will contain all input data
	//Node &input();

	/// return reference to 'output' node which will contain all output data
	//Node &output();

	/// return reference to 'xml script' node which will contain all Rosetta XML script for this task
	//File &script();

	// return current Task State
	QString state() const;

	Project *project() const { return project_; }

	// should not be needed, add/remove Tasks should be done on Project level
	// void project(Project *p)  { project_ = p; }

	/// Initiate submit procedure
	///   0. set state to _draft_
	///   1. sync all nodes to cloud
	///   2. set state to _queued_
	///   3. sync root node to cloud
	void submit();

	/// return task UUID if task was already submitted other wise return NULL UUID
	QUuid task_id() const;

	QString description() const { return description_; }
	void description(QString const &d);


	File const &input() const { return input_; }
	void input(File &&input) { if( input_ != input ) { input_ = std::move(input);  Q_EMIT changed(); } }

	File const &script() const { return script_; }
	void script(File &&script) { if( script_ != script ) { script_ = std::move(script); Q_EMIT changed(); } }

	File const &flags() const { return flags_; }
	void flags(File &&flags) { if( flags_ != flags ) { flags_ = std::move(flags); Q_EMIT changed(); } }


	bool operator ==(Task const &r) const;
	bool operator !=(Task const &r) const { return not (*this == r); }

	// serialization
	friend QDataStream &operator<<(QDataStream &, Task const&);
	friend QDataStream &operator>>(QDataStream &, Task &);

Q_SIGNALS:
	void submitted();

    void changed();

private Q_SLOTS:
	void draft_to_queued(void);

	void post_submit(void);


private:
	void create_sync_tree();

	void assign_input_node();
	void assign_script_node();
	void assign_flags_node();
	void assign_output_node();

private:
	QVariant task_data();
	void task_data(QVariant &&);

	/// Some project functions require access to project_
	friend void Project::add(Project::Key const &, TaskSP const &);
	friend bool Project::erase(Task *task);
	friend void Project::assign_ownership(TaskSP const &t);

	State state_ = State::_draft_;

	/// UUID for cloud syncing. Zero id by default (will be set to a new random value on 'submit' event).
	//QUuid task_id_;

	QString description_;

	File input_, script_, flags_;

	NodeWP input_node_, script_node_, flags_node_;
	NodeWP output_node_;

	/// root node of network syncing tree
	NodeSP root_;

	QPointer<Project> project_;
};



// class TaskBadFileFormatException : public ui::util::BadFileFormatException
// {
// public:
//     void raise() const { throw *this; }
//     TaskBadFileFormatException *clone() const { return new TaskBadFileFormatException(*this); }
// };


} // namespace task
} // namespace ui


#endif // TASK_H
