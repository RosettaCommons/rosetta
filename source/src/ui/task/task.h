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

//#include <ui/task/node.h>
#include <ui/task/project.h>
#include <ui/task/file.h>

#include <ui/util/exception.h>

//#include <QObject>
//#include <QJsonObject>

#include <ui/task/task_syncer.h>

namespace ui {
namespace task {


class Task : public QObject
{
	Q_OBJECT

public:
	enum class State {_none_, _draft_, _queued_, _running_, _finished_, _unknown_};

	struct Job {
		QString name, app, nstruct, flags;
	};
	using JobSP = std::shared_ptr<Job>;

public:
	explicit Task();
	//explicit Task(QString const &description);
	~Task();

	/// create clone of Task
	TaskSP clone() const;

	// return current Task State
	State state() const { return state_; }

	// return queue that will be used to run this Task
	QString queue() const { return queue_; }

	// return true if there is any outgoing/pending network operations
	bool is_syncing() const;

	// return <current_value, max_value> for syncing operation progress
	std::pair<int, int> syncing_progress() const;

	QAbstractTableModel * files_model() { return &files_model_; }

	Project *project() const { return project_; }

	/// add a new job at the end of jobs queue and return poiinter to it
	JobSP add_job();

	/// return jobs vector
	std::vector<JobSP> const & jobs() const { return jobs_; }

	/// delete given job
	void delete_job(JobSP const &job);

	/// try to rename job, return true on success
	bool rename_job(JobSP const&, QString const & name);

	/// swap jobs of given indexes
	void swap_jobs(int i, int j);

	/// add file to Task files collection
	void add_file(QString const &name, FileSP const &file);
	std::map<QString, FileSP> const & files() const { return files_; }

	/// delete file from task, return true if file was in task files
	bool delete_file(QString const &name);

	// should not be needed, add/remove Tasks should be done on Project level
	// void project(Project *p)  { project_ = p; }

	/// Initiate submit procedure
	void submit(QString const & queue);

	/// return QString representation of task_id (int) if task was already submitted other wise return empty string
	QString task_id() const;

	QString name() const { return name_; }
	void name(QString const &);

	QString version() const { return version_; }
	void version(QString const &);

	QString description() const { return description_; }
	void description(QString const &);

	// File const &input() const { return input_; }
	// void input(File &&input) { if( input_ != input ) { input_ = std::move(input);  Q_EMIT changed(); } }
	// File const &script() const { return script_; }
	// void script(File &&script) { if( script_ != script ) { script_ = std::move(script); Q_EMIT changed(); } }
	// File const &flags() const { return flags_; }
	// void flags(File &&flags) { if( flags_ != flags ) { flags_ = std::move(flags); Q_EMIT changed(); } }

	/// subscribe to network update stream
	void subscribe();

	// bool operator ==(Task const &r) const;
	// bool operator !=(Task const &r) const { return not (*this == r); }

	// serialization
	friend QDataStream &operator<<(QDataStream &, Task const&);
	friend QDataStream &operator>>(QDataStream &, Task &);


	static QString to_string(Task::State state);
	static Task::State from_string(QString const &s);

public Q_SLOTS:
	void rename_file(QString const &previous_value, QString const &new_value);

Q_SIGNALS:
	void submitted();

    void changed();
    void name_changed();
    void state_changed();

	void file_list_changed();
	void file_changed(QString const &);

	/// Emitted when node or sub-node syncing state changed
	void syncing();


	/// Emitted once when Task is fully downloaded and no further network operations will be issued
	void final();

private Q_SLOTS:


	//void output_topology_updated(Node const *, std::vector<QString> const & new_keys, std::vector<QString> const &  errased_keys);
	//void output_topology_updated(Node const *);


private:

	/// connect nodes and Task structre, assign callback's
	//void connect_task_and_nodes();

	//void output_topology_updated(Node const *, std::vector<QString> const & new_keys, std::vector<QString> const &  errased_keys);

private:
	QJsonValue task_data();
	void task_data(QJsonValue const &);

	// QJsonObject task_data();
	// void task_data(QJsonObject const &);

	// set <current_value, max_value> for syncing operation progress
	void syncing_progress(int value, int max);

	// move progress forward by `value` points
	void syncing_progress_advance(int value);

	/// Some project functions require access to project_
	friend void Project::add_task(TaskSP const &task);
	friend bool Project::erase(TaskSP const &task);
	friend void Project::assign_ownership(TaskSP const &t);

	friend class TaskSyncer_NodeStrategy;
	friend class TaskSyncer_TaskStrategy;
	//friend class TaskUpdateFunctor;

	State state_ = State::_none_;

	QString queue_;

	QString version_;

	QString name_, description_;

	QVariant task_id_;

	//File input_, script_, flags_;
	//std::map<QString, FileSP> output_;

	std::vector<JobSP> jobs_;

	std::map<QString, FileSP> files_;

	FileTableModel files_model_;

	//TaskSyncer_NodeStrategy syncer_;
	TaskSyncer_TaskStrategy syncer_;

	//NodeWP input_node_, script_node_, flags_node_;
	//NodeWP output_node_;

	/// root node of network syncing tree
	//NodeSP root_;

	QPointer<Project> project_;

	std::pair<int, int> syncing_progress_ = {0, 0};
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
