#ifndef TASK_SYNCER_H
#define TASK_SYNCER_H

#include <ui/task/task_syncer.fwd.h>

#include <ui/task/task.fwd.h>
#include <ui/task/node.h>
#include <ui/task/functor.h>
#include <ui/task/util.h>

#include <QTimer>

namespace ui {
namespace task {


class TaskSyncer_NodeStrategy : public QObject
{
    Q_OBJECT
public:
    //explicit TaskSyncer_NodeStrategy(QObject *parent = nullptr);
	explicit TaskSyncer_NodeStrategy(Task *);

	/// Initiate submit procedure
	///   0. set state to _draft_
	///   1. sync all nodes to cloud
	///   2. set state to _queued_
	///   3. sync root node to cloud
	void submit(QString const &queue);

	void subscribe();

	bool is_syncing() const;

	// return <current_value, max_value> for syncing operation progress
	std::pair<int, int> syncing_progress() const;

	// serialization
	friend QDataStream &operator<<(QDataStream &, TaskSyncer_NodeStrategy const&);
	friend QDataStream &operator>>(QDataStream &, TaskSyncer_NodeStrategy &);

Q_SIGNALS:
	void submitted();

public Q_SLOTS:


private Q_SLOTS:
	void draft_to_queued(void);

	void post_submit(void);

private:
	void create_sync_tree();
	void connect_task_and_nodes();
	void files_topology_updated(Node const *, std::vector<QString> const & new_keys, std::vector<QString> const & errased_keys);

	/// root node of network syncing tree
	NodeSP root_;

	NodeWP files_node_;

	Task * task_;
};



class TaskSyncer_TaskStrategy : public QObject
{
    Q_OBJECT
public:
    //explicit TaskSyncer_NodeStrategy(QObject *parent = nullptr);
	explicit TaskSyncer_TaskStrategy(Task *);

	/// Initiate submit procedure
	///   0. set state to _draft_
	///   1. sync all nodes to cloud
	///   2. set state to _queued_
	///   3. sync root node to cloud
	void submit(QString const &queue);

	void subscribe();

	bool is_syncing() const;

	// return <current_value, max_value> for syncing operation progress
	std::pair<int, int> syncing_progress() const;

	// serialization
	friend QDataStream &operator<<(QDataStream &, TaskSyncer_TaskStrategy const&);
	friend QDataStream &operator>>(QDataStream &, TaskSyncer_TaskStrategy &);

//Q_SIGNALS:
// 	void submitted();

private Q_SLOTS:
	//void update();

private:
	/// initiate task data upload and when its done initiate files upload
	void task_data_upload();

	/// initiate task files upload and when its done initiate files upload, then call task_queuing
	void task_files_upload();

	// Set state to `queued` and upload task data. When complete emit submitted
	void task_queuing();

	/// initiate files upload
	void files_upload();

	// serialization helper: check state of Task and resume syncing operations if any
	void resume();

private:
	FunctorSP functor_;

	Task * task_;

	//QTimer *timer_;
	//bool file_list_changed_ = false;
};


/*
/// Functor to pefromr Task update
class TaskUpdateFunctor : public Functor
{
    Q_OBJECT
public:
	TaskUpdateFunctor(Task *, QObject *parent = Q_NULLPTR);

private Q_SLOTS:
	void run() override;

private:
	QJsonDocument diff_;

	Task * task_;
	FunctorSP diff_functor_, file_downloader_functor_;
};
*/

} // namespace task
} // namespace ui

#endif // TASK_SYNCER_H
