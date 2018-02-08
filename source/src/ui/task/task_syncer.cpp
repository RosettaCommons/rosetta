#include <ui/task/task_syncer.h>

#include <ui/task/task.h>
#include <ui/task/functor.h>

#include <ui/util/serialization.h>

#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
//#include <QCryptographicHash>
#include <QTimer>

namespace ui {
namespace task {


TaskSyncer_NodeStrategy::TaskSyncer_NodeStrategy(Task * task) : task_(task)
{
}


void TaskSyncer_NodeStrategy::submit(QString const &queue)
{
	task_->queue_ = queue;
	task_->state_ = Task::State::_draft_;

	create_sync_tree();
	qDebug() << "TaskSyncer_NodeStrategy[" << task_->task_id() << "]: Task::submit() Name:" << task_->name_ << " queue:" << task_->queue_;

	connect(root_.get(), SIGNAL(tree_synced()), this, SLOT(draft_to_queued()));

	root_->data_is_fresh(true);

	Q_EMIT task_->changed();
}

void TaskSyncer_NodeStrategy::draft_to_queued(void)
{
	if( root_ )  {
		qDebug() << "TaskSyncer_NodeStrategy[" << task_->task_id() << "]: draft_to_queued()...";
		disconnect(root_.get(), SIGNAL(tree_synced()), this, SLOT(draft_to_queued()) );

		task_->state_ = Task::State::_queued_;

		connect(root_.get(), SIGNAL(synced()), this, SLOT(post_submit()));
		root_->data_is_fresh(false);

		if( task_->project_ and !task_->project_->file_name().isEmpty() ) save_project(*task_->project_, /* always_ask_for_file_name = */ false);
	}
}

void TaskSyncer_NodeStrategy::post_submit(void)
{
	qDebug() << "TaskSyncer_NodeStrategy[" << task_->task_id() << "]: post_submit()...";
	disconnect(root_.get(), SIGNAL(synced()), this, SLOT(post_submit()));
	Q_EMIT task_->submitted();
	Q_EMIT task_->changed();
	subscribe();
}


void TaskSyncer_NodeStrategy::subscribe()
{
	if(root_) {
		qDebug() << "TaskSyncer_NodeStrategy[" << task_->task_id() << "]: Name:" << task_->name_ << " TaskSyncer_NodeStrategy::subscribe";
		Updater::get().subscribe( root_.get() );
	}
}

void TaskSyncer_NodeStrategy::create_sync_tree()
{
	//QUuid task_id = QUuid::createUuid();  /// each distinct submitted Task _must_ have unique id

	root_ = std::make_shared<Node>("task", Node::Flags::data_in | Node::Flags::data_out | Node::Flags::topology_out/*, task_id*/);

	auto files_node = std::make_shared<Node>(Node::Flags::topology_in | Node::Flags::topology_out);
	root_->add("files", files_node);

	// auto input_node = std::make_shared<Node>(Node::Flags::data_out | Node::Flags::topology_out);
	// root_->add("input", input_node);
	// //assign_input_node();

	// auto script_node = std::make_shared<Node>(Node::Flags::data_out | Node::Flags::topology_out);
	// root_->add("script", script_node);
	// //assign_script_node();  //script_node_ = script_node;

	// auto flags_node = std::make_shared<Node>(Node::Flags::data_out | Node::Flags::topology_out);
	// root_->add("flags", flags_node);
	// //assign_flags_node();

	connect_task_and_nodes();
}

void TaskSyncer_NodeStrategy::connect_task_and_nodes()
{
	if(root_) {
		TaskSyncer_NodeStrategyQP qsyncer(this);

		root_->get_data_callback( [qsyncer]() { if(qsyncer) return qsyncer->task_->task_data(); else return QJsonValue(); } );
		root_->set_data_callback( [qsyncer](QJsonValue const &ba) { if(qsyncer) return qsyncer->task_->task_data(ba); } );

		files_node_ = root_->leaf("files");
		if(auto files_node = files_node_.lock() ) {

			if(task_->state_ == Task::State::_draft_) {
				for(auto & file : task_->files_) {
					auto fnode = std::make_shared<Node>(Node::Flags::data_out | Node::Flags::topology_out);
					files_node->add(file.first, fnode);
					auto file_sp = file.second;
					fnode->get_data_callback( [file_sp]() -> QJsonValue { return QJsonValue( QLatin1String( file_sp->data().toBase64() ) ); } );
				}
			}

			files_node->topology_updated_callback(
				[qsyncer](Node const *node, std::vector<QString> const & new_keys, std::vector<QString> const & errased_keys) -> void
				{ if(qsyncer) qsyncer->files_topology_updated(node, new_keys, errased_keys); }
		    );
			// connect(output_node.get(), SIGNAL(        topology_updated(Node const *, std::vector<QString> const & , std::vector<QString> const & ) ),
			// 		this,                SLOT( output_topology_updated(Node const *, std::vector<QString> const & , std::vector<QString> const & ) ));
		}


		// input_node_ = root_->leaf("input");
		// if(auto input_node = input_node_.lock() ) input_node->get_data_callback( [qsyncer]() -> QJsonValue { if(qsyncer) return QJsonValue( QLatin1String( qsyncer->input_.data().toBase64() ) ); else return QJsonValue(); } );

		// script_node_ = root_->leaf("script");
		// if(auto script_node = script_node_.lock() )  script_node->get_data_callback( [qsyncer]() -> QJsonValue { if(qsyncer) return QJsonValue( QLatin1String( qsyncer->script_.data().toBase64() ) ); else return QJsonValue(); } );

		// flags_node_ = root_->leaf("flags");
		// if(auto flags_node = flags_node_.lock() ) flags_node->get_data_callback( [qsyncer]() -> QJsonValue { if(qsyncer) return QJsonValue( QLatin1String( qsyncer->flags_.data().toBase64() ) ); else return QJsonValue(); } );

		// output_node_ = root_->leaf("output");
		// if(auto output_node = output_node_.lock() ) {
		// 	output_node->topology_updated_callback(
		// 		[qsyncer](Node const *node, std::vector<QString> const & new_keys, std::vector<QString> const & errased_keys) -> void
		// 		{ if(qsyncer) qsyncer->output_topology_updated(node, new_keys, errased_keys); }
		//     );
		// 	// connect(output_node.get(), SIGNAL(        topology_updated(Node const *, std::vector<QString> const & , std::vector<QString> const & ) ),
		// 	// 		this,                SLOT( output_topology_updated(Node const *, std::vector<QString> const & , std::vector<QString> const & ) ));
		// }

		connect(root_.get(), SIGNAL(tree_synced()), task_, SIGNAL(changed()));
		connect(root_.get(), SIGNAL(syncing()),     task_, SIGNAL(syncing()));

	} else {
		files_node_.reset();
	}
}

void TaskSyncer_NodeStrategy::files_topology_updated(Node const *, std::vector<QString> const & new_keys, std::vector<QString> const & errased_keys)
{
	qDebug() << "TaskSyncer_NodeStrategy::output_topology_updated: New keys: " << new_keys << " Errased keys:" << errased_keys;

	if(auto files_node = files_node_.lock() ) {
		TaskSyncer_NodeStrategyQP qsyncer(this);

		for(auto const &k : errased_keys) task_->files_.erase(k);

		for(auto const &k : new_keys) {
			if(NodeSP leaf = files_node->leaf(k) ) {
				FileSP file_sp = std::make_shared<File>();
				file_sp->file_name(k);
				task_->files_[k] = file_sp;

				qDebug() << "TaskSyncer_NodeStrategy::output_topology_updated: setting set-data-callback for node: " << leaf->node_id();
				leaf->set_data_callback( [qsyncer, file_sp](QJsonValue const &jv) {
						//qDebug() << "Data update for output key(jv):" << jv;
						auto ba = QByteArray::fromBase64( jv.toVariant().toByteArray() );

						qDebug() << "Data update for output key:" << file_sp->file_name() << "  data: " << ba.left(64) << ( ba.size() > 64 ? "..." : "");
						if(qsyncer) {
							file_sp->data(ba);
							Q_EMIT qsyncer->task_->changed();
						}
					} );

			}
		}
		if( not errased_keys.empty()  or  not new_keys.empty() ) task_->files_model_.update_from_task(*task_);
	}
}


bool TaskSyncer_NodeStrategy::is_syncing() const
{
	if(root_) return root_->is_syncing(true);
	return false;
}

// return <current_value, max_value> for syncing operation progress
std::pair<int, int> TaskSyncer_NodeStrategy::syncing_progress() const
{
	if(root_) {
		auto res = std::make_pair(root_->is_syncing(true), root_->tree_size());
		qDebug() << "TaskSyncer_NodeStrategy::syncing_progress():" << res;
		return res;
	}
	else return std::make_pair(0, 0);
}


quint32 const _TaskSyncer_NodeStrategy_stream_version_ = 0x00000001;

QDataStream &operator<<(QDataStream &out, TaskSyncer_NodeStrategy const&t)
{
	out.setVersion(QDataStream::Qt_5_6);

	out << ui::util::_TaskSyncer_NodeStrategy_magic_number_;
	out << _TaskSyncer_NodeStrategy_stream_version_;

	out << (t.root_ != nullptr);
    if( t.root_ != nullptr ) out << *t.root_;

	return out;
}


QDataStream &operator>>(QDataStream &in, TaskSyncer_NodeStrategy &t)
{
	in.setVersion(QDataStream::Qt_5_6);

	quint64 magic;
	in >> magic;
	if( magic != ui::util::_TaskSyncer_NodeStrategy_magic_number_ ) throw ui::util::BadFileFormatException( QString("Invalid _TaskSyncer_NodeStrategy_magic_number_: read %1, was expecting %2...").arg(magic).arg(ui::util::_TaskSyncer_NodeStrategy_magic_number_) );

	quint32 version;
	in >> version;
	if( version != _TaskSyncer_NodeStrategy_stream_version_ ) throw ui::util::BadFileFormatException( QString("Invalid _TaskSyncer_NodeStrategy_stream_version_: read %1, was expecting %2...").arg(magic).arg(_TaskSyncer_NodeStrategy_stream_version_) );

	bool root;
	in >> root;
	if(root) {
		t.root_ = std::make_shared<Node>();
		in >> *t.root_;

		t.connect_task_and_nodes();

		if( t.task_->state() == Task::State::_draft_ ) { // task was submitted but not yet synced, - trying to continue...
			qDebug() << "QDataStream &operator>>(QDataStream &in, TaskSyncer_NodeStrategy &t): task was submitted but not yet synced, - trying to continue...";
			//QObject::connect(t.root_.get(), SIGNAL(tree_synced()), &t, SLOT(draft_to_queued()));
			//t.root_->data_is_fresh(true);
			t.submit( t.task_->queue() );
		}
		else if( t.task_->state() != Task::State::_none_ ) t.subscribe();
	}

	return in;
}

TaskSyncer_TaskStrategy::TaskSyncer_TaskStrategy(Task * task) : task_(task)
{
	//connect(task_, &Task::submitted, this, &TaskSyncer_TaskStrategy::submitted);
}


void TaskSyncer_TaskStrategy::submit(QString const &queue)
{
	task_->queue_ = queue;
	task_->state_ = Task::State::_draft_;

	qDebug() << "TaskSyncer_TaskStrategy[" << task_->task_id() << "]: Task::submit() Name:" << task_->name_ << " queue:" << task_->queue_;

	task_data_upload();

	Q_EMIT task_->changed();
}

void TaskSyncer_TaskStrategy::task_data_upload()
{
	auto f = std::make_shared<FunctorNetworkCall>("TaskDataUploader for '" + task_->name() + "'");

	QPointer<Task> qtask = task_;
	QPointer<FunctorNetworkCall> qf = f.get();

	f->callback(
		[qtask](NetworkCall &nc) {
			if(qtask) {
				qtask->syncing_progress(0, 1);

				qDebug() << "TaskSyncer_TaskStrategy, initiating NetworkCall...";

				QJsonValue jvalue = qtask->task_data();
				QJsonObject jobject = jvalue.toObject();
				QJsonDocument jd = QJsonDocument(jobject); //QJsonDocument::fromVariant( QVariant(jv) );
				//qDebug() << "TaskSyncer_TaskStrategy: jobject:" << jobject << " jd:" << jd;
				nc.call(task_api_url() + "/task", QNetworkAccessManager::Operation::PostOperation, jd );
			} else {
				qDebug() << "TaskSyncer_TaskStrategy, TASK IS GONE!!!";
				// what to do??? Q_EMIT f->finish();
			}
		});
		// [](NetworkCall &nc) {
		// 	qDebug() << "TaskSyncer_TaskStrategy, initiating NetworkCall...";
		// 	nc.call(task_api_url() + "/task", QNetworkAccessManager::Operation::PostOperation, QJsonDocument() );
		// });

	connect(f.get(), &Functor::finished, [qf, qtask]() {
			if(qf  and  qtask) {
				QJsonDocument jd = qf->result_as_json();
				QJsonObject root = jd.object();
				qtask->task_data(root);

				// if( root["task_id"].isDouble() ) {
				// 	root["task_id"].toDouble();
				// }

				qtask->syncing_progress(0, 0);
				if( qtask->project_ ) save_project(*qtask->project_, /* always_ask_for_file_name = */ false);
			}
		});

	//connect(f.get(), &Functor::tick,     task_, &Task::syncing);
	//connect(f.get(), &Functor::finished, task_, &Task::syncing);

	connect(f.get(), &Functor::final, this, &TaskSyncer_TaskStrategy::task_files_upload);
	functor_ = f;
	functor_->execute();
}

void TaskSyncer_TaskStrategy::task_files_upload()
{
	QString task_id = task_->task_id();

	auto f = std::make_shared<FunctorSequence>("TaskFilesUploader for '" + task_->name() + "' task_id=" + task_id);

	QPointer<Task> qtask = task_;

	for(auto const & fl : task_->files() ) {
		QString file_name = fl.first;
		FileSP file = fl.second;

		auto u = std::make_shared<FunctorNetworkCall>("File uploader for '" + task_->name() + "' task_id=" + task_id + " file:" + file_name);
		QPointer<FunctorNetworkCall> qu = u.get();

		u->callback(
			[task_id, file_name, file](NetworkCall &nc) {
				nc.call(task_api_url() + "/file/" + task_id + "/" + file_name, QNetworkAccessManager::Operation::PostOperation, file->data());
			});
		f->add_functor(u);

		connect(u.get(), &Functor::finished,
				[qtask, qu, file]() {
					if(qu) {
						QJsonDocument jd = qu->result_as_json();
						QJsonObject o = jd.object();

						auto ho = o["hash"];
						if( ho.isString() ) {
							file->hash( ho.toString() );
						}
						qtask->syncing_progress_advance(1);
					}
				});
	}

	//connect(f.get(), &Functor::tick,  task_, &Task::syncing);
	//connect(f.get(), &Functor::final, task_, &Task::syncing);

	connect(f.get(), &Functor::final, this, &TaskSyncer_TaskStrategy::task_queuing);

	qtask->syncing_progress(0, task_->files().size());

	functor_ = f;
	functor_->execute();
}

void TaskSyncer_TaskStrategy::task_queuing()
{
	task_->state_ = Task::State::_queued_;

	QString task_id = task_->task_id();
	auto f = std::make_shared<FunctorNetworkCall>("TaskQueueing for '" + task_->name() + "' task_id=" + task_id);

	QPointer<Task> qtask = task_;
	QPointer<FunctorNetworkCall> qf = f.get();

	f->callback(
		[qtask, task_id](NetworkCall &nc) {
			if(qtask) {
				qtask->syncing_progress(0, 1);

				QJsonValue jvalue = qtask->task_data();
				QJsonObject jobject = jvalue.toObject();
				QJsonDocument jd = QJsonDocument(jobject);
				nc.call(task_api_url() + "/task/" + task_id, QNetworkAccessManager::Operation::PostOperation, jd );
			} else {
				qDebug() << "TaskSyncer_TaskStrategy, TASK IS GONE!!!";
			}
		});

	connect(f.get(), &Functor::finished, [qf, qtask]() {
			if(qf  and  qtask) {
				QJsonDocument jd = qf->result_as_json();
				QJsonObject root = jd.object();
				qtask->task_data(root);

				qtask->syncing_progress(0, 0);
				if( qtask->project_ ) save_project(*qtask->project_, /* always_ask_for_file_name = */ false);
			}
		});

	//connect(f.get(), &Functor::tick,  task_, &Task::syncing);
	//connect(f.get(), &Functor::final, task_, &Task::syncing);

	connect(f.get(), &Functor::final, task_, &Task::submitted);
	connect(f.get(), &Functor::final, task_, &Task::subscribe);
	functor_ = f;
	functor_->execute();

	Q_EMIT task_->changed();
}


void TaskSyncer_TaskStrategy::subscribe()
{
	QPointer<Task> qtask = task_;
	QString task_id = task_->task_id();

	struct TaskDiff {
		QJsonObject root;
	};

	auto diff = std::make_shared<TaskDiff>();

	auto diff_functor = std::make_shared<FunctorNetworkCall>("TaskUpdateFunctor::run(): '" + task_->name() + "' task_id=" + task_->task_id());
	QPointer<FunctorNetworkCall> diff_functor_qp = diff_functor.get();

	diff_functor->callback(
		[qtask, task_id](NetworkCall &nc) {
			if(qtask) {
				qDebug() << "TaskUpdateFunctor::run, initiating NetworkCall for getting diff...";

				qtask->syncing_progress(0, 1);

				QJsonObject o, files;

				for(auto const & f : qtask->files() ) files[f.first] = f.second->hash();

				o["files"] = files;

				//nc.call(task_api_url() + "/task_diff/" + task_id, QNetworkAccessManager::Operation::PostOperation, QJsonDocument(o));

				QByteArray payload( QJsonDocument(o).toJson(QJsonDocument::Compact) );
				payload = qCompress(payload).mid(4); // removing extra 32bit int (size of uncrompressed data from beginning of the data to make result compatible with raw zlib encoding

				//nc.call(task_api_url() + "/task_diff/" + task_id, QNetworkAccessManager::Operation::PostOperation, payload, "application/json; charset=utf-8", { {"Content-Encoding", "deflate"} });
				nc.call(task_api_url() + "/task_diff/" + task_id, QNetworkAccessManager::Operation::PostOperation, payload, "application/x-zlib-json");
			}
		});

	connect(diff_functor_qp, &Functor::finished,
			[qtask, diff_functor_qp, diff]() {
				if(diff_functor_qp and qtask) {
					diff->root = diff_functor_qp->result_as_json().object();

					QJsonObject files = diff->root["files"].toObject();
					QJsonArray deleted = files["deleted"].toArray();
					for(int i=0; i < deleted.size(); ++i) {
						QString file_name = deleted.at(i).toString();
						qtask->delete_file(file_name);
					}
				}
			});


	auto file_downloader_functor = std::make_shared<FunctorSequence>("Task-sequence-file-downloader<'" + task_->name() + "' task_id=" + task_id+">");
	QPointer<FunctorSequence> file_downloader_functor_qp = file_downloader_functor.get();

	connect(file_downloader_functor_qp, &Functor::started,
			[qtask, task_id, file_downloader_functor_qp, diff]() {
				if(file_downloader_functor_qp and qtask) {
					QJsonObject files = diff->root["files"].toObject();
					QJsonObject updated = files["updated"].toObject();
					QStringList file_names = updated.keys();

					qtask->syncing_progress(0, file_names.size());

					for (int i = 0; i < file_names.size(); ++i) {
						QString file_name = file_names.at(i);
						QString file_hash = updated[file_name].toString();

						//qDebug() << QString("Adding file %1 to download queue").arg( file_names.at(i) );
						auto d = std::make_shared<FunctorNetworkCall>("File downloader for '" + qtask->name() + "' task_id=" + task_id + " file:" + file_name);
						QPointer<FunctorNetworkCall> qd = d.get();

						d->callback(
							[task_id, file_name](NetworkCall &nc) {
								nc.set_termination_codes( {404} );
								nc.call(task_api_url() + "/file/" + task_id + "/" + file_name, QNetworkAccessManager::Operation::GetOperation);
							});
						file_downloader_functor_qp->add_functor(d);

						connect(d.get(), &Functor::finished,
								[qd, file_name, file_hash, qtask]() {
									if(qd and qtask) {
										if( qd->status_code() == 200 ) {

											QByteArray data = qd->result();

											// should we re-calculate hash locally instead? //QByteArray hash = QCryptographicHash::hash(data, QCryptographicHash::Md5);

											auto fl = std::make_shared<File>(data); fl->hash(file_hash);
											qtask->add_file(file_name, fl);
											Q_EMIT qtask->changed();
											//Q_EMIT qtask->syncing();

											qDebug() << QString("Adding file %1 to task:%2").arg(file_name).arg( qtask->task_id() );
										} else {
											qDebug() << QString("Could not retrive file %1... server replied: %2").arg(file_name).arg( qd->status_code() );
										}

										qtask->syncing_progress_advance(1);
									}
								});
					}
				}
			});


	connect(file_downloader_functor_qp, &Functor::final,
			[qtask, diff]() {
				if(qtask) {
					qtask->syncing_progress(0, 0);

					qtask->task_data(diff->root);
					if( qtask->state() == Task::State::_finished_) {
						Q_EMIT qtask->final();
					}
					else {
						if( qtask->state() == Task::State::_running_) QTimer::singleShot(1000*60*1, qtask, &Task::subscribe );
						else  QTimer::singleShot(1000*60*1, qtask, &Task::subscribe );
					}
					//Q_EMIT qtask->syncing();

					qDebug() << "Syncing state:" << qtask->is_syncing() << " progress:" << qtask->syncing_progress();
				}
			});


	auto diff_then_download = std::make_shared<FunctorSequence>("Task-diff-then-download<'" + task_->name() + "' task_id=" + task_id+">");

	diff_then_download->add_functor(diff_functor);
	diff_then_download->add_functor(file_downloader_functor);

	functor_ = diff_then_download;

	//connect(functor_.get(), &Functor::tick,  task_, &Task::syncing);
	//connect(functor_.get(), &Functor::final, task_, &Task::syncing);

	functor_->execute();
}

bool TaskSyncer_TaskStrategy::is_syncing() const
{
	//	return functor_  and  functor_->is_running();

	if(functor_) {
		auto p = functor_->progress();
		return p.first != p.second;
	}
	else return false;
}

// return <current_value, max_value> for syncing operation progress
std::pair<int, int> TaskSyncer_TaskStrategy::syncing_progress() const
{
	if(functor_) return functor_->progress();
	else return std::make_pair(1, 1);
}


void TaskSyncer_TaskStrategy::resume()
{
	if     ( task_->state() == Task::State::_none_ ) { // do nothing, Task was not yet submitted
	}
	else if( task_->state() == Task::State::_draft_ ) {
		if( task_->task_id().isEmpty() ) task_data_upload();
		else task_files_upload();
	}
	else if( task_->state() == Task::State::_queued_ ) {
		task_queuing(); // re-queueing in case state was changed only locally. The front-end handle will prevent Task from going back to `queue` state if it is already running
	}
	else if( task_->state() == Task::State::_running_ ) {
		subscribe();
	}
	else if( task_->state() == Task::State::_finished_ ) { // do nothing if Task could became finished only if it is already fully synced
	}
}


quint32 const _TaskSyncer_TaskStrategy_stream_version_ = 0x00000001;

QDataStream &operator<<(QDataStream &out, TaskSyncer_TaskStrategy const&)
{
	out.setVersion(QDataStream::Qt_5_6);

	out << ui::util::_TaskSyncer_TaskStrategy_magic_number_;
	out << _TaskSyncer_TaskStrategy_stream_version_;

	//out << (t.root_ != nullptr);
    //if( t.root_ != nullptr ) out << *t.root_;

	return out;
}


QDataStream &operator>>(QDataStream &in, TaskSyncer_TaskStrategy &t)
{
	in.setVersion(QDataStream::Qt_5_6);

	quint64 magic;
	in >> magic;
	if( magic != ui::util::_TaskSyncer_TaskStrategy_magic_number_ ) throw ui::util::BadFileFormatException( QString("Invalid _TaskSyncer_TaskStrategy_magic_number_: read %1, was expecting %2...").arg(magic).arg(ui::util::_TaskSyncer_TaskStrategy_magic_number_) );

	quint32 version;
	in >> version;
	if( version != _TaskSyncer_TaskStrategy_stream_version_ ) throw ui::util::BadFileFormatException( QString("Invalid _TaskSyncer_TaskStrategy_stream_version_: read %1, was expecting %2...").arg(magic).arg(_TaskSyncer_TaskStrategy_stream_version_) );

	t.resume();

	return in;
}

/*
TaskUpdateFunctor::TaskUpdateFunctor(Task *task, QObject *parent) : Functor( QString("TaskUpdateFunctor<task_id=%1>").arg( task->task_id() ), parent), task_(task)
{}

void TaskUpdateFunctor::run()
{
	QString task_id = task_->task_id();
	auto f = std::make_shared<FunctorNetworkCall>("TaskUpdateFunctor::run(): '" + task_->name() + "' task_id=" + task_->task_id());

	QPointer<Task> qtask = task_;
	QPointer<FunctorNetworkCall> qf = f.get();

	f->callback(
		[qtask, task_id](NetworkCall &nc) {
			if(qtask) {
				qDebug() << "TaskUpdateFunctor::run, initiating NetworkCall for getting diff...";

				QJsonObject o, files;

				for(auto const & f : qtask->files() ) files[f.first] = f.second->hash();

				o["files"] = files;

				nc.call(task_api_url() + "/task_diff/" + task_id, QNetworkAccessManager::Operation::PostOperation, QJsonDocument(o));
			}
		});

	connect(f.get(), &Functor::finished,
			[this, qf, qtask, task_id]() {
				if(qf  and  qtask) {
					this->diff_ = qf->result_as_json();

					QJsonObject update_root = this->diff_.object();
					QJsonObject files = update_root["files"].toObject();
					QJsonObject updated = files["updated"].toObject();
					QStringList file_names = updated.keys();

					auto fd = std::make_shared<FunctorSequence>("TaskFilesDownloader for '" + task_->name() + "' task_id=" + task_id);

					for (int i = 0; i < file_names.size(); ++i) {
						QString file_name = file_names.at(i);
						QString file_hash = updated[file_name].toString();

						//qDebug() << QString("Adding file %1 to download queue").arg( file_names.at(i) );
						auto d = std::make_shared<FunctorNetworkCall>("File uploader for '" + task_->name() + "' task_id=" + task_id + " file:" + file_name);
						QPointer<FunctorNetworkCall> qd = d.get();

						d->callback(
							[task_id, file_name](NetworkCall &nc) {
								nc.call(task_api_url() + "/file/" + task_id + "/" + file_name, QNetworkAccessManager::Operation::GetOperation);
							});
						fd->add_functor(d);

						connect(d.get(), &Functor::finished,
								[qd, file_name, file_hash, qtask]() {
									if(qd and qtask) {
										QByteArray data = qd->result();

										// should we re-calculate hash locally instead? //QByteArray hash = QCryptographicHash::hash(data, QCryptographicHash::Md5);

										auto fl = std::make_shared<File>(data); fl->hash(file_hash);
										qtask->add_file(file_name, fl);
										qDebug() << QString("Adding file %1 to task:%2").arg(file_name).arg( qtask->task_id() );
									}
								});
					}

					connect(fd.get(), &Functor::final,
							[qtask, update_root]() {
								if(qtask) {
									qtask->task_data(update_root);
									if( qtask->state() == Task::State::_finished_) {
										Q_EMIT qtask->final();
									}
									else QTimer::singleShot(1000*60, qtask, &Task::subscribe );
								}
							});

					//connect(fd.get(), &Functor::final, qtask, &Task::subscribe);
					file_downloader_functor_ = fd;
					file_downloader_functor_->execute();
				}
			});
	//connect(f.get(), &Functor::final, this, &TaskSyncer_TaskStrategy::task_files_upload);
	diff_functor_ = f;
	diff_functor_->execute();

}
*/

} // namespace task
} // namespace ui
