// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/task/task.cpp
/// @brief  Task class for ui library.
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu).

#include <ui/task/task.h>
#include <ui/task/util.h>
#include <ui/util/serialization.h>

#include <QJsonDocument>
#include <QJsonObject>
//#include <QVariant>
#include <QDataStream>
#include <QFile>
#include <QDir>
#include <QTimer>

#include <iterator>

namespace ui {
namespace task {

Node::Key const _input_key_ {"input"};
Node::Key const _output_key_{"output"};

Node::Key const _script_key_ {"script"};

File::File(QString const &file_name)
{
	init_from_file(file_name);
}

File::File(QByteArray const & file_data)
{
	data(file_data);
}


File::File(QString const & file_name, QByteArray const & file_data) : file_name_(file_name)
{
	data(file_data);
}


File& File::operator=(File&& other) noexcept
{
    if(this != &other) { // no-op on self-move-assignment
		file_name_.swap(other.file_name_);
		file_data_.swap(other.file_data_);
    }
    return *this;
}


void File::init_from_file(QString const & file_name)
{
	file_name_ = file_name;
	QFile file(file_name_);
    if( file.open(QIODevice::ReadOnly) ) {
		data( file.readAll() );
	}
}

QByteArray File::data() const
{
	return qUncompress(file_data_);
}

void File::data(QByteArray const &file_data)
{
	file_data_ = qCompress(file_data);
}


bool File::operator ==(File const &rhs) const
{
	return (file_name_ == rhs.file_name_) and (file_data_ == rhs.file_data_);
}

QDataStream &operator<<(QDataStream &out, File const&f)
{
	out << f.file_name_;
	out << f.file_data_;
	return out;
}

QDataStream &operator>>(QDataStream &in, File &f)
{
	in >> f.file_name_;
	in >> f.file_data_;
	return in;
}






/*
NodeTaskSyncer::NodeTaskSyncer(Task *task) : task_(task)
{
}

void NodeTaskSyncer::submit(QString const & queue)
{
	queue_ = queue;
	state_ = State::_draft_;

	create_sync_tree();
	qDebug() << "Task[" << task_id() << "]: Task::submit() Name:" << (project_ ? *project_->find(this) : "''") << " queue:" << queue_;

	connect(root_.get(), SIGNAL(tree_synced()), this, SLOT(draft_to_queued()));

	root_->data_is_fresh(true);

	Q_EMIT task_->changed();
}



void NodeTaskSyncer::create_sync_tree()
{
	//QUuid task_id = QUuid::createUuid();  /// each distinct submitted Task _must_ have unique id

	root_ = std::make_shared<Node>("task", Node::Flags::data_in | Node::Flags::data_out | Node::Flags::topology_out);

	auto files_node = std::make_shared<Node>(Node::Flags::all);
	root_->add("files", files_node);
	//assign_input_node();




	auto script_node = std::make_shared<Node>(Node::Flags::data_out | Node::Flags::topology_out);
	root_->add("script", script_node);
	//assign_script_node();  //script_node_ = script_node;

	auto flags_node = std::make_shared<Node>(Node::Flags::data_out | Node::Flags::topology_out);
	root_->add("flags", flags_node);
	//assign_flags_node();

	auto output_node = std::make_shared<Node>(Node::Flags::topology_in);
	root_->add("output", output_node);
	//assign_output_node();

	connect_task_and_nodes();
}
*/



QVariant FileTableModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (role == Qt::DisplayRole) {
		if (orientation == Qt::Horizontal) {
			switch (section) {
				case 0: return QString("name");
				case 1: return QString("local path");
			}
		} else if (orientation == Qt::Vertical) {
			return QString::number(section+1);
		}
	}
    return QVariant();
}


int FileTableModel::rowCount(const QModelIndex &/*parent*/) const
{
	return rows_.size();
}

int	FileTableModel::columnCount(const QModelIndex &/*parent*/) const
{
	return 2;
}

Qt::ItemFlags FileTableModel::flags(const QModelIndex &index) const
{
    if (!index.isValid()) return Qt::NoItemFlags;

	if( index.column() == 0  and  editable_ ) return Qt::ItemIsEditable | QAbstractItemModel::flags(index);

	return QAbstractItemModel::flags(index);
}

QVariant FileTableModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid()) return QVariant();
    if (role != Qt::DisplayRole && role != Qt::EditRole) return QVariant();

	//if( index.column() == 0 ) {
	//return QStringLiteral("r%1 c%2").arg(index.row()).arg(index.column());

	if( index.row() < static_cast<int>( rows_.size() ) ) {

		if( index.column() == 0 ) return rows_[index.row()].name;
		else if ( index.column() == 1 ) return rows_[index.row()].path;
	}
	return QVariant();
}

bool FileTableModel::setData(const QModelIndex &index, const QVariant &value, int role)
{
    if (role != Qt::EditRole) return false;

	QString s = value.toString();
	if( not s.isEmpty()  and  static_cast<std::size_t>(index.row()) < rows_.size() ) {
		QString previous_value = rows_[index.row()].name;
		//QTimer::singleShot(0, this, [=]() { Q_EMIT this->rename_file(previous_value, s); } );  //SLOT( rename(previous_value, s) ) );
		Q_EMIT this->rename_file(previous_value, s);
	}

    return false;
}


void FileTableModel::update_from_task(Task const &task)
{
	beginResetModel();

	editable_ = task.state() == Task::State::_none_;

	QDir project_path;
	if( auto project = task.project() ) project_path = QDir( QFileInfo( project->file_name() ).dir() );

	auto const & files = task.files();

	rows_.clear();
	rows_.reserve( files.size() );
	for(auto const & it : files) rows_.push_back( Row{it.first, project_path.relativeFilePath( it.second->file_name() ) } );

	endResetModel();
}


static std::map<QString, Task::State> const _String_to_Task_State_ = {
	{"none",     Task::State::_none_},
	{"draft",    Task::State::_draft_},
	{"queued",   Task::State::_queued_},
	{"running",  Task::State::_running_},
	{"finished", Task::State::_finished_},
	{"unknown",  Task::State::_unknown_},
};

QString Task::to_string(Task::State state)
{
	for(auto it = _String_to_Task_State_.begin(); it != _String_to_Task_State_.end(); ++it) {
		if( it->second == state ) return it->first;
	}
	//return State::_unknown_;
	throw std::out_of_range("");
}

Task::State Task::from_string(QString const &s)
{
	auto it = _String_to_Task_State_.find(s);
	if( it != _String_to_Task_State_.end() ) return it->second;
	return Task::State::_unknown_;
	//return _String_to_Task_State_.at(s);
}


Task::Task() : syncer_(this)
{
	app_ = "docking_protocol";
	version_ = "master";
	connect(&files_model_, SIGNAL( rename_file(QString const &, QString const &) ), this, SLOT( rename_file(QString const &, QString const &) ) );
}

//Task::Task(QUuid _node_id) : Node(_node_id)
// Task::Task(QString const &description) : description_(description)
// {
// 	//input(); // forcing creation in input node
// 	//output(); // forcing creation in output node
// }

Task::~Task()
{
	//qDebug() << "Task::~Task()";
}

void Task::add_file(QString const &name, FileSP const &file)
{
	files_[name] = file;
	files_model_.update_from_task(*this);
}

/// delete file from task, return true if file was in task files
bool Task::delete_file(QString const &name)
{
	auto it = files_.find(name);

	if( it == files_.end() ) return false;

	files_.erase(it);
	files_model_.update_from_task(*this);

	return true;
}


void Task::name(QString const &n)
{
	if( name_ != n ) { name_ = n; Q_EMIT changed(); }
}

void Task::app(QString const &a)
{
	if( app_ != a ) { app_ = a; Q_EMIT changed(); }
}


void Task::version(QString const &v)
{
	if( version_ != v ) { version_ = v; Q_EMIT changed(); }
}

void Task::flags(QString const &f)
{
	if( flags_ != f ) { flags_ = f; Q_EMIT changed(); }
}



void Task::description(QString const &d)
{
	if( description_ != d ) {
		description_ = d;
		Q_EMIT changed();
	}
}

void Task::nstruct(int n)
{
	if( nstruct_ != n ) {
		nstruct_ = n;
		Q_EMIT changed();
	}
}


void Task::subscribe()
{
	// if(root_) {
	// 	qDebug() << "Task[" << task_id() << "]: Name:" << name_ << " Task::subscribe";
	// 	Updater::get().subscribe( root_.get() );
	// }
	syncer_.subscribe();
}

/// return task UUID if task was already submitted other wise return NULL UUID
QString Task::task_id() const
{
	// if(root_) return root_->node_id();
	// else return QUuid();
	return task_id_.toString();
}

bool Task::is_syncing() const
{
	return syncing_progress_.first != syncing_progress_.second;
}

// return <current_value, max_value> for syncing operation progress
std::pair<int, int> Task::syncing_progress() const
{
	return syncing_progress_;
}

// set <current_value, max_value> for syncing operation progress
void Task::syncing_progress(int value, int max)
{
	auto syncing_progress = std::make_pair(value, max);
	if( syncing_progress_ != syncing_progress ) {
		syncing_progress_ = syncing_progress;
		Q_EMIT syncing();
	}
}

void Task::syncing_progress_advance(int value)
{
	syncing_progress(syncing_progress_.first + value, syncing_progress_.second);
}


void Task::submit(QString const & queue)
{
	// queue_ = queue;
	// state_ = State::_draft_;

	// create_sync_tree();
	// qDebug() << "Task[" << task_id() << "]: Task::submit() Name:" << name_ << " queue:" << queue_;

	// connect(root_.get(), SIGNAL(tree_synced()), this, SLOT(draft_to_queued()));

	// root_->data_is_fresh(true);

	// Q_EMIT changed();

	syncer_.submit(queue);
}


void Task::rename_file(QString const &previous_value, QString const &new_value)
{
	{ qDebug() << "Task::rename_file: " << previous_value << " -> " << new_value; }

	auto it = files_.find(previous_value);
	if( it != files_.end()  and  (not new_value.isEmpty() ) ) {
		auto file = it->second;
		files_.erase(it);
		files_.emplace(new_value, file);
		files_model_.update_from_task(*this);
	}
}

QJsonValue Task::task_data()
{
	//qDebug() << "Task[" << task_id() << "]: task_data()...";
	QJsonObject r;
	r["name"] = name_;
	r["state"] = to_string(state_);
	r["queue"] = queue_;
	r["description"] = description_;
	r["nstruct"] = nstruct_;
	r["flags"] = flags_;
	r["app"] = app_;
	r["version"] = version_;

	//QJsonDocument jd = QJsonDocument::fromVariant(r);
	//return jd.toJson(QJsonDocument::Compact);

	//qDebug() << "Task[" << task_id() << "]: task_data():" << r;

	return r;
}

void Task::task_data(QJsonValue const &jv)
{
	qDebug() << "Task[" << task_id() << "]: task_data(...): " << jv;

	bool changed = false;
	QJsonObject const o = jv.toObject();

	auto nm = o["name"];
	if( nm.isString() ) {
		auto name = nm.toString();
		if( name != name_) {
			name_ = name;
			changed = true;
		}
	}

	QJsonValue st = o["state"];
	if( st.isString() ) {
		State state = from_string( st.toString() );

		if( state != State::_unknown_  and  state != state_) {
			state_ = state;
			changed = true;
		}
	}

	auto ds = o["description"];
	if( ds.isString() ) {
		auto description = ds.toString();
		if( description != description_) {
			description_ = description;
			changed = true;
		}
	}

	auto q = o["queue"];
	if( q.isString() ) {
		auto queue = q.toString();
		if( queue != queue_) {
			queue_ = queue;
			changed = true;
		}
	}

	auto ns = o["nstruct"];
	if( ns.isDouble() ) {
		auto nstruct = ns.toInt();
		if( nstruct != nstruct_) {
			nstruct_ = nstruct;
			changed = true;
		}
	}

	auto fl = o["flags"];
	if( fl.isString() ) {
		auto flags = fl.toString();
		if( flags != flags_) {
			flags_ = flags;
			changed = true;
		}
	}

	auto ve = o["version"];
	if( ve.isString() ) {
		auto version = ve.toString();
		if( version != version_) {
			version_ = version;
			changed = true;
		}
	}

	auto ap = o["app"];
	if( ap.isString() ) {
		auto app = ap.toString();
		if( app != app) {
			app_ = app;
			changed = true;
		}
	}

	auto ti = o["task_id"];
	if( ti.isDouble() ) {
		auto task_id = ti.toInt();
		if( task_id != task_id_) {
			task_id_ = task_id;
			changed = true;
		}
	}


	if(changed) Q_EMIT this->changed();
}

// bool Task::operator ==(Task const &rhs) const
// {
// 	if( state_ != rhs.state_ ) return false;

// 	if( flags_ != rhs.flags_ ) return false;

// 	if( name_ != rhs.name_ ) return false;

// 	if( description_ != rhs.description_ ) return false;

// 	if( nstruct_ != rhs.nstruct_ ) return false;

// 	// if( bool(root_) xor bool(rhs.root_) ) return false;
// 	// if( bool(root_) and bool(rhs.root_)  and  *root_ != *rhs.root_ ) return false;

// 	if( syncer_ != rhs.syncer_ ) return false;

// 	return true;
// }

quint32 const _Task_stream_version_ = 0x00000001;

QDataStream &operator<<(QDataStream &out, Task const&t)
{
	out.setVersion(QDataStream::Qt_5_6);

	out << ui::util::_Task_magic_number_;
	out << _Task_stream_version_;

	using I = std::underlying_type<Task::State>::type;
	qint32 state = static_cast<I>(t.state_);
	out << state;

	out << t.task_id_;
	out << t.name_;
	out << t.app_;
	out << t.version_;
	out << t.nstruct_;
	out << t.flags_;
	out << t.description_;

	out << t.queue_;

	out << t.files_;

	out << t.syncer_;

	return out;
}


QDataStream &operator>>(QDataStream &in, Task &t)
{
	in.setVersion(QDataStream::Qt_5_6);

	quint64 magic;
	in >> magic;
	if( magic != ui::util::_Task_magic_number_ ) throw ui::util::BadFileFormatException( QString("Invalid _Task_magic_number_: read %1, was expecting %2...").arg(magic).arg(ui::util::_Task_magic_number_) );

	quint32 version;
	in >> version;
	if( version != _Task_stream_version_ ) throw ui::util::BadFileFormatException( QString("Invalid _Task_stream_version_: read %1, was expecting %2...").arg(magic).arg(_Task_stream_version_) );

	qint32 state;
	in >> state;
	t.state_ = static_cast<Task::State>(state);

	in >> t.task_id_;
	in >> t.name_;
	in >> t.app_;
	in >> t.version_;
	in >> t.nstruct_;
	in >> t.flags_;
	in >> t.description_;

	in >> t.queue_;

	in >> t.files_;
	t.files_model_.update_from_task(t);

	in >> t.syncer_;

	// moved to TaskSyncer operator>>  if( t.state_ != Task::State::_none_ ) t.subscribe();

	Q_EMIT t.changed();

	return in;
}



// QByteArray FileMixin::data() const
// {
// 	std::string fname = file_name();
// 	uint64_t size = fname.size();
// 	QByteArray r;
// 	r.append(reinterpret_cast<char*>(&size), sizeof(size));
// 	r.append(fname.c_str(), size);
// 	r.append( file_data() );
// 	return r;
// }
// void FileMixin::data(QByteArray const &)
// {
// }


} // namespace task
} // namespace ui
