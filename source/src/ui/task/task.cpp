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

#include <QJsonDocument>
#include <QJsonObject>
//#include <QVariant>
#include <QDataStream>
#include <QFile>

namespace ui {
namespace task {

Node::Key const _input_key_ {"input"};
Node::Key const _output_key_{"output"};

Node::Key const _script_key_ {"script"};

File::File(QString const &file_name)
{
	init_from_file(file_name);
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
		file_data_ = file.readAll();
	}
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

static std::map<QString, Task::State> const _String_to_Task_State_ = {
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


//Task::Task(QUuid _node_id) : Node(_node_id)
Task::Task(QString const &description) : description_(description)
{
	//input(); // forcing creation in input node
	//output(); // forcing creation in output node
}

Task::~Task()
{
	//qDebug() << "Task::~Task()";
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
	// if(auto sp = output_node_.lock() ) {
	// 	qDebug() << "Task[" << task_id() << "]: Name:" << (project_ ? *project_->find(this) : "''") << " Task::subscribe";
	// 	Updater::get().subscribe( sp.get() );
	// }
	if(root_) {
		qDebug() << "Task[" << task_id() << "]: Name:" << (project_ ? *project_->find(this) : "''") << " Task::subscribe";
		Updater::get().subscribe( root_.get() );
	}
}


// void Task::input(File &&input)
// {
// 	input_ = std::move(input);
// }

// File &Task::input()
// {
// 	return input_;
// }

// void Task::script(File &&script)
// {
// 	script_ = std::move(script);
// }

// File &Task::script()
// {
// 	return script_;
// }

// void Task::script(File &&script)
// {
// 	script_ = std::move(script);
// }

// File &Task::script()
// {
// 	return script_;
// }


// /// return pointer to 'input' node which will contain all input data
// Node &Task::input()
// {
// 	return at(_input_key_);
// }
// /// return pointer to 'output' node which will contain all input data
// Node &Task::output()
// {
// 	return at(_output_key_);
// }

/// return task UUID if task was already submitted other wise return NULL UUID
QUuid Task::task_id() const
{
	if(root_) return root_->node_id();
	else return QUuid();
}

bool Task::is_syncing() const
{
	if(root_) return root_->syncing(true);
	return false;
}

// return <current_value, max_value> for syncing operation progress
std::pair<int, int> Task::syncing_progress() const
{
	if(root_) {
		return std::make_pair(root_->syncing(true), root_->tree_size());
	}
	else return std::make_pair(0, 0);
}



void Task::submit(QString const & queue)
{
	queue_ = queue;
	state_ = State::_draft_;

	create_sync_tree();
	qDebug() << "Task[" << task_id() << "]: Task::submit() Name:" << (project_ ? *project_->find(this) : "''") << " queue:" << queue_;

	connect(root_.get(), SIGNAL(tree_synced()), this, SLOT(draft_to_queued()));

	root_->data_is_fresh(true);

	Q_EMIT changed();
}


void Task::draft_to_queued(void)
{
	if( root_ )  {
		qDebug() << "Task[" << task_id() << "]: draft_to_queued()...";
		disconnect(root_.get(), SIGNAL(tree_synced()), this, SLOT(draft_to_queued()) );

		state_ = State::_queued_;

		connect(root_.get(), SIGNAL(synced()), this, SLOT(post_submit()));
		root_->data_is_fresh(false);
	}
}

void Task::post_submit(void)
{
	qDebug() << "Task[" << task_id() << "]: post_submit()...";
	disconnect(root_.get(), SIGNAL(synced()), this, SLOT(post_submit()));
	Q_EMIT submitted();
	Q_EMIT changed();
	subscribe();
}


void Task::create_sync_tree()
{
	//QUuid task_id = QUuid::createUuid();  /// each distinct submitted Task _must_ have unique id

	root_ = std::make_shared<Node>("task", Node::Flags::data_in | Node::Flags::data_out | Node::Flags::topology_out/*, task_id*/);

	auto input_node = std::make_shared<Node>(Node::Flags::data_out | Node::Flags::topology_out);
	root_->add("input", input_node);
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

void Task::connect_task_and_nodes()
{
	if(root_) {
		TaskQP qtask(this);

		root_->get_data_callback( [qtask]() { if(qtask) return qtask->task_data(); else return QJsonValue(); } );
		root_->set_data_callback( [qtask](QJsonValue const &ba) { if(qtask) return qtask->task_data(ba); } );

		input_node_ = root_->leaf("input");
		if(auto input_node = input_node_.lock() ) input_node->get_data_callback( [qtask]() -> QJsonValue { if(qtask) return QJsonValue( QLatin1String( qtask->input_.data().toBase64() ) ); else return QJsonValue(); } );

		script_node_ = root_->leaf("script");
		if(auto script_node = script_node_.lock() )  script_node->get_data_callback( [qtask]() -> QJsonValue { if(qtask) return QJsonValue( QLatin1String( qtask->script_.data().toBase64() ) ); else return QJsonValue(); } );

		flags_node_ = root_->leaf("flags");
		if(auto flags_node = flags_node_.lock() ) flags_node->get_data_callback( [qtask]() -> QJsonValue { if(qtask) return QJsonValue( QLatin1String( qtask->flags_.data().toBase64() ) ); else return QJsonValue(); } );

		output_node_ = root_->leaf("output");
		if(auto output_node = output_node_.lock() ) {
			output_node->topology_updated_callback(
				[qtask](Node const *node, std::vector<QString> const & new_keys, std::vector<QString> const & errased_keys) -> void
				{ if(qtask) qtask->output_topology_updated(node, new_keys, errased_keys); }
		    );
			// connect(output_node.get(), SIGNAL(        topology_updated(Node const *, std::vector<QString> const & , std::vector<QString> const & ) ),
			// 		this,                SLOT( output_topology_updated(Node const *, std::vector<QString> const & , std::vector<QString> const & ) ));
		}

		connect(root_.get(), SIGNAL(tree_synced()), this, SIGNAL(changed()));

	} else {
		input_node_.reset();
		script_node_.reset();
		flags_node_.reset();
		output_node_.reset();
	}
}

/*void Task::assign_input_node()
{
	if(root_) input_node_ = root_->leaf("input");
	else input_node_.reset();
}
void Task::assign_script_node()
{
	if(root_) script_node_ = root_->leaf("script");
	else script_node_.reset();
}
void Task::assign_flags_node()
{
	if(root_) flags_node_ = root_->leaf("flags");
	else flags_node_.reset();
}
void Task::assign_output_node()
{
	if(root_) output_node_ = root_->leaf("output");
	else output_node_.reset();
}*/

QJsonValue Task::task_data()
{
	//qDebug() << "Task[" << task_id() << "]: task_data()...";
	QJsonObject r;
	r["state"] = to_string(state_);
	r["queue"] = queue_;
	r["description"] = description_;
	r["nstruct"] = nstruct_;

	//QJsonDocument jd = QJsonDocument::fromVariant(r);
	//return jd.toJson(QJsonDocument::Compact);
	return r;
}

void Task::task_data(QJsonValue const &jv)
{
	qDebug() << "Task[" << task_id() << "]: task_data(...): " << jv;

	bool changed = false;
	QJsonObject const o = jv.toObject();

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

	if(changed) Q_EMIT this->changed();
}

void Task::output_topology_updated(Node const *, std::vector<QString> const & new_keys, std::vector<QString> const & errased_keys)
{
	qDebug() << "Task::output_topology_updated: New keys: " << new_keys << " Errased keys:" << errased_keys;

	if(auto output_node = output_node_.lock() ) {
		TaskQP qtask(this);

		for(auto const &k : errased_keys) output_.erase(k);

		for(auto const &k : new_keys) {
			if(NodeSP leaf = output_node->leaf(k) ) {
				FileSP file_sp = std::make_shared<File>();
				file_sp->file_name(k);
				output_[k] = file_sp;

				qDebug() << "Task::output_topology_updated: setting set-data-callback for node: " << leaf->node_id();
				leaf->set_data_callback( [qtask, file_sp](QJsonValue const &jv) {
						//qDebug() << "Data update for output key(jv):" << jv;
						auto ba = QByteArray::fromBase64( jv.toVariant().toByteArray() );

						qDebug() << "Data update for output key:" << file_sp->file_name() << "  data: " << ba.left(64) << ( ba.size() > 64 ? "..." : "");
						if(qtask) {
							file_sp->data(ba);
							Q_EMIT qtask->changed();
						}
					} );

			}
		}
	}
}


bool Task::operator ==(Task const &rhs) const
{
	if( state_ != rhs.state_ ) return false;

	if( description_ != rhs.description_ ) return false;

	if( script_ != rhs.script_ ) return false;

	if( bool(root_) xor bool(rhs.root_) ) return false;
	if( bool(root_) and bool(rhs.root_)  and  *root_ != *rhs.root_ ) return false;

	return true;
}

quint64 const _Task_magic_number_   = 0xFFF2D9FACD279EE1;
quint32 const _Task_stream_version_ = 0x00000001;

QDataStream &operator<<(QDataStream &out, Task const&t)
{
	out.setVersion(QDataStream::Qt_5_6);

	out << _Task_magic_number_;
	out << _Task_stream_version_;

	using I = std::underlying_type<Task::State>::type;
	qint32 state = static_cast<I>(t.state_);
	out << state;

	out << t.description_;

	out << t.input_;
	out << t.script_;
	out << t.flags_;

	out << t.output_;

	out << (t.root_ != nullptr);
    if( t.root_ != nullptr ) out << *t.root_;

	return out;
}


QDataStream &operator>>(QDataStream &in, Task &t)
{
	in.setVersion(QDataStream::Qt_5_6);

	quint64 magic;
	in >> magic;
	if( magic != _Task_magic_number_ ) throw ui::util::BadFileFormatException();  //TaskFileFormatException();

	quint32 version;
	in >> version;
	if( version != _Task_stream_version_ ) throw ui::util::BadFileFormatException();  // TaskBadFileFormatException();

	qint32 state;
	in >> state;
	t.state_ = static_cast<Task::State>(state);

	in >> t.description_;

	in >> t.input_;
	in >> t.script_;
	in >> t.flags_;

	in >> t.output_;

	bool root;
	in >> root;
	if(root) {
		t.root_ = std::make_shared<Node>();
		in >> *t.root_;

		t.connect_task_and_nodes();
		// t.assign_input_node();
		// t.assign_script_node();
		// t.assign_flags_node();
		// t.assign_output_node();

		if( t.state_ != Task::State::_draft_ ) t.subscribe();
	}

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
