// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/task/node.h
/// @brief  Node class for ui library.
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu).

#ifndef UI_TASK_NODE_H
#define UI_TASK_NODE_H

#include <ui/util/exception.h>
#include <ui/task/util.h>

#include <QObject>
#include <QUuid>
#include <QString>
#include <QUuid>
#include <QPointer>
#include <QNetworkReply>
#include <QNetworkAccessManager>

#include <memory>
#include <type_traits>
#include <functional>
#include <map>

namespace ui {
namespace task {


QNetworkAccessManager * network_access_manager();

QString server_url();

class Node;
using NodeWP  = std::weak_ptr< Node >;
using NodeSP  = std::shared_ptr< Node >;
using NodeCSP = std::shared_ptr< Node const >;


class Updater final : public QObject
{
    Q_OBJECT

private:
	Updater() {}

public:

	static Updater & get();

	/// check/set if Updater is currently listen (connect to remote server) to updates
	bool listen() const { return listen_; }
	void listen(bool);

	/// subscribe and unsubscribe Node for updates
	void subscribe(Node *);
	//	void unsubscribe(Node *);


public Q_SLOTS:
	//void unsubscribe(QObject *);

private Q_SLOTS:
	void updates_finished();
	void update_data_ready();

	/// subscribe to server update channel. This signal can be used for all, even non-root nodes.
	void listen_to_updates();

private:
	// stop outgoing network operations if any
	void abort_network_operation();

	// erase deleted subscribers from watch list and return true if at least one subscribers was deleted
	bool erase_deleted_subscribers();

private:
	//QMap <QUuid, QPointer<Node> > subscribers_;
	std::multimap <QUuid, QPointer<Node> > subscribers_;

	bool listen_ = true;

	QPointer<QNetworkReply> reply_;
};


class Node final : public QObject
{
    Q_OBJECT

public:
	enum class Flags {
		none         = 0,
		data_in      = 1,
		data_out     = 2,
		topology_in  = 4,
		topology_out = 8,
		all          = data_in | data_out | topology_in | topology_out,
	};

	TimeStamp const _T_unknown_  = 0;
	TimeStamp const _T_outdated_ = 1;

	using Map = std::map< QString, NodeSP >;
	using Key = Map::key_type;

	using SetDataCallback = std::function< void (QJsonValue const &) >;
	using GetDataCallback = std::function< QJsonValue (void) >;

	using TopologyUpdatedCallback = std::function< void (Node const *, std::vector<QString> const & new_keys, std::vector<QString> const & errased_keys) >;

	//explicit Node(QUuid _node_id = QUuid());

	explicit Node(Flags _flags=Flags::none, QUuid _node_id = QUuid::createUuid());
	explicit Node(QString const &_type, Flags _flags, QUuid _node_id = QUuid::createUuid());

	~Node();


	/// return node type as string
	QString type() const { return type_; }

	QUuid node_id() const { return node_id_; }

	Flags flags() const { return flags_; }


	/// Add new leaf this node, If node was already added to some other Node then remove it first.
    void add(Key const &, NodeSP const &);

	// return number of elemnts stored in this node
	int size() const { return leafs_.size(); }

	// GUI helper function
	// return i'th leaf or nullptr if no such leaf exists
	NodeSP leaf(int i) const;

	// return leaf assosiate with key or nullptr if no such leaf exists
	NodeSP leaf(Key const &) const;

	// GUI helper function
	// return key of given leaf or nullptr if leaf could not be found
	Key const * find(Node *leaf) const;

	// erase given leaf from this node, do nothing if given node was not found. return true/false depending if element could be found
    bool erase(Node *leaf);

    Node * parent() const { return parent_; }

	// return index of node, return -1 if node could not be found
	int node_index(Node *node);


	TimeStamp local_modification_time() const { return local_modification_time_; }

	// execute given function on self and leafs
	void execute( std::function< void(Node &) > const &f );
	void execute( QUuid const & node_id, std::function< void(Node &) > const &f );


	// return a reference to Node assosiatate with Key. If not such node found - create it. If node found but type could not be cast to NodeType std::bad_cast is thrown
	// template <class NodeType=Node>
	// NodeType & at(Key const &key) {
	// 		Node *r = leaf(key);
	// 		if( r == nullptr ) {
	// 			NodeSP node = std::make_shared<NodeType>();
	// 			add(key, node);
	// 			r = node.get();
	// 		}
	// 		return dynamic_cast<NodeType&>(*r);
	// }

	void set_data_callback(SetDataCallback const &set_data) { set_data_ = set_data; }
	void get_data_callback(GetDataCallback const &get_data) { get_data_ = get_data; }

	void topology_updated_callback(TopologyUpdatedCallback const &topology_updated_callback) { topology_updated_ = topology_updated_callback; }


	bool operator ==(Node const &r) const;
	bool operator !=(Node const &r) const { return not (*this == r); }

	// serialization
	friend QDataStream &operator<<(QDataStream &, Node const&);
	friend QDataStream &operator>>(QDataStream &, Node &);


	// check if current node is syncing and if so return number of syncing nodes
	int syncing(bool recursive=false) const;

	// return number of nodes in tree
	int tree_size() const;

private:

	void abort_network_operation();

	void node_synced();


Q_SIGNALS:
	//void data(QByteArray const &);

	/// Emitted when node finished network operation
	void synced();

	/// Emitted when node and all its leafs finished synced operation
    void tree_synced();

	/// Emitted when node recieve new topology from remote server
	// moved to callback void topology_updated(Node const *, std::vector<QString> const & new_keys, std::vector<QString> const & errased_keys);

public Q_SLOTS:
	//QByteArray data() const;

	/// assume that current node data/topology is fresh and start upload process
	void data_is_fresh(bool recursive);


	/// assume that server have an updated version of data with (blob_id, blob_modification_time) and initiate sync
	void data_is_outdated();

	/// subscribe to server update channel. This signal can be used for all, even non-root nodes.
	//void listen_to_updates();


private Q_SLOTS:
	void data_upload_finished();
    void data_download_finished();
	void update_from_json(QJsonObject const &root, bool forced);

private:
	QString type_;

	Flags flags_ = Flags::none;

	/// UUID of Task to which this node belongs to
	//QUuid const task_id_;

	/// UUID of this node, use null UUID for a root node
	QUuid node_id_;

    Node * parent_ = nullptr;

    Map leafs_;

	SetDataCallback set_data_;
	GetDataCallback get_data_;

	TopologyUpdatedCallback topology_updated_;

	// sync state
	bool recursive_;
	QPointer<QNetworkReply> reply_;
	TimeStamp local_modification_time_  = _T_outdated_;
	TimeStamp server_modification_time_ = _T_unknown_;
};


inline Node::Flags operator&(Node::Flags lhs, Node::Flags rhs)
{
	using I = std::underlying_type<Node::Flags>::type;
	return static_cast<Node::Flags>(static_cast<I>(lhs) & static_cast<I>(rhs));
}
inline Node::Flags &operator&=(Node::Flags &lhs, Node::Flags rhs)
{
	using I = std::underlying_type<Node::Flags>::type;
	lhs = static_cast<Node::Flags>(static_cast<I>(lhs) & static_cast<I>(rhs));
	return lhs;
}

inline Node::Flags operator|(Node::Flags lhs, Node::Flags rhs)
{
	using I = std::underlying_type<Node::Flags>::type;
	return static_cast<Node::Flags>(static_cast<I>(lhs) | static_cast<I>(rhs));
}
inline Node::Flags &operator|=(Node::Flags &lhs, Node::Flags rhs)
{
	using I = std::underlying_type<Node::Flags>::type;
	lhs = static_cast<Node::Flags>(static_cast<I>(lhs) | static_cast<I>(rhs));
	return lhs;
}

inline bool has_permission(Node::Flags flag, Node::Flags permission)
{
	return (flag & permission) == permission;
}


// class NodeBadFileFormatException : public ui::util::BadFileFormatException
// {
// public:
//     void raise() const { throw *this; }
//     NodeBadFileFormatException *clone() const { return new NodeBadFileFormatException(*this); }
// };


} // namespace task
} // namespace ui

#endif // UI_TASK_NODE_H
