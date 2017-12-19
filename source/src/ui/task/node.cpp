// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/task/node.cpp
/// @brief  Node class for ui library.
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu).

#include <ui/task/node.h>

#include <ui/config/util.h>

#include <QJsonDocument>
#include <QJsonObject>
#include <QTimer>
#include <QDataStream>
#include <QAuthenticator>

#include <ctime>
#include <algorithm>

namespace ui {
namespace task {

QNetworkAccessManager * network_access_manager()
{
	static QSharedPointer<QNetworkAccessManager> manager = QSharedPointer<QNetworkAccessManager>(new QNetworkAccessManager);

	QAbstractSocket::connect(manager.data(), &QNetworkAccessManager::authenticationRequired, /*SIGNAL( authenticationRequired(QNetworkReply *, QAuthenticator *) )*/
							 [=](QNetworkReply * reply, QAuthenticator * authenticator) {
								 qDebug() << "authenticationRequired! " << reply->url().host();

								 //UserCredentials c = get_user_credentials( reply->url().host() );
								 config::UserCredentials c = config::get_user_credentials();
								 authenticator->setUser(c.user);
								 authenticator->setPassword(c.password);
							 });

	return manager.data();
}

// QString const data_url = "http://127.0.0.1:64078/data";
// QString const node_url = "http://127.0.0.1:64078/node";
// QString const meta_url = "http://127.0.0.1:64078/meta";
// QString const sync_url = "http://127.0.0.1:64078/sync";
// QString const updates_url = "http://127.0.0.1:64078/updates";

QString server_url()
{
	return "http://ui.graylab.jhu.edu/api";
}

Updater & Updater::get()
{
	static std::shared_ptr<Updater> o( new Updater() );

	return *o.get();
}


void Updater::subscribe(Node *node)
{
	QUuid id = node->node_id();

	auto range = subscribers_.equal_range(id);
    for(auto it = range.first; it != range.second; ++it)
		if( it->second == node ) return;

	subscribers_.emplace_hint(range.second, std::make_pair(id, node));
	//connect(node, SIGNAL( destroyed(QObject *) ), this, SLOT(unsubscribe(QObject *)));

	if(listen_) {
		abort_network_operation();
		listen_to_updates();
	}

	// if( not subscribers_.contains(id) ) {
	// 	subscribers_[id] = node;
	// 	connect(node, SIGNAL( destroyed(QObject *) ), this, SLOT(unsubscribe(QObject *)));
	// 	if(listen_) {
	// 		abort_network_operation();
	// 		listen_to_updates();
	// 	}
	// }
}

// void Updater::unsubscribe(Node *node)
// {
// 	qDebug() << "Updater::unsubscribe for node_id: " << node->node_id() << " node:" << node;
// 	auto range = subscribers_.equal_range( node->node_id() );
//     for(auto it = range.first; it != range.second; ++it) {
// 		if( it->second == node ) {
// 			subscribers_.erase(it);
// 			if(listen_) {
// 				abort_network_operation();
// 				listen_to_updates();
// 			}
// 			break;
// 		}
// 	}
// }

// void Updater::unsubscribe(QObject *obj)
// {
// 	// Fixme: refactor to create new list to prevent potential race conditions
// 	// we do not know which node this was (since Node part of the object is already destroyed) so we will just check all subscribers
//     for(auto it = subscribers_.rbegin(); it != subscribers_.rend(); ++it) {
// 		if( not it->second ) {
// 			qDebug() << "Updater::unsubscribe for obj: " << obj << " Node: " << it->first;
// 			subscribers_.erase( std::next(it++).base() );
// 		}
// 	}
// 	// keep current listener connection if more subscribers remains
// 	//if( listen_   and  subscribers_.empty() ) {
// 	if(listen_) {
// 		qDebug() << "Updater::unsubscribe: new subscribers list:" << subscribers_;
// 		abort_network_operation();
// 		listen_to_updates();
// 	}
// 	//if( Node *node = qobject_cast<Node *>(obj) ) unsubscribe(node); <-- does not work because Node instance is already destroyed at this point
// 	// for(auto &k : subscribers_.keys() ) {
// 	// 	if(not subscribers_[k]) {
// 	// 		qDebug() << "Updater::unsubscribe for obj: " << obj << " Node: " << k;
// 	// 		subscribers_.remove(k);
// 	// 		if(listen_) {
// 	// 			abort_network_operation();
// 	// 			listen_to_updates();
// 	// 		}
// 	// 		break;  // iterator will be invalidated after this call. However there is no point of further iterating since all other subscribers should be fine...
// 	// 	}
// 	// }
// }


bool Updater::erase_deleted_subscribers()
{
	bool flag = false;
	for(auto it = subscribers_.begin(); it != subscribers_.end(); ) {
		if( not it->second ) {
			qDebug() << "Updater::erase_deleted_subscribers erasing:" << it->first;
			it = subscribers_.erase(it);
			flag = true;
		}
		else ++it;
	}
	return flag;
}


void Updater::listen(bool l)
{
	if(l) {
		if(not listen_) {
			listen_to_updates();
			listen_ = true;
		}
	} else {
		if(listen_) {
			listen_ = false;
			abort_network_operation();
		}
	}
}

void Updater::abort_network_operation()
{
	if(reply_) {
		//qDebug() << "Updater::abort_network_operation ABORTING previous request!";
		disconnect(reply_, SIGNAL(finished()), 0, 0); // we have to disconnet 'finished' signal because call to abort() will issue 'finished'
		Q_EMIT reply_->abort();
		reply_->deleteLater();
		reply_ = nullptr;
	}
}

void Updater::listen_to_updates()
{
	abort_network_operation();
	if( subscribers_.empty() ) return;
	//

	//qDebug() << "Updater::listen_to_updates...";

	// QVariantMap r;
	// r["node_id"]   = node_id_.toByteArray();
	//QList<QUuid> keys;
	QVariantList l;
	for(auto it = subscribers_.begin(); it != subscribers_.end(); ++it) {
		l << it->first.toString().mid(1, 36);
		if( it != subscribers_.begin() ) {
			auto jt = it;
			++it;
			for(; it != subscribers_.end(); ++it) {
				if( jt->first != it->first ) break;
			}
			if( it == subscribers_.end() ) break;
		}
	}


    //for(auto const & e : keys) l << e.toString().mid(1, 36);;
	QJsonDocument jd = QJsonDocument::fromVariant(l);
	//qDebug() << "Updater::listen_to_updates jd: " << jd;

    QUrl url( server_url() + "/updates" );
	QNetworkRequest request(url);
	request.setHeader( QNetworkRequest::ContentTypeHeader, "application/json; charset=utf-8" );
	//request.setHeader( QNetworkRequest::ContentTypeHeader, "application/octet-stream" );
    reply_ = network_access_manager()->post( request, jd.toJson(QJsonDocument::Compact));

    connect(reply_, SIGNAL(finished()), this, SLOT(updates_finished()));
    if( reply_->isFinished() ) { updates_finished(); }

	connect(reply_, SIGNAL(readyRead()), this, SLOT(update_data_ready()));
}

void Updater::updates_finished()
{
	if(reply_) {
		reply_->deleteLater();
		reply_ = nullptr;
		qDebug() << "Updater::updates_finished! Initiating retry...";
		//if(listen_) QTimer::singleShot(64*1000, this, SLOT( listen_to_updates() ));
		if(listen_) QTimer::singleShot(5*1000, this, SLOT( listen_to_updates() ));
	}
}

void Updater::update_data_ready()
{
	//qDebug() << "Updater::update_data_ready...";

	while( reply_->canReadLine() ) {
		QByteArray d = reply_->readLine();
		//qDebug() << "Node::update_data_ready: got raw data: " << d;
		if( not d.isEmpty() ) d.resize( d.size() -1 ); // removing '\n' at the end
		if( d.size() >= 36+1+ 36+1 +1  and  d.at(36) == ' '  and  d.at(36+36+1) == ' ' ) {
			QUuid tree_id( d.left(36) );
			QUuid node_id( d.mid(37, 36) );
			TimeStamp time = d.mid(36+36+2).toInt();
			//qDebug() << "Node::update_data_ready: got data: " << tree_id << node_id << " time:" << time;

			auto range = subscribers_.equal_range(tree_id);
			for(auto it = range.first; it != range.second; ++it) {
				QPointer<Node> node(it->second);
				if( node ) {
					node->execute( node_id, [time](Node &n) {
							if( n.local_modification_time() < time  and
								(has_permission(n.flags(), Node::Flags::topology_in) or has_permission(n.flags(), Node::Flags::data_in) ) ) {
								qDebug() << "Node:" << n.node_id() << " update_if_outdatad: honoring update request: local_modification_time=" << n.local_modification_time() << " server_modification_time=" << time;
								n.data_is_outdated();
							} else {
								//qDebug() << "Node::update_if_outdatad: ignoring update request: local_modification_time=" << n.local_modification_time() << " server_modification_time=" << time;
							}
						});
				}
			}

			// auto it = subscribers_.find(tree_id);
			// if( it != subscribers_.end() ) {
			// 	QPointer<Node> node (*it);
			// 	if( node ) {
			// 		node->execute( node_id, [time](Node &n) {
			// 				if( n.local_modification_time() < time ) {
			// 					qDebug() << "Node:" << n.node_id() << " update_if_outdatad: honoring update request: local_modification_time=" << n.local_modification_time() << " server_modification_time=" << time;
			// 					n.data_is_outdated();
			// 				} else {
			// 					//qDebug() << "Node::update_if_outdatad: ignoring update request: local_modification_time=" << n.local_modification_time_ << " server_modification_time=" << time;
			// 				}
			// 			});
			// 	}
			// 	else { // QPointer<Node> is nullptr so Node object was destroyed
			// 		subscribers_.remove(tree_id);
			// 		listen(listen_);
			// 	}
			// }
		}
		else {
			if( erase_deleted_subscribers() ) {
				if(listen_) {
					qDebug() << "Updater::update_data_ready: new subscribers list:" << subscribers_;
					abort_network_operation();
					listen_to_updates();
				}
				break; // abort_network_operation() will set reply_ to nullptr, so we have to break
			}
		}
	}
	// else {
	// 		qDebug() << "Project::update_data_ready: canReadLine() is false!";
	// }
}



Node::Node(Flags _flags, /*QUuid task_id,*/ QUuid _node_id) : Node("node", _flags, _node_id)
{
}

Node::Node(QString const &_type, Flags _flags, QUuid _node_id) : type_(_type), flags_(_flags), node_id_(_node_id), parent_(nullptr)
{
}

Node::~Node()
{
	//qDebug() << "Node::~Node(): " << node_id_;
	abort_network_operation();
}

// std::string Node::type() const
// {
// 	return type_; // type for generic node that does act as simple container with no extra functionality
// }


// QByteArray Node::data() const
// {
// 	return QByteArray();
// }


// void Node::data(QByteArray const &)
// {
// }


void Node::add(Key const &key, NodeSP const &node)
{
	if( node->parent_ ) node->parent_->erase( node.get() );

	leafs_[key] = node;
    node->parent_ = this;
}


// GUI helper function
// return i'th leaf or nullptr if no such leaf exists
NodeSP Node::leaf(int i) const
{
	if( i < 0  or  i >= size() ) return NodeSP();

    return std::next( leafs_.begin(), i)->second;
}


// return leaf assosiate with key or nullptr if no such leaf exists
NodeSP Node::leaf(Key const &key) const
{
	auto it = leafs_.find(key);
    if (it == leafs_.end()) return NodeSP();
	else return it->second;
}


// GUI helper function
// return key of given leaf or nullptr if leaf could not be found
Node::Key const * Node::find(Node *leaf) const
{
    auto it = find_if(leafs_.begin(), leafs_.end(), [&](Map::value_type const &p) { return p.second.get() == leaf; } );

	if( it == leafs_.end() ) return nullptr;
	else return &it->first;
}


// erase given leaf from this node
bool Node::erase(Node *leaf)
{
    auto it = find_if(leafs_.begin(), leafs_.end(), [&](Map::value_type const &p) { return p.second.get() == leaf; } );

	if( it == leafs_.end() ) return false;

    leafs_.erase(it);
	leaf->parent_ = nullptr;
   	return true;
}

// return index of node, return -1 if node could not be found
int Node::node_index(Node *node)
{
	auto it = find_if(leafs_.begin(), leafs_.end(), [&](Map::value_type const &p) { return p.second.get() == node; } );

	if( it != leafs_.end() ) return distance(leafs_.begin(), it);
	else return -1;
}

// execute given function on self and leafs
void Node::execute( std::function< void(Node &) > const &f )
{
	f(*this);
	for(auto & l : leafs_) l.second->execute(f);
}

// execute given function on node whos id is _node_id
void Node::execute( QUuid const & node_id, std::function< void(Node &) > const &f )
{
	if( node_id_ == node_id) f(*this);
	else for(auto & l : leafs_) l.second->execute(node_id, f);
}


// return number of nodes in tree
int Node::tree_size() const
{
	int r = 0;
	for(auto const &it : leafs_) r += 1 + it.second->tree_size();
	return r;
}


void Node::abort_network_operation()
{
	if(reply_) {
		qDebug() << "Node::abort_network_operation ABORTING previous request!";
		disconnect(reply_, SIGNAL(finished()), 0, 0); // we have to disconnet 'finished' signal because call to abort() will issue 'finished'
		Q_EMIT reply_->abort();
		reply_->deleteLater();
		reply_ = nullptr;
	}
}

int Node::syncing(bool recursive) const
{
	int r = 0;
	if(reply_) {
		//qDebug() << "Node " << node_id_ << " is still syncing...";
		++r;
	}

	if(recursive) {
		for(auto const &l : leafs_) {
			r += l.second->syncing(true);
		}
	}
	return r;
}

void Node::node_synced()
{
	if( parent_ ) parent_->node_synced();
	else {
		//qDebug() << "node_synced: root:" << node_id_;
		if( not syncing(true) ) Q_EMIT tree_synced();
	}
}


bool Node::operator ==(Node const &rhs) const
{
	//qDebug() << "Node::operator ==: flags_: " << int(flags_) << " " << int(rhs.flags_);
	if( flags_ != rhs.flags_ ) return false;

	//qDebug() << "Node::operator ==: type_: " << QString::fromStdString(type_) << " " << QString::fromStdString(rhs.type_);
	if( type_ != rhs.type_ ) return false;

	if( node_id_ != rhs.node_id_ ) return false;

	if( local_modification_time_ != rhs.local_modification_time_ ) return false;

	//qDebug() << "Node::operator ==: node_id_ equal";

	//parent_

	if( leafs_.size() != rhs.leafs_.size() ) return false;

	if( not std::equal( leafs_.begin(), leafs_.end(), rhs.leafs_.begin(),
					[](Map::value_type const &l, Map::value_type const &r) {
						if( l.first  != r.first  ) return false;
						//qDebug() << "Node::operator ==: second: " << l.second.get() << " == " << r.second.get() << " --> " << (*l.second == *r.second);

						if( *l.second != *r.second ) return false;
						return true;
					})
		) {
		return false;
	}

	// if( leafs_.size() != rhs.leafs_.size() ) return false;
	// auto l =     leafs_.begin();
	// auto r = rhs.leafs_.begin();
	// while( l != leafs_.end() ) {
	// 	if( l->first != r->first ) return false;
	// 	if( l->second != r->second ) return false;
	// 	++l; ++r;
	// }

	return true;
}


void Node::data_is_fresh(bool recursive)
{
	abort_network_operation();

    local_modification_time_ = get_utc_timestamp();

	recursive_ = recursive;

	//qDebug() << "Node::data_is_fresh node=" << node_id << " timestamp:" <<  data_local_modification_time << " initiating uploading...";

    QVariantMap r;
	r["node_id"]   = node_id_.toByteArray();
	r["parent_id"] = parent_ ? parent_->node_id_ : QVariant();
    r["type"] = type_;

	if( has_permission(flags_, Flags::topology_out ) ) {
		QVariantMap jleafs;
		for(auto const &l : leafs_) {
			QVariantMap leaf;
			leaf["type"] = l.second->type_;
			leaf["leaf_id"] = l.second->node_id_.toString().mid(1, 36);
			jleafs[l.first] = leaf;
		}
		r["leafs"]     = jleafs;
	}

	if( has_permission(flags_, Flags::data_out ) ) {
		QJsonValue jv;
		if(get_data_) {
			jv = get_data_();
			//if( static_cast<QMetaType::Type>( qvd.type() ) == QMetaType::QByteArray ) qvd = qvd.toByteArray().toBase64();
			//qDebug() << "Node:" << node_id_ << " Setting custom data: " << qvd;
		}
		r["data"] = jv;
	}

	// QByteArray d;
	// if(get_data_) d = get_data_();
	// r["data"] = d.toBase64();

	QJsonDocument jd = QJsonDocument::fromVariant(r);
	//qDebug() << "jd: " << jd;
	//qDebug() << "jd: " << jd.toBinaryData();

	QUrl url( server_url() + "/node" + /*'/' + task_id_.toString().mid(1, 36) + */ '/' + node_id_.toString().mid(1, 36) + '/' + QString::number(local_modification_time_) );

	// //QUrlQuery q;  q.addQueryIte ("modification_time", QString::number(data_local_modification_time) );
	// //url.setQuery(q);
	// //qDebug() << "Using query: " << url;

	QNetworkRequest request(url);
	//request.setHeader( QNetworkRequest::ContentTypeHeader, "application/octet-stream" );
	request.setHeader( QNetworkRequest::ContentTypeHeader, "application/json; charset=utf-8" );



	// todo: add Content-Encoding gzip: request.setRawHeader(QByteArray("Content-Encoding"), QByteArray("gzip"));
	// QByteArray d = qCompress(jd.toJson(QJsonDocument::Compact), 5);

    reply_ = network_access_manager()->post( request, jd.toJson(QJsonDocument::Compact));

    connect(reply_, SIGNAL(finished()), this, SLOT(data_upload_finished()));

    if( reply_->isFinished() ) { data_upload_finished(); }

	/* moved to data_upload_finished
	if(recursive) {
		int delay = 1;
		for(auto &l : leafs_) {
			QTimer::singleShot(delay*100.0, l.second.get(), SLOT(data_is_fresh()));
			//Node * p = l.second.get();
			//QTimer::singleShot(delay*100.0, p, [p]() {p->data_is_fresh(true); } );
			++delay;
		}
	}
	*/
}

void Node::data_upload_finished(void)
{
	if(reply_) {
 		reply_->deleteLater();
		//qDebug() << "Node::data_upload_finished for " << node_id_;

		QVariant qv_status_code = reply_->attribute( QNetworkRequest::HttpStatusCodeAttribute );
		if( qv_status_code.isValid() ) {
			// if( int status = qv_status_code.toInt()) {
			// 	//qDebug() << "Node::upload_finished for " << node_id_ << " status:" << status;
			// }
		}

        if( reply_->error() == QNetworkReply::NoError ) {
			server_modification_time_ = local_modification_time_;

			if(recursive_) {
				for(auto &l : leafs_) {
					if( not l.second->syncing() ) {
						Q_EMIT l.second->data_is_fresh(true);
					}
				}
			}

			reply_ = nullptr;
			Q_EMIT synced();
			node_synced();

		} else {
            qDebug() << "Node::upload_finished with error: " << reply_->error() << " reply:" << reply_->readAll() << " initiating retry...";
			reply_ = nullptr;

			//QTimer::singleShot(get_retry_interval()*1000.0, this, SLOT(data_is_fresh()));
			QTimer::singleShot(get_retry_interval()*1000.0, [this]() { data_is_fresh(recursive_); });
		}
 	}
}


/// assume that server have an updated version of data with (blob_id, blob_modification_time) and initiate sync
void Node::data_is_outdated()
{
	qDebug() << "Node::data_is_outdated for" << node_id_;
	abort_network_operation();

        //server_modification_time = modification_time;

    QUrl url( server_url() + "/node/"  +/* project_id().toString().mid(1, 36) + '/' + */ node_id_.toString().mid(1, 36) );

	QNetworkRequest request(url);
	//request.setHeader( QNetworkRequest::ContentTypeHeader, "application/octet-stream" );

    reply_ = network_access_manager()->get(request);
    reply_->setReadBufferSize(0);  // setting unlimited buffere size

    connect(reply_, SIGNAL(finished()), this, SLOT(data_download_finished()));

    if( reply_->isFinished() ) { data_download_finished(); }
}


void Node::data_download_finished(void)
{
	if(reply_) {
 		reply_->deleteLater();

		qDebug() << "Node::data_download_finished for" << node_id_;

        if( reply_->error() == QNetworkReply::NoError ) {
			//reply_->deleteLater();
			QByteArray data = reply_->readAll();
			reply_ = nullptr;
			//qDebug() << "Node::data_download_finished data: " << data;

			QJsonDocument jd = QJsonDocument::fromJson(data);
			//qDebug() << "Node::data_download_finished json: " << jd;

			QJsonObject root = jd.object();

			update_from_json(root, true);

			node_synced();
			Q_EMIT synced();
		} else {
			if( reply_->error() == QNetworkReply::ContentGoneError ) { // node is no longer there...
				qDebug() << "Node::data_download_finished with error: " << reply_->error() << " aborting...";
				reply_ = nullptr;

				node_synced();
				Q_EMIT synced();
			} else {
				qDebug() << "Node::data_download_finished with error: " << reply_->error() << " initiating retry...";
				reply_ = nullptr;

				QTimer::singleShot(0, this, SLOT( data_is_outdated() ));
			}
		}
 	}
}


void Node::update_from_json(QJsonObject const &root, bool forced)
{
	if( root["modification_time"].isDouble() ) {
		server_modification_time_ = root["modification_time"].toDouble();
		//qDebug() << "Node::download_finished modification_time:" << server_modification_time_;

		if( local_modification_time_ < server_modification_time_  or  forced ) {

			if( has_permission(flags_, Flags::topology_in)  and  root["leafs"].isObject() ) {
				//qDebug() << "Node::download_finished updating topology for node " << node_id_;

				auto jo_leafs = root["leafs"].toObject();

				struct LeafInfo {
                    QUuid node_id;
					QString type;
				};

				std::map< QString, LeafInfo> new_leafs;

				for(auto  it = jo_leafs.constBegin(); it != jo_leafs.constEnd(); ++it) {
					//std::string key = it.key().toUtf8().constData();
					QString key = it.key();
                    LeafInfo li {
							QUuid( it.value().toObject()["leaf_id"].toString() ),
							it.value().toObject()["type"].toString().toUtf8().constData()
					};

					//qDebug() << it.key() << " : " << it.value();
					//qDebug() << key.c_str() << " : { " << li.node_id << ", " << li.type.c_str() << "}";
					new_leafs.insert( std::make_pair(key, li) );
				}

				std::vector<QString> new_keys, errased_keys;
				//std::transform(new_leafs.begin(), new_leafs.end(), back_inserter(new_keys), [](std::pair<const QString, LeafInfo> const &p) { return p.first; } );

                std::vector<NodeSP> new_nodes;

                //begin_topology_update(this);
				//begin_topology_update();

				auto it = leafs_.begin();

				for(auto nt = new_leafs.begin(); nt != new_leafs.end();) {
					//qDebug() << "Considering keys: " << it->first.c_str() << " " << nt->first.c_str();
					if( it == leafs_.end()  or  Map::key_compare()(nt->first, it->first) ) {
                        qDebug() << "Inserting new node: " << nt->first;

						//new_nodes.push_back( create_node(nt->second.type, nt->second.node_id, this) );
						NodeSP l = std::make_shared<Node>(Node::Flags::all, nt->second.node_id);
						l->type_ = nt->second.type;
						l->parent_ = this;
						new_nodes.push_back(l);

                        it = leafs_.insert( std::make_pair(nt->first, new_nodes.back() ) ).first;

						new_keys.push_back(nt->first);

						++it;
						++nt;
					}
					else if( not Map::key_compare()(it->first, nt->first) ) {
                        if( it->second->node_id_ != nt->second.node_id  /*or it->type != nt->type */) { // replacing node type without chaning node_id should not be legal so we avoid this testing
							qDebug() << "Replacing node: " << it->first;
							it = leafs_.erase(it);

							//new_nodes.push_back( create_node(nt->second.type, nt->second.node_id, this) );
							NodeSP l = std::make_shared<Node>(Node::Flags::all, nt->second.node_id);
							l->type_ = nt->second.type;
							l->parent_ = this;
							new_nodes.push_back(l);

							it = leafs_.insert( std::make_pair(nt->first, new_nodes.back()  ) ).first;

							new_keys.push_back(nt->first);
						}
						++it;
						++nt;
					}
					else {
                        qDebug() << "Errasing node: " << it->first;
						errased_keys.push_back(it->first);
						it = leafs_.erase(it);
					}
				}

				while( it != leafs_.end() ) {
					std::transform(it, leafs_.end(), back_inserter(errased_keys), [](typename Map::value_type const &p) { return p.first; } );
					it = leafs_.erase(it);
				}

                //end_topology_update();

				//Q_EMIT topology_updated(this, new_keys, errased_keys);
				if(topology_updated_) topology_updated_(this, new_keys, errased_keys);

                for(auto n : new_nodes) n->data_is_outdated();
			}

			if( has_permission(flags_, Flags::data_in) ) {
				//begin_topology_update(this);


				//data( QByteArray::fromBase64( root["data"].toVariant().toByteArray() ) );
				if(set_data_) {
					// QJsonValue d = root["data"];
					// if( not d.isUndefined() ) set_data_(d);

					QJsonValue d = root["data"];
					if( not d.isUndefined() ) {
						set_data_( d/*.toVariant()*/ );
					} else {
						//qDebug() << "Data-setting-callback COULD NOT BE CALLED because isUndefined() is true for 'data' key!!! node_id:" << node_id_;
					}
				} else {
					//qDebug() << "Data-setting-callback is not set for node " << node_id_;
				}


				//end_topology_update();
				//qDebug() << "Setting data: " << QByteArray::fromBase64( root["data"].toVariant().toByteArray() );
			}

			local_modification_time_ = server_modification_time_;
		}
	}
}






quint64 const _Node_magic_number_   = 0x1A5CD6BF94D5E57F;
quint32 const _Node_stream_version_ = 0x00000001;


QDataStream &operator<<(QDataStream &out, Node const&n)
{
	out.setVersion(QDataStream::Qt_5_6);

	out << _Node_magic_number_;
	out << _Node_stream_version_;

	out << n.type_;

	using I = std::underlying_type<Node::Flags>::type;
	qint32 flags = static_cast<I>(n.flags_);
	out << flags;

	out << n.node_id_;

	out << n.leafs_;

	out << static_cast<qint64>( n.local_modification_time_ );
	out << static_cast<qint64>( n.server_modification_time_ );

	return out;
}

QDataStream &operator>>(QDataStream &in, Node &n)
{
	//in.setVersion(QDataStream::Qt_5_6);

	quint64 magic;
	in >> magic;
	if( magic != _Node_magic_number_ ) throw ui::util::BadFileFormatException( QString("Invalid _Node_magic_number_: read %1, was expecting %2...").arg(magic).arg(_Node_magic_number_) );  // NodeFileFormatException();

	quint32 version;
	in >> version;
	if( version != _Node_stream_version_ ) throw ui::util::BadFileFormatException( QString("Invalid _Node_stream_version_: read %1, was expecting %2...").arg(version).arg(_Node_stream_version_) );  // NodeBadFileFormatException();

	in >> n.type_;

	qint32 flags;
	in >> flags;
	n.flags_ = static_cast<Node::Flags>(flags);

	in >> n.node_id_;

	in >> n.leafs_;

	for( auto & l : n.leafs_ ) l.second->parent_ = &n;

	qint64 local_modification_time;
	in >> local_modification_time;
	n.local_modification_time_ = local_modification_time;

	qint64 server_modification_time;
	in >> server_modification_time;
	n.server_modification_time_ = server_modification_time;

	// we do not touch other network related variables assuming that `n` was freshly created

	return in;
}


} // namespace task
} // namespace ui
