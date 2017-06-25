// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   ui/task/project_model.cpp
/// @brief  Project model class for ui library.
/// @author Sergey Lyskov (sergey.lyskov@jhu.edu).

#include <ui/task/project_model.h>

#include <ui/task/task.h>

#include <QDebug>


namespace ui {
namespace task {


ProjectModel::ProjectModel(QObject *parent)
    : QAbstractItemModel(parent)
{
	root_ = new Project( QUuid("AABCDEF0-1234-0000-0000-123456846422") );  // temp debug project with a known ID
	root_->setParent(this);

	root_->add( "Task-1", std::make_shared<Task>( QUuid::createUuid(), nullptr) );
	root_->add( "Task-2", std::make_shared<Task>( QUuid::createUuid(), nullptr) );

	qDebug() << "ProjectModel(...)";
}


QVariant ProjectModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (role == Qt::DisplayRole) {
		if (orientation == Qt::Horizontal) {
			switch (section) {
				case 0: return QString("key");
				case 1: return QString("type");
				case 2: return QString("value");
				case 3: return QString("uuid");
			}
		}
	}
    return QVariant();
}


Node *ProjectModel::get_item(const QModelIndex &index) const
{
    if (index.isValid()) {
        Node *item = static_cast<Node*>(index.internalPointer());
        if (item) return item;
    }
    return root();
}


QModelIndex ProjectModel::index(int row, int column, const QModelIndex &parent) const
{
    if (parent.isValid() && parent.column() != 0) return QModelIndex();

	Node *parent_node = get_item(parent);

    Node *leaf = parent_node->leaf(row);

    if (leaf) return createIndex(row, column, leaf);
    else return QModelIndex();
}


QModelIndex ProjectModel::parent(const QModelIndex &index) const
{
    if (!index.isValid()) return QModelIndex();

    Node *leaf = get_item(index);
    Node *parent = leaf->parent();

    if( parent == nullptr  or  leaf == nullptr ) return QModelIndex();

    return createIndex(parent->node_index(leaf), 0, parent);
}


int ProjectModel::rowCount(const QModelIndex &parent_index) const
{
    //if (!parent.isValid()) return 0;
	if( auto parent = get_item(parent_index) ) {
		qDebug() << "ProjectModel::rowCount: " << parent->size();
		return parent->size();
	}
	return 0;
}


int ProjectModel::columnCount(const QModelIndex &/*parent*/) const
{
    //if (!parent.isValid()) return 0;

    return 2;
}


QVariant ProjectModel::data(const QModelIndex &index, int role) const
{
    if (!index.isValid()) return QVariant();

    if (role != Qt::DisplayRole && role != Qt::EditRole) return QVariant();
	//if (role != Qt::DisplayRole ) return QVariant();

    Node *node = get_item(index);

    if( node->parent() == nullptr ) return QString("root"); // this branch is never executed in practice

	if( index.column() == 0 ) {
		Node::Key const * key = node->parent()->find(node);
		if( key ) {
                    //return QString("<%1, '%2'>") .arg(key->first) .arg( QString::fromStdString(key->second) );
                    //return QString::fromStdString(*key);
                    return *key;
		}
	}
    else if( index.column() == 1 ) {
		//qDebug() << "data: " << QString(node->data().toHex() );
		//return QString("%1 {%2}") .arg( QString::fromStdString( node->type() ) ).arg( QString(node->data().toHex() ) );
		return QString::fromStdString( node->type() );
	}
        /*
	else if( index.column() == 2 ) {
		return QString::fromStdString( node->to_string() );
	}

	else if( index.column() == 3 ) {
		return node->get_node_id().toString().mid(1, 36);
        }*/

	return QString("unknown");
}


Qt::ItemFlags ProjectModel::flags(const QModelIndex &index) const
{
    if (!index.isValid()) return Qt::NoItemFlags;

	//if( index.column() == 2 ) return Qt::ItemIsEditable | QAbstractItemModel::flags(index);

	return QAbstractItemModel::flags(index);
}

} // namespace task
} // namespace ui
