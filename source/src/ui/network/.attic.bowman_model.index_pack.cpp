#include "bowman_model.h"

#include <protocols/network/util.hh>

#include <QDebug>

using namespace protocols::network;


namespace ui {
namespace network {

BowmanModel::BowmanModel(QObject *parent) : QAbstractItemModel(parent)
{
	connect(&Bowman::bowman(), &Bowman::back_ends_changed, this, &BowmanModel::on_bowman_back_ends_changed);
}

QVariant BowmanModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (role == Qt::DisplayRole) {
		if (orientation == Qt::Horizontal) {
			switch (section) {
				case 0: return QString("back-end");
				case 1: return QString("mover");
			}
		} else if (orientation == Qt::Vertical) {
			return QString::number(section+1);
		}
	}
    return QVariant();
}

//quintptr const client_id_offset = 0x00010000  // 640k should be enough for anyone
struct IndexPack {

    // will usually occupy 4 bytes, so build on 32bit platform where quintptr is quint32 should be ok
    unsigned int level : 1, back_end : 8, function: 16+7;


	IndexPack(unsigned int level=0, unsigned int back_end=0, unsigned int function=0) : level(level), back_end(back_end), function(function) {}

	IndexPack(quintptr v) {
		static_assert( sizeof(IndexPack) <= sizeof(quintptr), "Platform is not compatible with IndexPack implementation!");

		*this = * reinterpret_cast<IndexPack*>(&v);
	}

	operator quintptr() {
		quintptr r = 0;

		* reinterpret_cast<IndexPack*>(&r) = *this;
		return r;
	}
};



QModelIndex BowmanModel::index(int row, int column, QModelIndex const & parent) const
{
	if ( !hasIndex(row, column, parent) ) return QModelIndex();


	//if (parent.isValid() && parent.column() != 0) return QModelIndex();

	if ( parent.isValid() ) { // child
		if( parent.column() == 0 ) {
			IndexPack i = parent.internalId();
			if( i.level == 0  and  row < 0x0FFFFF ) {
				i.level = 1;
				i.function = row;
				return createIndex(row, column, i);
			}
		}
	}
    else {  // root
		if( row < 1 ) return createIndex(row, column, IndexPack(0, row, 0) );
	}

	return QModelIndex();

	// return createIndex(row, column);
    // // if ( !parent.isValid() ) return createIndex(row, column); // root
    // // else
    // //     parentItem = static_cast<TreeItem*>(parent.internalPointer());
    // // TreeItem *childItem = parentItem->child(row);
    // // if (childItem)
    // //     return createIndex(row, column, childItem);
    // // else
    // //     return QModelIndex();

	// if (!index.isValid())
    //     return QModelIndex();

    // TreeItem *childItem = static_cast<TreeItem*>(index.internalPointer());
    // TreeItem *parentItem = childItem->parentItem();

    // if (parentItem == rootItem)
    //     return QModelIndex();

    // return createIndex(parentItem->row(), 0, parentItem);
}


QModelIndex BowmanModel::parent(const QModelIndex &index) const
{
	if (!index.isValid()) return QModelIndex();

	IndexPack i = index.internalId();

	if( i.level == 0 ) return QModelIndex();
	else {
		i.level = 0;
		i.function = 0;
		return createIndex(i.back_end, 0, i);
	}
}

int BowmanModel::rowCount(const QModelIndex &parent) const
{
    if ( parent.isValid() ) {
		IndexPack i = parent.internalId();

		if( i.level == 0 ) return functions_.size();
		return 0;
	}

	return 4;
	// IndexPack i = parent.internalId();

	// if( i.level ) return 10;
	// else return 5;


    // return functions_.size();
}

int BowmanModel::columnCount(const QModelIndex &parent) const
{
    //if( parent.isValid() ) return 0;

	return 1;
}

QVariant BowmanModel::data(const QModelIndex &index, int role) const
{
    if( !index.isValid() ) return QVariant();

	if (role != Qt::DisplayRole ) return QVariant();

	IndexPack i = index.internalId();

	return QString("level=%1, back_end=%2, function=%3").arg(i.level).arg(i.back_end).arg(i.function);

	int row = index.row();
	if( row < static_cast<int>( functions_.size() )  and  row >= 0 ) {
		if      ( index.column() == 0 ) return QString::fromStdString(functions_[row].name);
		else if ( index.column() == 1 ) return QString::fromStdString(functions_[row].hal_id);
	}

    return QVariant();
}

/// return FunctionIdentifier for given index, return nullptr if index is invalid
BowmanModel::FunctionIdentifier * BowmanModel::get_identifier(QModelIndex const &index)
{
	if( index.isValid()  and  index.row() < static_cast<int>( functions_.size() ) ) return &functions_[index.row()];
	else return nullptr;
}

void BowmanModel::on_bowman_back_ends_changed(Bowman const *bowman)
{
	//qDebug() << "BowmanModel::on_bowman_back_ends_changed()";

	beginResetModel();

	functions_.resize(0);

	//int hal_unique_index_ = 0;
	for(auto back_end_it = bowman->begin(); back_end_it != bowman->end(); ++back_end_it ) {

		for(auto const & function : back_end_it->second ) {
			functions_.emplace_back( FunctionIdentifier{ function.first, back_end_it->first } );
		}

		// json const & j = * back_end_it->second;
		// if( j.count(_f_functions_) ) {
		// 	for(auto const & it : j[_f_functions_] ) {
		// 		if( it.count(_f_name_) ) {
		// 			std::string name = it[_f_name_];
		// 			functions_.emplace_back( FunctionIdentifier{ name, back_end_it->first } );
		// 		}
		// 	}
		// }
	}

	endResetModel();
}



} // namespace network
} // namespace ui
