#include "bowman_model.h"

#include <QDebug>

namespace ui {
namespace network {

BowmanModel::BowmanModel(QObject *parent) : QAbstractTableModel(parent)
{
}

QVariant BowmanModel::headerData(int section, Qt::Orientation orientation, int role) const
{
	if (role == Qt::DisplayRole) {
		if (orientation == Qt::Horizontal) {
			switch (section) {
				case 0: return QString("mover");
				case 1: return QString("back-end");
			}
		} else if (orientation == Qt::Vertical) {
			return QString::number(section+1);
		}
	}
    return QVariant();
}

int BowmanModel::rowCount(const QModelIndex &parent) const
{
    if (parent.isValid())
        return 0;

    return functions_.size();
}

int BowmanModel::columnCount(const QModelIndex &parent) const
{
    if( parent.isValid() ) return 0;

	return 2;
}

QVariant BowmanModel::data(const QModelIndex &index, int role) const
{
    if( !index.isValid() ) return QVariant();

	if (role != Qt::DisplayRole ) return QVariant();


	if( index.row() < static_cast<int>( functions_.size() ) ) {
		if      ( index.column() == 0 ) return QString("zero-") + QString::number( index.row() );
		else if ( index.column() == 1 ) return QString("one-") + QString::number( index.row() );
	}

    return QVariant();
}

/// return FunctionIdentifier for given index, return nullptr if index is invalid
BowmanModel::FunctionIdentifier * BowmanModel::get_identifier(QModelIndex const &index)
{
	if( index.isValid()  and  index.row() < static_cast<int>( functions_.size() ) ) return &functions_[index.row()];
	else return nullptr;
}

void BowmanModel::update_from_bowman(Bowman const &bowman)
{
	qDebug() << "BowmanModel::update_from_bowman():" << bowman.size();

	beginResetModel();

	functions_.resize(0);

	for(auto it = bowman.begin(); it != bowman.end(); ++it ) {
		functions_.emplace_back( FunctionIdentifier{"name", it->first} );
	}

	endResetModel();
}



} // namespace network
} // namespace ui
