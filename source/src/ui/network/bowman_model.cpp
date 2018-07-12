#include "bowman_model.h"

#include <protocols/network/util.hh>

#include <QDebug>

using namespace protocols::network;


namespace ui {
namespace network {

BowmanModel::BowmanModel(QObject *parent) : QAbstractTableModel(parent)
{
	connect(&Bowman::bowman(), &Bowman::back_ends_changed, this, &BowmanModel::on_bowman_back_ends_changed);
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
	qDebug() << "BowmanModel::on_bowman_back_ends_changed()";


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
