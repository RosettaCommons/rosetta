#include "bowman_model.h"

#include <protocols/network/util.hh>

#include <QDebug>

using namespace protocols::network;

using std::string;
using std::vector;


namespace ui {
namespace network {

/// struct to help represent BowmanModel as a tree. This is local datatype used only by BowmanModel
// struct BowmanModel::Node
// {
// 	BowmanModelNode *parent = nullptr;

// 	QString data;

// 	vector <PNodeSP> leafs;


// 	PNode(Project const &p) : name("root"), type("project") {
// 		for(auto it=p.tasks_.begin(); it!=p.tasks_.end(); ++it) leafs.push_back( std::make_shared<PNode>(this, it->first, it->second) );
// 	}

// 	PNode(PNode *_parent, QString const &_name, TaskSP const &t) : parent(_parent), name(_name), type("task"), task(t) {
// 		if( !t->input() .empty() ) leafs.push_back( std::make_shared<PNode>(this, "input",  t, static_cast<File const & (Task::*)() const>( &Task::input)  ) );
// 		if( !t->script().empty() ) leafs.push_back( std::make_shared<PNode>(this, "script", t, static_cast<File const & (Task::*)() const>( &Task::script) ) );
// 		if( !t->flags() .empty() ) leafs.push_back( std::make_shared<PNode>(this, "flags",  t, static_cast<File const & (Task::*)() const>( &Task::flags)  ) );
// 	}

// 	PNode(PNode *_parent, QString const &_name, TaskSP const &t, File const & (Task::*/* f */)() const) : parent(_parent), name(_name), type("file"), task(t) {
// 	}


// 	// return i'th leaf or nullptr if no such leaf exists
// 	PNode * leaf(unsigned int i) const {
// 		if( i < leafs.size() ) return leafs[i].get();
// 		return nullptr;
// 	}

// 	// return index of node, return -1 if node could not be found
// 	int node_index(PNode *node) {
// 		for(unsigned int i=0; i<leafs.size(); ++i) {
// 			if( leafs[i].get() == node ) return i;
// 		}

// 		return -1;
// 	}

// 	int size() { return leafs.size(); }
// };


// Node *ProjectModel::Node::get_node(const QModelIndex &index) const
// {
//     if (index.isValid()) {
//         Node *item = static_cast<Node*>(index.internalPointer());
//         if (item) return item;
//     }
//     return &root;
// }




BowmanModel::BowmanModel(QObject *parent) : util::TreeNodeModel<BowmanModelNodeData>(parent)  //QAbstractItemModel(parent)
{
	connect(&Bowman::bowman(), &Bowman::back_ends_changed, this, &BowmanModel::on_bowman_back_ends_changed);

	// using Node = ModelBase::Node;
	// auto b1 = std::make_shared<Node>();  b1->name = "back-end-1";
	// auto b2 = std::make_shared<Node>();  b2->name = "back-end-2";
	// auto b1_1 = std::make_shared<Node>();  b1_1->name = "1-f1";
	// auto b1_2 = std::make_shared<Node>();  b1_2->name = "1-f2";
	// b1->push_back(b1_1);
	// b1->push_back(b1_2);
	// auto b2_1 = std::make_shared<Node>();  b2_1->name = "2-f1";
	// b2->push_back(b2_1);
	// root_.push_back(b1);
	// root_.push_back(b2);
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




// QModelIndex BowmanModel::index(int row, int column, QModelIndex const & parent) const
// {
// 	return ModelBase::index(row, column, parent);
// }


// QModelIndex BowmanModel::parent(const QModelIndex &index) const
// {
// 	return ModelBase::parent(index);
// }

// int BowmanModel::rowCount(const QModelIndex &parent) const
// {
// 	return ModelBase::rowCount(parent);
// }

int BowmanModel::columnCount(const QModelIndex &/*parent*/) const
{
    //if( parent.isValid() ) return 0;

	return 1;
}

QVariant BowmanModel::data(const QModelIndex &index, int role) const
{
    if( !index.isValid() ) return QVariant();

	if (role != Qt::DisplayRole ) return QVariant();

	//IndexPack i = index.internalId();
	if( auto node = get_node(index) ) {
		return QString::fromStdString( node->name );
	}

	return "";

	// int row = index.row();
	// if( row < static_cast<int>( functions_.size() )  and  row >= 0 ) {
	// 	if      ( index.column() == 0 ) return QString::fromStdString(functions_[row].name);
	// 	else if ( index.column() == 1 ) return QString::fromStdString(functions_[row].hal_id);
	// }

    // return QVariant();
}

/// return FunctionID for given index, return nullptr if index is invalid
FunctionID BowmanModel::function_id(QModelIndex const &index)
{
	if( index.isValid() ) {

		auto node = get_node(index);
		if( node  and  node->parent != &root_ ) return *node;  //qDebug() << "BowmanModel::function_id():" << node->hal_id.c_str() << node->name.c_str();// << " parent.hal_id:" << node->parent->hal_id.c_str(); }

		// if ( index.row() < static_cast<int>( node->leafs.size() ) ) {
		// 	return *node->leafs[index.row()];
		// 	//and  index.row() < static_cast<int>( functions_.size() ) ) return &functions_[index.row()];
		// }
	}
	return FunctionID();
}

void BowmanModel::on_bowman_back_ends_changed(Bowman const *bowman)
{
	//qDebug() << "BowmanModel::on_bowman_back_ends_changed()";

	beginResetModel();

	root_.leafs.resize(0);

	using Node = ModelBase::Node;

	string filter = filter_.toStdString();

	for(auto back_end_it = bowman->begin(); back_end_it != bowman->end(); ++back_end_it ) {

		auto b = std::make_shared<Node>();
		b->name = back_end_it->second.name;
		b->hal_id = back_end_it->first;

		for(auto const & function : back_end_it->second.functions ) {
			//functions_.emplace_back( FunctionID{ function.first, back_end_it->first } );


			//if ( function.first.find(filter) != std::string::npos ) {

			auto pos = std::search(function.first.begin(), function.first.end(), filter.begin(), filter.end(), [](char l, char r) { return (std::tolower(l) == std::tolower(r)); });
			if ( pos != function.first.end() ) {

				auto f = std::make_shared<Node>();
				f->name = function.first;
				f->hal_id = back_end_it->first;

				b->push_back(f);

				//qDebug() << "BowmanModel::on_bowman_back_ends_changed(): " << function.first.c_str();
			}

		}

		if( not b->leafs.empty() ) root_.push_back(b);


		// json const & j = * back_end_it->second;
		// if( j.count(_f_functions_) ) {
		// 	for(auto const & it : j[_f_functions_] ) {
		// 		if( it.count(_f_name_) ) {
		// 			std::string name = it[_f_name_];
		// 			functions_.emplace_back( FunctionID{ name, back_end_it->first } );
		// 		}
		// 	}
		// }
	}
	// auto b1 = std::make_shared<Node>();  b1->label = "back-end-1";
	// auto b2 = std::make_shared<Node>();  b2->label = "back-end-2";
	// auto b1_1 = std::make_shared<Node>();  b1_1->label = "1-f1";
	// auto b1_2 = std::make_shared<Node>();  b1_2->label = "1-f2";
	// b1->push_back(b1_1);
	// b1->push_back(b1_2);
	// auto b2_1 = std::make_shared<Node>();  b2_1->label = "2-f1";
	// b2->push_back(b2_1);
	// root_.push_back(b1);
	// root_.push_back(b2);

	endResetModel();
}


void BowmanModel::set_filter(QString const &filter)
{
	filter_ = filter;

	on_bowman_back_ends_changed( &Bowman::bowman() );
}



} // namespace network
} // namespace ui
