#pragma once

#include <QModelIndex>
#include <QAbstractItemModel>

#include <vector>

namespace ui {
namespace util {

///
/// auxiliary class to generate tree-like structures as a helper when creating QAbstractItemModel-like types
///
/// Usage:
///   - define struct with data memebers that will be hold on each node: struct MyData { int a, b, c; };
///   - and then use TreeNodeModel<MyData>::Node as base for your tree structure
///
/// Example of usage could be seen in: ui/network/bowman_model.h
///
template <typename T>
struct TreeNodeModel : public QAbstractItemModel
{

public:

	explicit TreeNodeModel(QObject *parent = nullptr) : QAbstractItemModel(parent) {};

	//template <typename T>
	struct Node : public T
	{
		using NodeSP  = std::shared_ptr< Node >;

		Node *parent = nullptr;

		std::vector <NodeSP> leafs;

		// return i'th leaf or nullptr if no such leaf exists
		Node * leaf(unsigned int i) const;

		// return index of node, return -1 if node could not be found
		int node_index(Node *node);

		void push_back(NodeSP const &);
	};


	mutable Node root_;  // mutable due to QAbstractItemModel usage of `void *` instead of `void const *`

	Node *get_node(QModelIndex const & index) const;

	QModelIndex index(int row, int column, const QModelIndex &parent) const override;
	QModelIndex parent(const QModelIndex &index) const override;

	int rowCount(QModelIndex const &parent) const override;
};


// return i'th leaf or nullptr if no such leaf exists
template <typename T>
typename TreeNodeModel<T>::Node * TreeNodeModel<T>::Node::leaf(unsigned int i) const
{
	if( i < leafs.size() ) return leafs[i].get();
	return nullptr;
}


// return index of node, return -1 if node could not be found
template <typename T>
int TreeNodeModel<T>::Node::node_index(Node *node)
{
	for(unsigned int i=0; i<leafs.size(); ++i) {
		if( leafs[i].get() == node ) return i;
	}
	return -1;
}


template <typename T>
void TreeNodeModel<T>::Node::push_back(NodeSP const &n)
{
	n->parent = this;
	leafs.push_back(n);
}



template <typename T>
typename TreeNodeModel<T>::Node * TreeNodeModel<T>::get_node(QModelIndex const & index) const
{
	if ( index.isValid() ) return static_cast<Node*>( index.internalPointer() );
    else return &root_;


	// if ( index.isValid() ) {
	// 	Node *node = static_cast<Node*>(index.internalPointer());
    //     if (node) return node;
	// }
    // return &root_;
}


template <typename T>
QModelIndex TreeNodeModel<T>::index(int row, int column, const QModelIndex &parent) const
{
    if (parent.isValid() and parent.column() != 0 ) return QModelIndex();

	Node *parent_node = get_node(parent);

    Node *leaf = parent_node->leaf(row);

    if (leaf) return createIndex(row, column, leaf);
    else return QModelIndex();
}


template <typename T>
QModelIndex TreeNodeModel<T>::parent(const QModelIndex &index) const
{
    if (!index.isValid()) return QModelIndex();

    Node *leaf = get_node(index);
    Node *parent = leaf->parent;

    if( parent == nullptr  or  leaf == nullptr ) return QModelIndex();

    return createIndex(parent->node_index(leaf), 0, parent);
}


template <typename T>
int TreeNodeModel<T>::rowCount(QModelIndex const &parent_index) const
{
	if( auto parent = get_node(parent_index) ) return parent->leafs.size();
	else return 0;
}


// int TreeNodeModel<T>::columnCount(QModelIndex const &parent) const
// {
// }



} // namespace util
} // namespace ui
