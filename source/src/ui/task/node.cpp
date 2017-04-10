#include <ui/task/node.h>


namespace ui {
namespace task {

Node::Node(QUuid _node_id, Node *_parent) : node_id(_node_id), parent_(_parent)
{

}

std::string Node::type() const
{
	return "node"; // type for generic node that does act as simple container with no extra functionality
}


QByteArray Node::data() const
{
	return QByteArray();
}


void Node::data(QByteArray const &)
{
}


void Node::add(Key const &key, NodeSP const &node)
{
	if( node->parent_ ) node->parent_->erase( node.get() );

	leafs_[key] = node;
    node->parent_ = this;
}


// GUI helper function
// return i'th leaf or nullptr if no such leaf exists
Node * Node::leaf(int i) const
{
	if( i < 0  or  i >= size() ) return nullptr;

    return std::next( leafs_.begin(), i)->second.get();
}

// GUI helper function
// return key of given leaf or nullptr if leaf could not be found
Node::Key const * Node::find(Node *leaf) const
{
    auto it = find_if(leafs_.begin(), leafs_.end(), [&](Map::value_type const &p) { return p.second.get() == leaf; } );

	if( it == leafs_.end() ) return nullptr;
	else return &it->first;
}


// erase given lead from this node
bool Node::erase(Node *leaf)
{
    auto it = find_if(leafs_.begin(), leafs_.end(), [&](Map::value_type const &p) { return p.second.get() == leaf; } );

	if( it == leafs_.end() ) return false;

    leafs_.erase(it);
   	return true;
}

// return index of node, return -1 if node could not be found
int Node::node_index(Node *node)
{
	auto it = find_if(leafs_.begin(), leafs_.end(), [&](Map::value_type const &p) { return p.second.get() == node; } );

	if( it != leafs_.end() ) return distance(leafs_.begin(), it);
	else return -1;
}


} // namespace task
} // namespace ui
