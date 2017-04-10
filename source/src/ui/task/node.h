#ifndef UI_TASK_NODE_H
#define UI_TASK_NODE_H

#include <QObject>
#include <QUuid>
#include <QString>
#include <QUuid>

#include <memory>

namespace ui {
namespace task {

class Node;
using NodeSP  = std::shared_ptr< Node >;
using NodeCSP = std::shared_ptr< Node const >;

class Node : public QObject
{
    Q_OBJECT

public:
	using Map = std::map< QString, NodeSP >;
	using Key = Map::key_type;

	explicit Node(QUuid _node_id, Node *_parent);

	virtual ~Node() {}


	/// return node type as string
	virtual std::string type() const;

	/// Add new leaf this node, If node was already added to some other Node then remove it first.
    void add(Key const &, NodeSP const &);

	// return number of elemnts stored in this node
	int size() const { return leafs_.size(); }

	// GUI helper function
	// return i'th leaf or nullptr if no such leaf exists
	Node * leaf(int i) const;

	// GUI helper function
	// return key of given leaf or nullptr if leaf could not be found
	Key const * find(Node *leaf) const;

	// erase given leaf from this node, do nothing if given node was not found. return true/false depending if element could be found
    bool erase(Node *leaf);

    Node * parent() const { return parent_; }

	// return index of node, return -1 if node could not be found
	int node_index(Node *node);


protected:
	virtual QByteArray data() const;
	virtual void data(QByteArray const &);


signals:

public slots:



private slots:
    //void data_upload_finished();
    //void data_download_finished();
	//void update_from_json(QJsonObject const &root, bool forced);


private:

	QUuid const node_id;

    Node * parent_ = nullptr;

    Map leafs_;



};

} // namespace task
} // namespace ui

#endif // UI_TASK_NODE_H
