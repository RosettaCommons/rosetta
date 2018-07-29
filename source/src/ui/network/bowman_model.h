#ifndef BOWMANMODEL_H
#define BOWMANMODEL_H

#include <ui/network/bowman.h>


#include <ui/util/tree_node.h>

#include <QAbstractItemModel>

namespace ui {
namespace network {

// struct BowmanModelNodeData {
// 	QString label;
// 	std::string hal_id;
// };

using BowmanModelNodeData = FunctionID;

class BowmanModel : public util::TreeNodeModel<BowmanModelNodeData>
{
    Q_OBJECT

	using ModelBase = util::TreeNodeModel<BowmanModelNodeData>;

public:
    explicit BowmanModel(QObject *parent = nullptr);

    // Header:
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;


    // Basic functionality:

	//QModelIndex index(int row, int column, QModelIndex const & parent = QModelIndex()) const override;
	//QModelIndex parent(const QModelIndex &index) const override;

    //int rowCount(const QModelIndex &parent = QModelIndex()) const override;
    int columnCount(const QModelIndex &parent = QModelIndex()) const override;

    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;

	/// return FunctionIdentifier for given index, return empty identifier if index is invalid
	FunctionID function_id(QModelIndex const &index);

public Q_SLOTS:
	void set_filter(QString const &);

private Q_SLOTS:
	void on_bowman_back_ends_changed(Bowman const *bowman);

private:

	//std::vector<FunctionIdentifier> functions_;


private:
	QString filter_;

	//using Node = util::TreeNode<NodeData>;

	//Node root_;
};

} // namespace network
} // namespace ui

#endif // BOWMANMODEL_H
