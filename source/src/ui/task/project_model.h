#ifndef UI_TASK_PROJECTMODEL_H
#define UI_TASK_PROJECTMODEL_H

#include <QAbstractItemModel>

#include <ui/task/project.h>


namespace ui {
namespace task {


class ProjectModel : public QAbstractItemModel
{
    Q_OBJECT

public:
    explicit ProjectModel(QObject *parent = 0);

    // Header:
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;

    // Basic functionality:
    QModelIndex index(int row, int column, QModelIndex const &parent = QModelIndex()) const override;
    QModelIndex parent(const QModelIndex &index) const override;

    int rowCount(const QModelIndex &parent = QModelIndex()) const override;
    int columnCount(const QModelIndex &parent = QModelIndex()) const override;

    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;

	Qt::ItemFlags flags(const QModelIndex &index) const override;


	Node * root() const { return root_; }



private:
	Node *get_item(const QModelIndex &index) const;

	Project* root_ = nullptr;

};

} // namespace task
} // namespace ui

#endif // UI_TASK_PROJECTMODEL_H
