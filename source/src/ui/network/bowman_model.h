#ifndef BOWMANMODEL_H
#define BOWMANMODEL_H

#include <ui/network/bowman.h>


#include <QAbstractTableModel>

namespace ui {
namespace network {

class BowmanModel : public QAbstractTableModel
{
    Q_OBJECT

public:
    explicit BowmanModel(QObject *parent = nullptr);

    // Header:
    QVariant headerData(int section, Qt::Orientation orientation, int role = Qt::DisplayRole) const override;

    // Basic functionality:
    int rowCount(const QModelIndex &parent = QModelIndex()) const override;
    int columnCount(const QModelIndex &parent = QModelIndex()) const override;

    QVariant data(const QModelIndex &index, int role = Qt::DisplayRole) const override;

	struct FunctionIdentifier {
		std::string name, hal_id;
	};

	/// return FunctionIdentifier for given index, return nullptr if index is invalid
	FunctionIdentifier * get_identifier(QModelIndex const &index);

	void update_from_bowman(Bowman const &bowman);

private:

	std::vector<FunctionIdentifier> functions_;
};

} // namespace network
} // namespace ui

#endif // BOWMANMODEL_H
