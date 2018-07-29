#ifndef BOWMAN_VIEW_H
#define BOWMAN_VIEW_H

#include <ui/network/bowman.h>
#include <ui/network/bowman_model.h>

#include <QGroupBox>
#include <QSortFilterProxyModel>

namespace Ui {
class BowmanView;
}


namespace ui {
namespace network {


// class BowmanViewFilter : public QSortFilterProxyModel
// {
//     Q_OBJECT
// public:
// 	using QSortFilterProxyModel::QSortFilterProxyModel;
// protected:
// 	bool filterAcceptsRow(int source_row, const QModelIndex &source_parent) const override;
// };


class BowmanView : public QGroupBox
{
    Q_OBJECT

public:
    explicit BowmanView(QWidget *parent = 0);
    ~BowmanView();


	FunctionID get_select_item();


Q_SIGNALS:
	void double_clicked(ui::network::FunctionID const &);


private Q_SLOTS:
	// void on_bowman_specification_received(std::string const &, JSON_CSP const &);
	// void on_bowman_client_disconnected(std::string const &);

	// void on_filter_textChanged(const QString &text);

	void on_functions_doubleClicked(QModelIndex const & index);

	void when_model_changed();


private:
	ui::network::BowmanModel bowman_model_;
	//BowmanViewFilter filter_;

    Ui::BowmanView *ui;
};

} // namespace network
} // namespace ui


#endif // BOWMAN_VIEW_H
