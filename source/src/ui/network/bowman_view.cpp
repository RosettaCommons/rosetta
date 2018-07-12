#include "bowman_view.h"
#include "ui_bowman_view.h"

#include <QDebug>

namespace ui {
namespace network {


BowmanView::BowmanView(QWidget *parent) :
    QGroupBox(parent),
    ui(new Ui::BowmanView)
{
    ui->setupUi(this);

	// Bowman & bowman = Bowman::bowman();

	//connect(&bowman, &Bowman::client_connected,       this, &BowmanView::client_connected);
	// connect(&bowman, &Bowman::client_disconnected,    this, &BowmanView::client_disconnected);
	// connect(&bowman, &Bowman::specification_received, this, &BowmanView::specification_received);
	// connect(&bowman, &Bowman::result_received,        this, &BowmanView::result_received);
	// connect(&bowman, &Bowman::progress_data_received, this, &BowmanView::progress_data_received);

	filter_.setSourceModel( &bowman_model_ );
	ui->functions->setModel( &filter_ );

	//setAttribute(Qt::WA_DeleteOnClose);
}

BowmanView::~BowmanView()
{
	//qDebug() << "BowmanView::~BowmanView()";
    delete ui;
}


void BowmanView::on_filter_textChanged(const QString &text)
{
	//qDebug() << "BowmanView::on_filter_textChanged(...)";

	filter_.setFilterRegExp(QRegExp(text, Qt::CaseInsensitive, QRegExp::FixedString));
	filter_.setFilterKeyColumn(0);
}


QString BowmanView::get_select_item()
{
	// QVariant qv = model->data(index, Qt::DisplayRole);
	// if( qv.isValid() ) {
	// 	index.isValid()

	QModelIndex index = ui->functions->currentIndex();
	QString name = index.data(Qt::DisplayRole).toString();

	return name;
}


void BowmanView::on_functions_doubleClicked(const QModelIndex &index)
{
	//qDebug() << "BowmanView::on_functions_doubleClicked()";

	QString name = index.data(Qt::DisplayRole).toString();

	if( not name.isEmpty() ) Q_EMIT double_clicked(name);
}


// void BowmanView::on_bowman_specification_received(std::string const &, JSON_CSP const &)
// {
// 	bowman_model_.update_from_bowman()
// }
// void BowmanView::on_bowman_client_disconnected(std::string const &)
// {
// }


} // namespace network
} // namespace ui
