#include "bowman_view.h"
#include "ui_bowman_view.h"

#include <QTimer>
#include <QDebug>

namespace ui {
namespace network {


// bool BowmanViewFilter::filterAcceptsRow(int source_row, const QModelIndex &source_parent) const
// {
// 	if( not source_parent.isValid() ) return true; // always accept tree roots
// 	return QSortFilterProxyModel::filterAcceptsRow(source_row, source_parent);
// }


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


	connect(ui->filter, &QLineEdit::textChanged, &bowman_model_, &BowmanModel::set_filter);

	//connect(&bowman_model_, &BowmanModel::modelReset, ui->functions, &QTreeView::expandAll);
	connect(&bowman_model_, &BowmanModel::modelReset, this, &BowmanView::when_model_changed);


	//filter_.setSourceModel( &bowman_model_ );
	//ui->functions->setModel( &filter_ );
	ui->functions->setModel( &bowman_model_ );

	ui->functions->expandAll();

	//setAttribute(Qt::WA_DeleteOnClose);
}

BowmanView::~BowmanView()
{
	//qDebug() << "BowmanView::~BowmanView()";
    delete ui;
}


// void BowmanView::on_filter_textChanged(const QString &text)
// {
// 	//qDebug() << "BowmanView::on_filter_textChanged(...)";
// 	filter_.setFilterRegExp(QRegExp(text, Qt::CaseInsensitive, QRegExp::FixedString));
// 	filter_.setFilterKeyColumn(0);
// }


void BowmanView::when_model_changed()
{
	//qDebug() << "BowmanView::on_action_model_changed";
	QTimer::singleShot(500, ui->functions, SLOT(expandAll()));
}



FunctionID BowmanView::get_select_item()
{
	//fix me
	// QVariant qv = model->data(index, Qt::DisplayRole);
	// if( qv.isValid() ) {
	// 	index.isValid()

	// QModelIndex index = ui->functions->currentIndex();
	// QString name = index.data(Qt::DisplayRole).toString();


	return bowman_model_.function_id( ui->functions->currentIndex() );
}


void BowmanView::on_functions_doubleClicked(const QModelIndex &index)
{
	//qDebug() << "BowmanView::on_functions_doubleClicked()";

	auto fi = bowman_model_.function_id(index);

	//qDebug() << "BowmanView::on_functions_doubleClicked():" << fi.name.c_str() << as_hexadecimal(fi.hal_id).c_str();

	if ( not fi.name.empty() ) Q_EMIT double_clicked(fi);
	//QString name = index.data(Qt::DisplayRole).toString();
	//if( not name.isEmpty() ) Q_EMIT double_clicked(name);
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
