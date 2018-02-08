#include "config_dialog.h"
#include "ui_config_dialog.h"

#include <ui/config/util.h>

#include <QFileDialog>
#include <QDebug>

namespace ui {
namespace config {

ConfigDialog::ConfigDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::ConfigDialog)
{
    ui->setupUi(this);

	update_ui_from_user_settings();

	connect(this, SIGNAL(accepted()), this, SLOT(update_user_settings_from_ui()));
}

ConfigDialog::~ConfigDialog()
{
    delete ui;
}


void ConfigDialog::update_ui_from_user_settings()
{
	UserCredentials c = get_user_credentials();

	ui->user_name->setText(c.user);
	ui->password->setText(c.password);

	ui->pdb_viewer->setText( get_pdb_viewer_path() );
}

void ConfigDialog::update_user_settings_from_ui()
{
	qDebug() << "ConfigDialog::update_user_settings_from_ui()";

	UserCredentials c;
	c.user = ui->user_name->text();
	c.password = ui->password->text();

	set_user_credentials(c);

	set_pdb_viewer_path( ui->pdb_viewer->text() );
}

void ConfigDialog::on_pdb_viewer_browse_clicked()
{
	//qDebug() << "ConfigDialog::on_pdb_viewer_browse_clicked()";

	ui->pdb_viewer->setText( QFileDialog::getOpenFileName(this, tr("Select PDB viewer App") ) );
}




} // namespace config
} // namespace ui
