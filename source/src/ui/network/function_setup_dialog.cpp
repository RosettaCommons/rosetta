#include "function_setup_dialog.h"
#include "ui_function_setup_dialog.h"

#include <ui/network/argument.h>
#include <protocols/network/util.hh>

#include <QLabel>
#include <QFormLayout>
#include <QValidator>
#include <QLineEdit>

#include <QDebug>

using namespace protocols::network;


namespace ui {
namespace network {

FunctionSetupDialog::FunctionSetupDialog(QString const & function_name, ArgumentsSP const &args, QWidget *parent) :
    QDialog(parent),
    ui(new Ui::FunctionSetupDialog)
{
    ui->setupUi(this);

	setWindowTitle(function_name + " setup");

	//qDebug() << "ui->scrollArea->widget():" << ui->scrollArea->widget();

	auto w = ui->scrollArea->widget();

	QFormLayout *layout = new QFormLayout;

	if( args->empty()  or  (args->size() == 1  and  args->count(_f_pose_) == 1 ) ) {
		auto l = new QLabel( QString("function %1 has no arguments\nclick \"OK\" button to start execution...").arg(function_name), w);
		layout->addWidget(l);
		setMinimumWidth(424);
	}
	else {
		//w->setMinimumHeight(1024);
		//w->setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::Expanding);

		for(auto const & a : *args) {
			if( a.first != protocols::network::_f_pose_ ) {
				auto name = new QLabel( QString("%1:").arg( a.first.c_str() ), w);
				auto input = a.second->create_input_widget(w);
				input->setToolTip( QString::fromStdString(a.second->description() ) );
				layout->addRow(name, input);
			}
		}
	}
	// for(int i=0; i<100; ++i) {
	// 	auto l1 = new QLabel("label-", w);
	// 	layout->addWidget(l1);
	// }
	// QValidator *validator = new QIntValidator(100, 999, this);
	// QLineEdit *edit = new QLineEdit(this);
	// edit->setValidator(validator);

	w->setLayout(layout);

	adjustSize();  //setSizeGripEnabled(true); <-- does not work on Mac
}

FunctionSetupDialog::~FunctionSetupDialog()
{
    delete ui;
}

void FunctionSetupDialog::ok_clicked()
{
	//qDebug() << "FunctionSetupDialog::on_ok_clicked";

	// todo: add validation code here
	Q_EMIT accept();
}


} // namespace network
} // namespace ui
