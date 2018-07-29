#include "pose_editor.h"
#include "ui_pose_editor.h"

#include <ui/network/bowman.h>
#include <ui/network/function_setup_dialog.h>


#include <protocols/network/util.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <QFileDialog>
#include <QMessageBox>
#include <QApplication>
#include <QDebug>

using namespace ui::network;
using namespace protocols::network;

namespace ui {
namespace widgets {

using std::string;

quint64 as_number(string const &s)
{
	quint64 r = 0;
	for(unsigned char c : s) r += c*256; // cast to `unsigned char` is intentional here
	return r;
}


PoseEditor::PoseEditor(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::PoseEditor)

{
    ui->setupUi(this);

	Bowman & bowman = Bowman::bowman();

	//connect(&bowman, &Bowman::client_connected,       this, &PoseEditor::client_connected);
	//connect(&bowman, &Bowman::client_disconnected,    this, &PoseEditor::client_disconnected);
	//connect(&bowman, &Bowman::specification_received, this, &PoseEditor::specification_received);
	connect(&bowman, &Bowman::result_received,        this, &PoseEditor::result_received);
	connect(&bowman, &Bowman::progress_data_received, this, &PoseEditor::progress_data_received);

	//connect(ui->functions, &BowmanView::double_clicked, this,  &PoseEditor::on_apply_clicked);

	// ui->movers->setModel(&bowman_model_);
	// connect(&bowman, &Bowman::specification_received, [this]() { bowman_model_.update_from_bowman(bowman); } );
	// bowman_model_.update_from_bowman(bowman);

	//setAttribute(Qt::WA_DeleteOnClose);  created on the stack in pose_viewer/main.cpp main
}

PoseEditor::~PoseEditor()
{
    delete ui;
}


void PoseEditor::on_action_open_pose_triggered()
{
	//qDebug() << "PoseEditor::on_action_open_pose_triggered()";

	core::pose::PoseOP pose = std::make_shared<core::pose::Pose>();

	QString file_name = QFileDialog::getOpenFileName(this, tr("Open Pose"), "", tr("Rosetta Pose (*.pose, *.pdb)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
	if( not file_name.isEmpty() ) {
		QFile file(file_name);

		if (!file.open(QIODevice::ReadOnly) ) return;

		QByteArray data = file.readAll();

		qDebug() << "pose_from_pdbstring...";
		core::import_pose::pose_from_pdbstring(*pose, data.toStdString(), file_name.toUtf8().constData() );
		qDebug() << "pose_from_pdbstring... OK";

		ui->pose->set_pose(pose);
		file_name_ = file_name;

		ui->pose->update_pose_draw();
	}
}

void PoseEditor::on_action_save_pose_triggered()
{
	qDebug() << "PoseEditor::on_action_save_pose_triggered()";
	save_pose(false);
}

void PoseEditor::on_action_save_pose_as_triggered()
{
	qDebug() << "PoseEditor::on_action_save_pose_as_triggered()";
	save_pose(true);
}

bool PoseEditor::save_pose(bool always_ask_for_file_name)
{
	QString file_name = always_ask_for_file_name ? "" : file_name_;

	if( file_name.isEmpty() ) {
		file_name = QFileDialog::getSaveFileName(nullptr, QObject::tr("Save Pose to PDB File"), "", QObject::tr("PDB files (*.pdb)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
	}

	if( not file_name.isEmpty() ) {
		// QFile file(file_name);
		// if( !file.open(QIODevice::WriteOnly) ) return false;
		// QByteArray data = ...;
		// file.write(data)
		// file.close();

		core::io::pdb::dump_pdb( *ui->pose->pose(), file_name.toStdString() );

		file_name_ = file_name;

		return true;
	}
	else return false;
}



// void PoseEditor::client_connected(std::string const & /*id*/)
// {
// 	//qDebug() << "PoseEditor::client_connected(): " << id;
// 	// basically ignore until we got specification (see specification_received)
// }

// void PoseEditor::client_disconnected(std::string const & _id)
// {
// 	auto id = as_number(_id);
// 	qDebug() << "PoseEditor::client_disconnected(): " << id;

// 	auto it = back_ends_.find(_id);

// 	if( it != back_ends_.end() ) {
// 		auto i = ui->back_ends->indexOf(it->second);
// 		if( i >= 0 ) ui->back_ends->removeTab(i);

// 		delete it->second;
// 		back_ends_.erase(it);
// 	}
// }

// void PoseEditor::specification_received(std::string const & _id, JSON_CSP const &_j)
// {
// 	auto id = as_number(_id);
// 	json const & j = *_j;

// 	qDebug() << "PoseEditor::specification_received(): id=" << id << "\n" << j.dump(2).c_str(); //QString::fromStdString(j->dump(2));

// 	auto movers = new QListWidget();

// 	//new QListWidgetItem("Test", movers);
// 	if( j.count(_f_functions_) ) {
// 		for(auto const & it : j[_f_functions_] ) {
// 			if( it.count(_f_name_) ) {
// 				std::string name = it[_f_name_];
// 				new QListWidgetItem(name.c_str(), movers);
// 			}
// 		}
// 	}

// 	ui->back_ends->addTab(movers, QString::number(hal_unique_index_++) );

// 	back_ends_.emplace(_id, movers);
// }

void PoseEditor::result_received(core::pose::PoseOP const &pose, JSON_CSP const &/*_result*/)
{
	qDebug() << "PoseEditor::result_received():";// << result.dump(2).c_str(); //QString::fromStdString(j->dump(2));
	//qDebug() << "PoseEditor::result_received():" << result["meta"].dump(2).c_str(); //QString::fromStdString(j->dump(2));
	//qDebug() << "PoseEditor::result_received(): result[\"bin\"]:" << as_hexadecimal( result["bin"], true ).c_str();

	// json const & result = *_result;
	// auto it_pose = result.find(_f_pose_);
	// if( it_pose != result.end() ) {
	// 	core::pose::PoseOP pose = protocols::network::bytes_to_pose( it_pose.value() );
	// 	//std::cout << "received pose: " << *pose << std::endl;
	// 	ui->pose->set_pose(pose);
	// 	ui->pose->update_pose_draw();
	// }

	if(pose) {
        ui->pose->set_pose(pose);
		ui->pose->update_pose_draw();
	}
}

void PoseEditor::progress_data_received(core::pose::PoseOP const &pose, JSON_CSP const &/*_iresult*/)
{
	//qDebug() << "PoseEditor::progress_data_received():";// << result.dump(2).c_str(); //QString::fromStdString(j->dump(2));

	// json const & result = *_iresult;
	// auto it_pose = result.find(_f_pose_);
	// if( it_pose != result.end() ) {
	// 	core::pose::PoseOP pose = protocols::network::bytes_to_pose( it_pose.value() );
	// 	//std::cout << "received pose: " << *pose << std::endl;
	// 	ui->pose->set_pose(pose);
	// 	ui->pose->update_pose_draw();
	// }

	if(pose) {
		ui->pose->set_pose(pose);
        // First check if there's a FMI maybe
        //ui->pose->set_colour_mode( ui::ui_core::pose_draw::SPDOGLW_colour_by_input_domain ) ;
		ui->pose->update_pose_draw();
	}
}

void PoseEditor::on_apply_clicked()
{
	Q_EMIT on_functions_double_clicked( ui->functions->get_select_item() );
}


void PoseEditor::on_functions_double_clicked(FunctionID const &fid)
{
	auto pose = ui->pose->pose();

	//QString qfunction_name = ui->functions->get_select_item();
	//string function_name = qfunction_name.toStdString();
	ArgumentsSP args = network::get_function_arguments(Bowman::bowman(), fid);

	if(args) {
		if( not pose ) {
			if( args->find(_f_pose_) != args->end() ) {
				QMessageBox box(QMessageBox::Warning, "Structure is Empty!", QString("Current structure is Empty!\n\nFunction `%1`\nrequire non-empty Pose object\n\nPlease load the Pose first!").arg(fid.name.c_str()) );
				box.exec();
				return;
			}
		}

		FunctionSetupDialog setup(fid.name.c_str(), args, this);
		if( setup.exec() == QDialog::Accepted ) Bowman::bowman().execute( ui->pose->pose(), fid, *args);
	}
	else {
		qDebug() << "PoseEditor::on_functions_double_clicked(...): supplied function_name is empty!";
		QApplication::beep();
	}
}

void PoseEditor::on_abort_clicked()
{
	QMessageBox box(QMessageBox::Critical, "Abort HAL Call?", "This will terminate Rosetta run.\nAre you sure you want to proceed?", QMessageBox::Yes | QMessageBox::No);
	box.setDefaultButton(QMessageBox::No);
	if( box.exec() == QMessageBox::Yes ) Bowman::bowman().abort();
}

void PoseEditor::on_color_toggled(bool state)
{
	//qDebug() << "PoseEditor::on_color_toggled(" << state << ")";
    ui->pose->set_color_mode(state ? ui_core::pose_draw::SimplePoseDrawOpenGLWidget::ColorMode::input_domain : ui_core::pose_draw::SimplePoseDrawOpenGLWidget::ColorMode::none );
	ui->pose->update_pose_draw();
}


void PoseEditor::on_pause_toggled(bool state)
{
	Bowman::bowman().pause(state);

	QPalette pal = ui->pause->palette();

	if(state) {
		ui->pause->setText("paused");

		pal.setColor(QPalette::ButtonText, QColor(Qt::red) );
		//pal.setColor(QPalette::WindowText,  QColor(Qt::yellow) );
		// ui->pause->setAutoFillBackground(true);
		//ui->pause->update();

		//
		// pal.setColor(QPalette::Button, QColor(Qt::blue));
		// ui->pause->setAutoFillBackground(true);
		// ui->pause->setPalette(pal);
		// ui->pause->update();
		//ui->pause->setText("resume");
	}
	else {
		pal.setColor(QPalette::ButtonText, QColor(Qt::black) );
		ui->pause->setText("pause");
		//ui->pause->setAutoFillBackground(false);
	}
	ui->pause->setPalette(pal);

}

} // namespace widgets
} // namespace ui
