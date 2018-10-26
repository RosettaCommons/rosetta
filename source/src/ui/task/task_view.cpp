#include <ui/task/task_view.h>
#include "ui_task_view.h"

#include <ui/task/task_submit.h>
#include <ui/task/job_view.h>
#include <ui/task/task.h>
#include <ui/task/file.h>
#include <ui/task/task_cancel_dialog.h>

#include <ui/config/util.h>
#include <ui/config/config_dialog.h>
#include <ui/viewers/score_file_view.h>

#include <QFileDialog>
#include <QStringListModel>
#include <QProcess>
#include <QMenu>
#include <QJsonObject>
#include <QJsonArray>
#include <QMessageBox>
#include <QTemporaryFile>

namespace ui {
namespace task {

TaskView::TaskView(TaskSP const &task,QWidget *parent) :
	QWidget(parent),
	ui(new Ui::TaskView),
	task_(task),
	previewed_file_(FileID::Kind::none, "")
{
    ui->setupUi(this);

	auto model = new QStringListModel(this);
	ui->output->setModel(model);

	//bool old = task->blockSignals(true); // we now handling this on Model (ie Task layer) by ommiting change() emiition if data is the same
	update_ui_from_task();
	//task->blockSignals(old);

	connect(task_.get(), &Task::changed,           this, &TaskView::update_ui_from_task);
	connect(task_.get(), &Task::file_list_changed, this, &TaskView::update_ui_file_list_from_task);
	connect(task_.get(), &Task::syncing,           this, &TaskView::update_syncing_progress);
	connect(task_.get(), &Task::file_changed,      this, &TaskView::file_changed);

	ui->output->setContextMenuPolicy(Qt::CustomContextMenu);
	connect(ui->output, SIGNAL(customContextMenuRequested(const QPoint &)), this, SLOT(create_output_context_menu(const QPoint &)));

	for(auto const & job : task_->jobs()) {
		auto jv = new JobView(job, task_);
		ui->jobs->addTab(jv, job->name);
	}

	setAttribute(Qt::WA_DeleteOnClose);
}

TaskView::~TaskView()
{
    //qDebug() << "TaskView::~TaskView";
    delete ui;
}

void TaskView::update_ui_from_task()
{
    //qDebug() << "TaskView::update_ui_from_task";

	ui->name->setText( task_->name() );
	ui->version->setText( task_->version() );

	ui->queue->setText( task_->queue() );

	ui->state->setText( Task::to_string( task_->state() ) );

	ui->description->document()->setPlainText( task_->description() );

	ui->task_id->setText( task_->task_id() );

	// ui->input_file_name->setText( task_->input().file_name() );
	//ui->input->document()->setPlainText( task_->input().data() );

	// ui->script_file_name->setText( task_->script().file_name() );
	// ui->script->document()->setPlainText( task_->script().data() );

	// if( task_->input().data().size() and (not ui->input_preview->pose()) ) {
	// 	core::pose::PoseOP pose = std::make_shared<core::pose::Pose>();
	// 	core::import_pose::pose_from_pdbstring(*pose, std::string( task_->input().data().constData(), task_->input().data().length() ), task_->input().file_name().toUtf8().constData() );
	// 	ui->input_preview->set_pose(pose);
	// }

	// QWidget * inputs[] = {ui->name, ui->description, ui->input_file_name, ui->script_file_name, ui->script, ui->flags_file_name, ui->flags, ui->nstruct};

	// if( task_->state() == Task::State::_draft_  and (not syncing) ) {
	// 	ui->submit->show();

	// 	for(auto w : inputs) w->setEnabled(true);
	// }
	// else {
	// 	ui->submit->hide();

	// 	for(auto w : inputs) w->setEnabled(false);
	// }

	update_ui_file_list_from_task();
	update_syncing_progress();
}

void TaskView::update_ui_file_list_from_task()
{
	if( auto model = qobject_cast<QStringListModel*>( ui->output->model() ) ) {
		QStringList output_files;
		for(auto const &file_sp : task_->files() ) {
			//qDebug() << "TaskView::update_ui_from_task adding output file: " << k_fsp.first;
			//if( not k_fsp.second->data().isEmpty() ) output_files << k_fsp.first;
			output_files << file_sp->name();
		}
		model->setStringList(output_files);
	} else {
		qDebug() << "TaskView::update_ui_from_task: no QStringListModel!!!";
	}

	ui->files_count->setText( QString::number( task_->files().size() ) );
}

void TaskView::update_syncing_progress()
{
    //qDebug() << "TaskView::update_syncing_progress()...";

	if( task_->is_syncing() ) {
		auto sp = task_->syncing_progress();
		ui->syncing_progress->setMaximum(sp.second);
		ui->syncing_progress->setValue(sp.first);

		//qDebug() << "TaskView::update_ui_from_task: progress:" << sp;
		ui->syncing->show();

		//ui->syncing_progress->update();
	}
	else {
		ui->syncing->hide();
		//ui->syncing->update();
	}
}


// void TaskView::on_input_set_from_file_clicked()
// {
//     qDebug() << "TaskView::on_input_set_from_file_clicked";

// 	QString file_name = QFileDialog::getOpenFileName(this, tr("Open PDB file"), "", tr("PDB (*.pdb)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
// 	if( not file_name.isEmpty() ) {
// 		// File f(file_name);
// 		// task_->input( std::move(f) );

// 		ui->input_preview->set_pose( core::pose::PoseOP() );
// 		ui->input_preview->update_pose_draw();

// 		//task_->input( File(file_name) );
// 	}
// }


// void TaskView::on_script_set_from_file_clicked()
// {
//     qDebug() << "TaskView::on_script_set_from_file_clicked";

// 	QString file_name = QFileDialog::getOpenFileName(this, tr("Open Rosetta XML script"), "", tr("XML (*.xml)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
// 	if( not file_name.isEmpty() ) {
// 		//task_->script( File(file_name) );
// 	}
// }

// // void TaskView::on_script_textChanged()
// // {
// //     //qDebug() << "TaskView::on_description_textChanged";
// // 	task_->script( ui->script->document()->toPlainText() );
// // }

// void TaskView::on_flags_set_from_file_clicked()
// {
//     qDebug() << "TaskView::on_flags_set_from_file_clicked";

// 	QString file_name = QFileDialog::getOpenFileName(this, tr("Open Rosetta flags script"), "", tr("(*.flags *.*)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
// 	if( not file_name.isEmpty() ) {
// 		//task_->flags( File(file_name) );
// 	}
// }


// void TaskView::on_submit_clicked()
// {
//     qDebug() << "TaskView::on_submit_clicked queue:" << ui->queue->currentText();

// 	Project * project = task_->project();

// 	if( project  and  check_submit_requirements(*project) ) {
// 		ui->submit->setEnabled(false);
// 		task_->submit( ui->queue->currentText() );

// 		save_project(*project, /* always_ask_for_file_name = */ false);
// 	}
// }

void TaskView::create_output_context_menu(const QPoint &pos)
{
	if( /*QStringListModel *model =*/ qobject_cast<QStringListModel*>( ui->output->model() ) ) {
		QPoint item = ui->output->mapToGlobal(pos);

		QMenu submenu;
		submenu.addAction("Open...", this, &TaskView::action_output_open);
		submenu.addAction("Save as...", this, SLOT( action_output_save_as() ));
		/*QAction* rightClickItem =*/ submenu.exec(item);
		// if(rightClickItem && rightClickItem->text().contains("Open...") ) {
		// 	//qDebug() << ui->output->indexAt(pos).row();
		// 	auto index = ui->output->indexAt(pos);
		// 	QVariant qv = model->data(index, Qt::DisplayRole);
		// 	if( qv.isValid() ) {
		// 		QString line = qv.toString();
		// 		std::map<QString, FileSP> const & output = task_->output();
		// 		auto it = output.find(line);
		// 		if( it != output.end() ) {
		// 			qDebug() << "TaskView::create_output_context_menu: open file:" << it->second->file_name();
		// 		}
		// 	}
		// }
	}
}


void TaskView::action_output_open()
{
	if( QStringListModel *model = qobject_cast<QStringListModel*>( ui->output->model() ) ) {
		QModelIndexList indexes = ui->output->selectionModel()->selectedIndexes();

		//QString dir;
		//if( indexes.size() > 1 ) dir = QFileDialog::getExistingDirectory(this, tr("Save Output file to dir...") );

		auto const & files = task_->files();

		for (int i = 0; i <indexes.size(); ++i) {
			QModelIndex index = indexes[i];

			QVariant qv = model->data(index, Qt::DisplayRole);
			if( qv.isValid() ) {
				auto r = index.row();

				if( r < static_cast<int>(files.size()) ) {
					auto it = std::next(files.begin(), r);
					open_file_viewer(*it, task_, this);
				}

				// {
				// 	QString line = qv.toString();
				// 	auto const & files = task_->files();
				// 	auto it = task_->files_find(line);
				// 	if( it != files.end() ) {
				// 		open_file_viewer(*it, task_, this);
				// 	}
				// }
			}
		}
	}
}


void TaskView::action_output_save_as()
{
    //qDebug() << "TaskView::action_output_save_as";

	if( QStringListModel *model = qobject_cast<QStringListModel*>( ui->output->model() ) ) {
		QModelIndexList indexes = ui->output->selectionModel()->selectedIndexes();

		QString dir;
		if( indexes.size() > 1 ) dir = QFileDialog::getExistingDirectory(this, tr("Save Output file to dir...") );

		auto const & files = task_->files();

		for (int i = 0; i <indexes.size(); ++i) {
			QModelIndex index = indexes[i];

			QVariant qv = model->data(index, Qt::DisplayRole);
			if( qv.isValid() ) {
				auto r = index.row();

				if( r < static_cast<int>(files.size()) ) {
					FileSP file_sp = * std::next(files.begin(), r);

					QString file_name;
					if( indexes.size() > 1 ) {
						QFileInfo fi( file_sp->name() );
						file_name = dir + '/' + fi.fileName();
					}
					else {
						QDir qd(".");
						QString suggested_name = qd.absolutePath() + "/" + file_sp->name();
						//qDebug() << "dir: " << suggested_name;
						file_name = QFileDialog::getSaveFileName(this, tr("Save Output file as..."), suggested_name);
					}
					QFile file(file_name);
					if( file.open(QIODevice::WriteOnly) ) {
						qDebug() << "TaskView::action_output_save_as: file:" << file_name;
						file.write( file_sp->data() );
						file.close();
					}
				}

				// {
				// 	QString line = qv.toString();
				// 	Files const & output = task_->files();
				// 	auto it = output.find(line);
				// 	if( it != output.end() ) {
				// 		QString file_name;
				// 		if( indexes.size() > 1 ) {
				// 			QFileInfo fi( it->second->name() );
				// 			file_name = dir + '/' + fi.fileName();
				// 		}
				// 		else {
				// 			QDir qd(".");
				// 			QString suggested_name = qd.absolutePath() + "/" + QFileInfo(it->first).fileName();
				// 			//qDebug() << "dir: " << suggested_name;
				// 			file_name = QFileDialog::getSaveFileName(this, tr("Save Output file as..."), suggested_name);
				// 		}
				// 		QFile file(file_name);
				// 		if( file.open(QIODevice::WriteOnly) ) {
				// 			qDebug() << "TaskView::action_output_save_as: file:" << file_name;
				// 			file.write( it->second->data() );
				// 			file.close();
				// 		}
				// 	}
				// }
			}
		}
	}
}


void TaskView::on_export_all_files_clicked()
{
    qDebug() << "TaskView::on_export_all_files_clicked()...";

	QString dir;
	dir = QFileDialog::getExistingDirectory(this, tr("Save Output file to dir...") );

	for(auto const & it : task_->files() ) {
		auto file_info = QFileInfo( dir + '/' + it->name() );
		auto dir = file_info.dir();
		if( not dir.exists() ) dir.mkpath(".");

		QString file_name = file_info.filePath();  //dir + '/' + it.first;
		QFile file(file_name);
		if( file.open(QIODevice::WriteOnly) ) {
			qDebug() << "TaskView::on_export_all_files_clicked: saving file:" << file_name;
			file.write( it->data() );
			file.close();
		}
	}
}

/*
QWidget * TaskView::create_viewer_for_file(FileSP const &f)
{
	QStringList list = f->file_name().split('.');

	if( not list.isEmpty() ) {
		QString ext = list.last();

		if( false  and  ext == "pdb" ) {
			ui_core::pose_draw::SimplePoseDrawOpenGLWidget *pv = new ui_core::pose_draw::SimplePoseDrawOpenGLWidget(this->ui->preview);

			core::pose::PoseOP pose = std::make_shared<core::pose::Pose>();

			qDebug() << "pose_from_pdbstring...";
			core::import_pose::pose_from_pdbstring(*pose, std::string( f->data().constData(), f->data().length() ), f->file_name().toUtf8().constData() );
			qDebug() << "pose_from_pdbstring... OK";

			pv->set_pose(pose);

			// QProcess * p = new QProcess(this);
			// p->start("MacPyMOL");

			return pv;
		}

		if( ext == "log"  or  ext == "sf"  or true ) { // default viewer
			QPlainTextEdit *e = new QPlainTextEdit(f->data(), this->ui->preview);
			e->setReadOnly(true);

			e->setLineWrapMode(QPlainTextEdit::NoWrap);

			QFont font("monospace");
			e->setFont(font);

			return e;
		}

	}

	return nullptr;
}
*/
void TaskView::on_output_clicked(QModelIndex const &index)
{
	if(auto file = index_to_file(index) ) {
		qDebug() << "TaskView::on_output_clicked: file:" << file->name();
		preview_file(file);
	}
}


void TaskView::on_output_doubleClicked(const QModelIndex &index)
{
	if(auto file = index_to_file(index) ) {
		qDebug() << "TaskView::on_output_doubleClicked: file:" << file->name();
		open_file_viewer(file, task_, this);
	}
}


void TaskView::file_changed(FileID const &file_id)
{
	if(file_id == previewed_file_) {
		Files const & files = task_->files();

		auto it = task_->files_find(file_id);
		if( it != files.end() ) preview_file(*it);
	}
}


void TaskView::preview_file(FileSP const & file)
{
	previewed_file_ = file->file_id();

	if(viewer_) viewer_->deleteLater();
	viewer_ = create_viewer_for_file(ui->preview, file);

	if(viewer_) {
		viewer_->resize( this->ui->preview->size() );
		viewer_->show();
	}
}


/// return file object that correspond to given model index, return null if index could not be mapped to file
FileSP TaskView::index_to_file(QModelIndex const &index)
{
	if( QStringListModel *model = qobject_cast<QStringListModel*>( ui->output->model() ) ) {
		QVariant qv = model->data(index, Qt::DisplayRole);
		if( qv.isValid() ) {
			auto r = index.row();

			auto const & files = task_->files();

			if( r < static_cast<int>(files.size()) ) {
				FileSP file = * std::next(files.begin(), r);
				return file;
			}
		}
	}
	return FileSP();
}


void TaskView::on_cancel_task_clicked()
{
	TaskCancelDialog tcd;
	tcd.run(task_);
}


void open_file_viewer(FileSP const &file, TaskSP const &task, QWidget *parent)
{
	//QFileInfo fi( name_and_file.first );
	//QString file_name = fi.fileName();
	QString name = file->name();

	if( file->data().isEmpty() ) return;

	//qDebug() << "temp: " << it->first << file_name;

	if( name.endsWith(".pdb") ) {
		QTemporaryFile * tmp_file = new QTemporaryFile("XXXXXX." + name, parent);
		if( tmp_file->open() ) {
			QString file_system_name = tmp_file->fileName();
			qDebug() << "file_system_name:" << file_system_name;
			tmp_file->write( file->data() );
			tmp_file->close();

			QString program = config::get_pdb_viewer_path();

			if( program.isEmpty() ) {

				QMessageBox msgBox;
				msgBox.setText( QString("External PDB viewer app is not set, would you like to set it now?") );
				msgBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
				msgBox.setDefaultButton(QMessageBox::Cancel);

				if( msgBox.exec() == QMessageBox::Ok) {
					config::ConfigDialog * preferences = new config::ConfigDialog(parent);
					preferences->show();
				}
			}
			else {
				QStringList arguments;
				arguments << file_system_name; //arguments << "-style" << "fusion";
				QProcess * process = new QProcess(parent);
				process->start(program, arguments);

				//QProcess *myProcess = new QProcess(this);
				//myProcess->start(program, arguments);
				//process->setProgram(program);
				//process->setArguments(arguments);
				//process->startDetached()
			}
		}
	} else {
		for(auto const & ending : {"score", ".sc", ".out", ".silent"} ) {
			if( name.endsWith(ending) ) {
				auto sfv = new viewers::ScoreFileView(file, task, nullptr);
				sfv->show();
				break;
			}
		}
		for(auto const & contain : {".out."} ) {
			if( name.contains(contain) ) {
				auto sfv = new viewers::ScoreFileView(file, task, nullptr);
				sfv->show();
				break;
			}
		}

		// if( file_name.endsWith(".score") or file_name.endsWith(".sc")  or  file_name.endsWith(".out")  or  file_name.endsWith(".silent") ) {
		// auto sfv = new viewers::ScoreFileView(name_and_file, task, nullptr);
		// sfv->show();
	}
}


} // namespace task
} // namespace ui
