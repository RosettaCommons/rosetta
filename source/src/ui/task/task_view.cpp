#include "task_view.h"
#include "ui_task_view.h"

#include <ui/task/task.h>

#include <ui/ui_core/pose_draw/SimplePoseDrawOpenGLWidget.h>

#include <QFileDialog>
#include <QStringListModel>
#include <QProcess>
#include <QMenu>
#include <QJsonObject>
#include <QJsonArray>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>

namespace ui {
namespace task {

TaskView::TaskView(TaskSP const &task,QWidget *parent) :
	QWidget(parent),
	ui(new Ui::TaskView),
	task_(task)
{
    ui->setupUi(this);

	auto model = new QStringListModel(this);
	ui->output->setModel(model);

	//bool old = task->blockSignals(true); // we now handling this on Model (ie Task layer) by ommiting change() emiition if data is the same
	update_ui_from_task();
	//task->blockSignals(old);

	connect(task_.get(), SIGNAL(changed()), this, SLOT(update_ui_from_task()));

	ui->output->setContextMenuPolicy(Qt::CustomContextMenu);
	connect(ui->output, SIGNAL(customContextMenuRequested(const QPoint &)), this, SLOT(create_output_context_menu(const QPoint &)));

	setAttribute(Qt::WA_DeleteOnClose);

	connect(&queues, &NetworkCall::finished, this, &TaskView::on_queues_finished);
	queues.call("queues");
}

TaskView::~TaskView()
{
    //qDebug() << "TaskView::~TaskView";
    delete ui;
}

void TaskView::update_ui_from_task()
{
    //qDebug() << "TaskView::update_ui_from_task";

	ui->name->setText( task_->project() ? *task_->project()->find(task_.get()) : "no-project");

	ui->state->setText( Task::to_string( task_->state() ) );

	ui->description->document()->setPlainText( task_->description() );

	ui->input_file_name->setText( task_->input().file_name() );
	//ui->input->document()->setPlainText( task_->input().data() );

	ui->script_file_name->setText( task_->script().file_name() );
	ui->script->document()->setPlainText( task_->script().data() );

	ui->flags_file_name->setText( task_->flags().file_name() );
	ui->flags->document()->setPlainText( task_->flags().data() );

	if( auto model = qobject_cast<QStringListModel*>( ui->output->model() ) ) {
		QStringList output_files;
		for(auto const &k_fsp : task_->output() ) {
			//qDebug() << "TaskView::update_ui_from_task adding output file: " << k_fsp.first;
			if( not k_fsp.second->data().isEmpty() ) output_files << k_fsp.first;
		}
		model->setStringList(output_files);
	} else {
		qDebug() << "TaskView::update_ui_from_task: no QStringListModel!!!";
	}

	if( task_->input().data().size() and (not ui->input_preview->pose()) ) {
		core::pose::PoseOP pose = std::make_shared<core::pose::Pose>();

		core::import_pose::pose_from_pdbstring(*pose, std::string( task_->input().data().constData(), task_->input().data().length() ), task_->input().file_name().toUtf8().constData() );

		ui->input_preview->set_pose(pose);
	}

	auto syncing = task_->is_syncing();
	if( task_->state() == Task::State::_draft_  and (not syncing) ) ui->submit->show();
	else ui->submit->hide();

	if(syncing) {
		auto sp = task_->syncing_progress();
		ui->syncing_progress->setMaximum(sp.second);
		ui->syncing_progress->setValue(sp.second - sp.first);

		qDebug() << "TaskView::update_ui_from_task: progress:" << sp;
		ui->syncing->show();
	}
	else ui->syncing->hide();
}

void TaskView::on_description_textChanged()
{
    //qDebug() << "TaskView::on_description_textChanged";
	task_->description( ui->description->document()->toPlainText() );
}


void TaskView::on_input_set_from_file_clicked()
{
    qDebug() << "TaskView::on_input_set_from_file_clicked";

	QString file_name = QFileDialog::getOpenFileName(this, tr("Open PDB file"), "", tr("PDB (*.pdb)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
	if( not file_name.isEmpty() ) {
		// File f(file_name);
		// task_->input( std::move(f) );

		ui->input_preview->set_pose( core::pose::PoseOP() );
		ui->input_preview->update_pose_draw();

		task_->input( File(file_name) );
	}
}


void TaskView::on_script_set_from_file_clicked()
{
    qDebug() << "TaskView::on_script_set_from_file_clicked";

	QString file_name = QFileDialog::getOpenFileName(this, tr("Open Rosetta XML script"), "", tr("XML (*.xml)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
	if( not file_name.isEmpty() ) {
		task_->script( File(file_name) );
	}
}

// void TaskView::on_script_textChanged()
// {
//     //qDebug() << "TaskView::on_description_textChanged";
// 	task_->script( ui->script->document()->toPlainText() );
// }

void TaskView::on_flags_set_from_file_clicked()
{
    qDebug() << "TaskView::on_flags_set_from_file_clicked";

	QString file_name = QFileDialog::getOpenFileName(this, tr("Open Rosetta flags script"), "", tr("(*.flags *.*)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
	if( not file_name.isEmpty() ) {
		task_->flags( File(file_name) );
	}
}


void TaskView::on_submit_clicked()
{
    qDebug() << "TaskView::on_submit_clicked queue:" << ui->queue->currentText();
	ui->submit->setEnabled(false);
	task_->submit( ui->queue->currentText() );
}

void TaskView::create_output_context_menu(const QPoint &pos)
{
	if( /*QStringListModel *model =*/ qobject_cast<QStringListModel*>( ui->output->model() ) ) {
		QPoint item = ui->output->mapToGlobal(pos);

		QMenu submenu;
		submenu.addAction("Open...");
		submenu.addAction("Save as...", this, SLOT( on_action_output_save_as() ));
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


void TaskView::on_action_output_save_as()
{
    //qDebug() << "TaskView::on_action_output_save_as";

	if( QStringListModel *model = qobject_cast<QStringListModel*>( ui->output->model() ) ) {
		QModelIndexList indexes = ui->output->selectionModel()->selectedIndexes();

		QString dir;
		if( indexes.size() > 1 ) dir = QFileDialog::getExistingDirectory(this, tr("Save Output file to dir...") );

		for (int i = 0; i <indexes.size(); ++i) {
			QModelIndex index = indexes[i];

			QVariant qv = model->data(index, Qt::DisplayRole);
			if( qv.isValid() ) {
				QString line = qv.toString();

				std::map<QString, FileSP> const & output = task_->output();

				auto it = output.find(line);
				if( it != output.end() ) {

					QString file_name;

					if( indexes.size() > 1 ) {
						QFileInfo fi( it->second->file_name() );
						file_name = dir + '/' + fi.fileName();
					}
					else file_name = QFileDialog::getSaveFileName(this, tr("Save Output file as..."), it->second->file_name());

					QFile file(file_name);
					if( file.open(QIODevice::WriteOnly) ) {
						qDebug() << "TaskView::on_action_output_save_as: file:" << file_name;
						file.write( it->second->data() );
						file.close();
					}
				}
			}
		}
	}
}


QWidget * TaskView::create_viewer_for_file(FileSP const &f)
{
	QStringList list = f->file_name().split('.');

	if( not list.isEmpty() ) {
		QString ext = list.last();

		if( ext == "pdb" ) {
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

void TaskView::on_output_clicked(const QModelIndex &index)
{
	if( QStringListModel *model = qobject_cast<QStringListModel*>( ui->output->model() ) ) {
		QVariant qv = model->data(index, Qt::DisplayRole);
		if( qv.isValid() ) {
			QString line = qv.toString();

			std::map<QString, FileSP> const & output = task_->output();

			auto it = output.find(line);
			if( it != output.end() ) {
				qDebug() << "TaskView::on_output_clicked: file:" << it->second->file_name();

				if(viewer_) viewer_->deleteLater();
				viewer_ = create_viewer_for_file(it->second);

				if(viewer_) {
					viewer_->resize( this->ui->preview->size() );
					viewer_->show();
				}
			}
		}
	}
}


void TaskView::on_queues_finished()
{
	qDebug() << "TaskView::on_queues_finished data:" << queues.result();
	QJsonDocument jd = queues.json();
	QJsonArray root = jd.array();

    for (QJsonArray::const_iterator it = root.constBegin(); it != root.constEnd(); ++it) {
		//qDebug() << "TaskView::on_queues_finished: " << *it;
		ui->queue->addItem( it->toString() );
    }

	if( not root.isEmpty() ) ui->submit->setEnabled(true);
}


} // namespace task
} // namespace ui
