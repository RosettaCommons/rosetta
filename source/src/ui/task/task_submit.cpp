#include "task_submit.h"
#include "ui_task_submit.h"

#include <ui/task/task.h>
#include <ui/task/task_view.h>
#include <ui/task/job_submit.h>

#include <ui/util/apps.gen.h>

#include <QFileDialog>
#include <QKeyEvent>
#include <QDebug>
#include <QJsonArray>

#include <set>

namespace ui {
namespace task {

bool FileTableWatcher::eventFilter(QObject * obj, QEvent * event)
{
	//auto table = qobject_cast<QTableWidget*>(receiver);
	//if (table  and  event->type() == QEvent::KeyPress) {
	if (event->type() == QEvent::KeyPress) {
		auto keyEvent = static_cast<QKeyEvent*>(event);
		if( keyEvent->key() == Qt::Key_Backspace ) {
			Q_EMIT backspace_pressed(); // table->currentRow(), table->currentColumn());
			return true;
		}
	}
	// standard event processing
	return QObject::eventFilter(obj, event);
}


TaskSubmit::TaskSubmit(TaskSP const &task, QWidget *parent) :
    QWidget(parent),
    ui(new Ui::TaskSubmit),
	task_(task)
{
    ui->setupUi(this);

	ui->files->setModel( task_->files_model() );
	ui->files->horizontalHeader()->setStretchLastSection(true);
	auto fw = new FileTableWatcher(this);
	ui->files->installEventFilter(fw);
	connect(fw, SIGNAL(backspace_pressed()), this, SLOT(backspace_pressed_on_files()));

	//auto apps = {"AbinitioRelax", "docking_protocol", "relax"};
	//for(auto const & a : util::rosetta_app_list) ui->app->addItem(a);

	update_ui_from_task();
	connect(task_.get(), SIGNAL(changed()), this, SLOT(update_ui_from_task()));

	connect(ui->jobs->tabBar(), &QTabBar::tabMoved, this, &TaskSubmit::s_on_jobs_tabMoved);

	connect(&queues, &NetworkCall::finished, this, &TaskSubmit::s_on_queues_finished);
	queues.call(task_api_url() + '/' + "queues");

	setAttribute(Qt::WA_DeleteOnClose);
}

TaskSubmit::~TaskSubmit()
{
	//qDebug() << "TaskSubmit::~TaskSubmit()";
    delete ui;
}


void TaskSubmit::update_ui_from_task()
{
    //qDebug() << "TaskSubmit::update_ui_from_task";

	ui->name->setText( task_->name() );

	//ui->app->setCurrentText( task_->app() );

	ui->version->setText( task_->version() );

	//ui->state->setText( Task::to_string( task_->state() ) );

	ui->description->document()->setPlainText( task_->description() );


	//ui->flags->document()->setPlainText( task_->flags() );


	// ui->input_file_name->setText( task_->input().file_name() );
	//ui->input->document()->setPlainText( task_->input().data() );

	// ui->script_file_name->setText( task_->script().file_name() );
	// ui->script->document()->setPlainText( task_->script().data() );

	// ui->flags_file_name->setText( task_->flags().file_name() );
	// ui->flags->document()->setPlainText( task_->flags().data() );

	//ui->nstruct->setValue( task_->nstruct() );

	for(int i=ui->jobs->count(); i>0 ; --i) ui->jobs->removeTab(i-1);

	for(auto const & job : task_->jobs() ) {
		auto js = new JobSubmit(job, task_);
		ui->jobs->addTab(js, job->name);
	}
}

void TaskSubmit::on_name_textChanged(QString const &text)
{
    //qDebug() << "TaskSubmit::on_name_textChanged";
	task_->name( text /*ui->name->text() */);
}

// void TaskSubmit::on_app_activated(const QString &text)
// {
//     qDebug() << "TaskSubmit::on_app_currentIndexChanged";
// 	task_->app( text /*ui->app->currentText()*/ );
// }

void TaskSubmit::on_version_textChanged(QString const &text)
{
    //qDebug() << "TaskSubmit::on_name_textChanged";
	task_->version( text /*ui->version->text()*/ );
}

// void TaskSubmit::on_nstruct_valueChanged(int)
// {
//     //qDebug() << "TaskSubmit::on_nstruct_valueChanged";
// 	task_->nstruct( ui->nstruct->value() );
// }

void TaskSubmit::on_description_textChanged()
{
    //qDebug() << "TaskSubmit::on_description_textChanged";
	task_->description( ui->description->document()->toPlainText() );
}

// void TaskSubmit::on_flags_textChanged()
// {
// 	task_->flags( ui->flags->document()->toPlainText() );
// }

void TaskSubmit::on_add_input_structure_clicked()
{
    //qDebug() << "TaskSubmit::add_input_structure";

	QString file_name = QFileDialog::getOpenFileName(this, tr("Open PDB file"), "", tr("PDB (*.pdb)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
	if( not file_name.isEmpty() ) {
		task_->add_file( "input.pdb", std::make_shared<File>(file_name) );
	}
}

void TaskSubmit::on_add_native_structure_clicked()
{
    //qDebug() << "TaskSubmit::on_add_native_structure_clicked";

	QString file_name = QFileDialog::getOpenFileName(this, tr("Open PDB file"), "", tr("PDB (*.pdb)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
	if( not file_name.isEmpty() ) {
		task_->add_file( "native.pdb", std::make_shared<File>(file_name) );
	}
}

// void TaskSubmit::on_flags_from_file_clicked()
// {
//     //qDebug() << "TaskSubmit::on_flags_from_file_clicked";
// 	QString file_name = QFileDialog::getOpenFileName(this, tr("Open Rosetta flags script"), "", "", Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
// 	if( not file_name.isEmpty() ) {
// 		QFile file(file_name);
// 		if( file.open(QIODevice::ReadOnly) ) {
// 			task_->flags( file.readAll() );
// 		}
// 	}
// }

void TaskSubmit::backspace_pressed_on_files()
{
	qDebug() << "TaskSubmit::backspace_pressed_on_files...";

	QItemSelectionModel *select = ui->files->selectionModel();

	if( select->hasSelection() ) {
		auto rows = select->selectedIndexes();

		std::set<int> selection;
		for(auto const & index : rows) selection.insert( index.row() );

		for(auto r : selection) qDebug() << "row:" << r;
	}
}

void TaskSubmit::on_add_files_clicked()
{
	qDebug() << "TaskSubmit::on_add_files_clicked...";

	QStringList const file_names = QFileDialog::getOpenFileNames(this, tr("Add Files"), "", tr(""), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
	if( not file_names.isEmpty() ) {
		for(auto it = file_names.constBegin(); it != file_names.constEnd(); ++ it) {
			qDebug() << "file: " << *it;
			QString name = QFileInfo(*it).fileName();

			// if( name == "flags" ) {
			// 	QFile file(*it);
			// 	if( file.open(QIODevice::ReadOnly) ) {
			// 		task_->flags( file.readAll() );
			// 	}
			// }
			// else task_->add_file(name, std::make_shared<File>(*it) );

			task_->add_file(name, std::make_shared<File>(*it) );

			// // File f(file_name);
			// // task_->input( std::move(f) );

			// ui->input_preview->set_pose( core::pose::PoseOP() );
			// ui->input_preview->update_pose_draw();

			// task_->input( File(file_name) );
		}
	}
}

void TaskSubmit::on_add_job_clicked()
{
	//qDebug() << "TaskSubmit::on_add_job_clicked...";

	auto job = task_->add_job();
	job->app = *util::rosetta_app_list.begin();

	auto js = new JobSubmit(job, task_);

	ui->jobs->setCurrentIndex( ui->jobs->addTab(js, job->name) );
}

void TaskSubmit::s_on_jobs_tabMoved(int from, int to)
{
	qDebug() << "TaskSubmit::on_jobs_tabMoved:" << from << " ->" << to;

	task_->swap_jobs(from, to);
}


void TaskSubmit::on_files_clicked(const QModelIndex &index)
{
	if( FileTableModel *model = qobject_cast<FileTableModel*>( ui->files->model() ) ) {
		// QVariant qv = model->data(index, Qt::DisplayRole);
		// if( qv.isValid() ) {
		// 	QString line = qv.toString();
		// 	std::map<QString, FileSP> const & files = task_->files();
		// 	auto it = files.find(line);
		// 	if( it != files.end() ) {
		// 		qDebug() << "TaskView::on_output_clicked: file:" << it->second->file_name();
		// 		if(viewer_) viewer_->deleteLater();
		// 		viewer_ = create_viewer_for_file(ui->preview, it->second);
		// 		if(viewer_) {
		// 			viewer_->resize( this->ui->preview->size() );
		// 			viewer_->show();
		// 		}
		// 	}
		// }
		std::map<QString, FileSP> const & files = task_->files();
		if( index.isValid()  and  index.row() >=0  and  index.row() < static_cast<int>( files.size() ) ) {
			auto it = files.begin();
			std::advance(it, index.row() );

			if(viewer_) viewer_->deleteLater();
			viewer_ = create_viewer_for_file(ui->preview, it->second);
			if(viewer_) {
				viewer_->resize( this->ui->preview->size() );
				viewer_->show();
			}

		}
	}
}

	void TaskSubmit::s_on_queues_finished()
{
	qDebug() << "TaskSubmit::on_queues_finished data:" << queues.result();
	QJsonDocument jd = queues.result_as_json();
	QJsonArray root = jd.array();

    for (QJsonArray::const_iterator it = root.constBegin(); it != root.constEnd(); ++it) {
		//qDebug() << "TaskView::on_queues_finished: " << *it;
		ui->queue->addItem( it->toString() );
    }

	if( not root.isEmpty() ) ui->submit->setEnabled(true);
}

void TaskSubmit::on_submit_clicked()
{
    qDebug() << "TaskSubmit::on_submit_clicked queue:" << ui->queue->currentText();

	//Project * project = task_->project();

	ui->submit->setEnabled(false);

	auto const & files = task_->files();

	std::vector<QString> names;
	for(auto const & it : files) names.push_back(it.first);
	for(auto const & n : names) task_->rename_file(n, "input/" + n);

	task_->submit( ui->queue->currentText() );

	// if( project  and  check_submit_requirements(*project) ) {
	// 	ui->submit->setEnabled(false);
	// 	task_->submit( ui->queue->currentText() );
	// 	save_project(*project, /* always_ask_for_file_name = */ false);
	// }

	auto tv = new TaskView(task_);
	tv->show();

	this->close();
}


} // namespace task
} // namespace ui
