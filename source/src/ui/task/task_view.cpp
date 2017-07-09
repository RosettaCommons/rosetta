#include "task_view.h"
#include "ui_task_view.h"

#include <ui/task/task.h>

#include <QFileDialog>

namespace ui {
namespace task {

TaskView::TaskView(TaskSP const &task,QWidget *parent) :
	QWidget(parent),
	ui(new Ui::TaskView),
	task_(task)
{
    ui->setupUi(this);

	//bool old = task->blockSignals(true); // we now handling this on Model (ie Task layer) by ommiting change() emiition if data is the same
	update_ui_from_task();
	//task->blockSignals(old);

	connect(task_.get(), SIGNAL(changed()), this, SLOT(update_ui_from_task()));

	setAttribute(Qt::WA_DeleteOnClose);
}

TaskView::~TaskView()
{
    qDebug() << "TaskView::~TaskView";
    delete ui;
}

void TaskView::update_ui_from_task()
{
    //qDebug() << "TaskView::update_ui_from_task";

	ui->name->setText( task_->project() ? *task_->project()->find(task_.get()) : "no-project");

	ui->state->setText( task_->state() );

	ui->description->document()->setPlainText( task_->description() );

	ui->input_file_name->setText( task_->input().file_name() );
	//ui->input->document()->setPlainText( task_->input().data() );

	ui->script_file_name->setText( task_->script().file_name() );
	ui->script->document()->setPlainText( task_->script().data() );

	ui->flags_file_name->setText( task_->flags().file_name() );
	ui->flags->document()->setPlainText( task_->flags().data() );
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
    qDebug() << "TaskView::on_submit_clicked";
	task_->submit();
}



} // namespace task
} // namespace ui
