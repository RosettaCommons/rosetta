#include "job_submit.h"
#include "ui_job_submit.h"

#include <ui/util/apps.gen.h>

#include <QIntValidator>
#include <QInputDialog>
#include <QMessageBox>
#include <QFileDialog>

namespace ui {
namespace task {

JobSubmit::JobSubmit(Task::JobSP const &job, TaskWP const &task, QWidget *parent) :
    QWidget(parent),
	job_(job),
	task_(task),
    ui(new Ui::JobSubmit)
{
    ui->setupUi(this);

	for(auto const & a : util::rosetta_app_list) ui->app->addItem(a);

	ui->nstruct->setValidator( new QIntValidator(1, 10000, this) );

	update_ui_from_job();

	setAttribute(Qt::WA_DeleteOnClose);
}


JobSubmit::~JobSubmit()

{
	//qDebug() << "JobSubmit::~JobSubmit():" << job_->name;
    delete ui;
}


void JobSubmit::update_ui_from_job()
{
	//ui->name->setText( job_->name );
	ui->nstruct->setText( job_->nstruct );

	//ui->nstruct->setText( QString::number(job_->nstruct) );

	//qDebug() << "setting app to:" << job_->app;
	ui->app->setCurrentText( job_->app );

	ui->flags->document()->setPlainText( job_->flags );
}


// void JobSubmit::on_name_editingFinished()
// {
// 	qDebug() << "JobSubmit::on_name_editingFinished()...";
// 	//job_->name = ui->name->text();
// 	if( auto task = task_.lock() ) {
// 		if( task->rename_job(job_, ui->name->text() ) ) {
// 			// qDebug() << "JobSubmit::parent:" << parentWidget() << " parent-of-parent:" << parentWidget()->parentWidget();
// 			// JobSubmit::parent: QStackedWidget(0x7f987d615590, name="qt_tabwidget_stackedwidget")  parent-of-parent: QTabWidget(0x7f987d615310, name="jobs")
// 			if( QTabWidget *tw = qobject_cast<QTabWidget *>( parentWidget()->parentWidget() ) ) {
// 				//QTabBar * tb = tw->tabBar();
// 				//tb->setTabText( tb->currentIndex(), job_->name );
// 				tw->setTabText( tw->currentIndex(), job_->name );
// 			}
// 		}
// 	}
// 	ui->name->setText( job_->name );
// }

void JobSubmit::on_nstruct_textChanged()
{
	job_->nstruct = ui->nstruct->text();
}



void JobSubmit::on_app_activated(QString const &text)
{
	//qDebug() << "JobSubmit::on_app_activated():" << text;
	job_->app = text;
}

void JobSubmit::on_flags_textChanged()
{
	job_->flags = ui->flags->document()->toPlainText();
}

void JobSubmit::on_flags_from_file_clicked()
{
    //qDebug() << "JobSubmit::on_flags_from_file_clicked";
	QString file_name = QFileDialog::getOpenFileName(this, tr("Open Rosetta flags script"), "", "", Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
	if( not file_name.isEmpty() ) {
		QFile file(file_name);
		if( file.open(QIODevice::ReadOnly) ) {
			job_->flags = file.readAll();
			ui->flags->document()->setPlainText( job_->flags );
		}
	}
}


void JobSubmit::on_rename_job_clicked()
{
	//qDebug() << "JobSubmit::on_rename_job_clicked()";
	if( auto task = task_.lock() ) {
		bool ok;
		QString text = QInputDialog::getText(this, tr("QInputDialog::getText()"),
											 tr("New job name (must be unique):"), QLineEdit::Normal,
											 job_->name, &ok);

		if( ok  and  !text.isEmpty()  and  task->rename_job(job_, text ) ) {
			qDebug() << "JobSubmit::on_rename_job_clicked(): renaming to:" << text;
			//qDebug() << "JobSubmit::parent:" << parentWidget() << " parent-of-parent:" << parentWidget()->parentWidget();
			// JobSubmit::parent: QStackedWidget(0x7f987d615590, name="qt_tabwidget_stackedwidget")  parent-of-parent: QTabWidget(0x7f987d615310, name="jobs")
			if( QTabWidget *tw = qobject_cast<QTabWidget *>( parentWidget()->parentWidget() ) ) {
				//QTabBar * tb = tw->tabBar();
				//tb->setTabText( tb->currentIndex(), job_->name );
				tw->setTabText( tw->currentIndex(), job_->name );
			}
		}
	}
}

void JobSubmit::on_delete_job_clicked()
{
	if( auto task = task_.lock() ) {
		QMessageBox msgBox;
		//msgBox.setText("The document has been modified.");
		//msgBox.setInformativeText( QString("Do you want to delete job %1?").arg(job_->name) );
		msgBox.setText( QString("Do you want to delete job %1?").arg(job_->name) );
		msgBox.setStandardButtons(QMessageBox::Ok | QMessageBox::Cancel);
		msgBox.setDefaultButton(QMessageBox::Cancel);

		if( msgBox.exec() == QMessageBox::Ok) {
			qDebug() << "JobSubmit::on_delete_job_clicked(): deleting job:" << job_->name;
			task->delete_job(job_);
			if( QTabWidget *tw = qobject_cast<QTabWidget *>( parentWidget()->parentWidget() ) ) {
				tw->removeTab( tw->currentIndex() ); // technically we do not need this since `close()` below will remove tab
				close();
			}
		}
	}
	qDebug() << "JobSubmit::on_delete_job_clicked(): Done!";
}

} // namespace task
} // namespace ui
