#include "task_cancel_dialog.h"
#include "ui_task_cancel_dialog.h"


#include <ui/task/task.h>

#include <QMessageBox>
#include <QJsonObject>

namespace ui {
namespace task {

TaskCancelDialog::TaskCancelDialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::TaskCancelDialog)
{
    ui->setupUi(this);
}

TaskCancelDialog::~TaskCancelDialog()
{
    delete ui;
}


int TaskCancelDialog::run(TaskSP const &task)
{
	if( task->task_id().isEmpty() ) {
		QMessageBox box(QMessageBox::Critical, "Error cancaleing Task!", QString("Can not cancel Task because it has not being submitted yet!") );
		box.exec();
		return QDialog::Rejected;
	}
	else if( task->state() == Task::State::finished ) {
		QMessageBox box(QMessageBox::Critical, "Error cancaleing Task!", QString("Can not cancel Task because it already finished!") );
		box.exec();
		return QDialog::Rejected;
	}
	else if( task->state() != Task::State::canceled ) {
		QMessageBox box(QMessageBox::Critical, QString("Cancal task№%1").arg(task->task_id()), QString("This will abort all jobs of Task№%1.\nAre you sure you want to proceed?").arg(task->task_id()), QMessageBox::Yes | QMessageBox::No);
		box.setDefaultButton(QMessageBox::No);
		if( box.exec() == QMessageBox::Yes ) {
			task_ = task;

			auto state = task->state_;

			task->state_ = Task::State::canceled;
			QJsonValue jvalue = task->task_data();
			task->state_ = state;

			QJsonObject jobject = jvalue.toObject();
			QJsonDocument jd = QJsonDocument(jobject);
			network_call_.call(task_api_url() + "/task/" + task->task_id(), QNetworkAccessManager::Operation::PostOperation, jd );

			connect(&network_call_, &NetworkCall::finished, this, &TaskCancelDialog::network_call_finished);

			//network_call_.call(task_api_url() + "/task/" + task->task_id() + "/output/" + file_name, QNetworkAccessManager::Operation::GetOperation);
			//task['state'] = S_canceled

			return exec(); //QDialog::Rejected;
		}
	}

	return QDialog::Rejected;
}

void TaskCancelDialog::network_call_finished()
{
	if(task_) {
		task_->state_ = Task::State::canceled;
		Q_EMIT task_->state_changed();
		Q_EMIT task_->changed();

		Q_EMIT done(QDialog::Accepted);
	}
}


} // namespace task
} // namespace ui
