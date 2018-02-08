#include <ui/task/job_view.h>
#include "ui_job_view.h"

namespace ui {
namespace task {

JobView::JobView(Task::JobSP const &job, TaskWP const &task, QWidget *parent) :
    QWidget(parent),
	job_(job),
	task_(task),
    ui(new Ui::JobView)
{
    ui->setupUi(this);

	ui->app->setText( job_->app );
	ui->nstruct->setText( job_->nstruct );
	ui->flags->document()->setPlainText( job_->flags );

	setAttribute(Qt::WA_DeleteOnClose);
}

JobView::~JobView()
{
    delete ui;
}

} // namespace task
} // namespace ui
