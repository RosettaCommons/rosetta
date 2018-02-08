#ifndef JOB_SUBMIT_H
#define JOB_SUBMIT_H

#include <ui/task/task.h>

#include <QWidget>

namespace Ui {
class JobSubmit;
}

namespace ui {
namespace task {

class JobSubmit : public QWidget
{
    Q_OBJECT

public:
    explicit JobSubmit(Task::JobSP const &, TaskWP const &task, QWidget *parent = nullptr);
    ~JobSubmit();

private:
	void update_ui_from_job();

private Q_SLOTS:
	//void on_name_editingFinished();
	void on_nstruct_textChanged();

	void on_app_activated(const QString &text);

	void on_flags_textChanged();

	void on_flags_from_file_clicked();

	void on_rename_job_clicked();
	void on_delete_job_clicked();



private:
	Task::JobSP job_;
	TaskWP task_;

    Ui::JobSubmit *ui;
};

} // namespace task
} // namespace ui

#endif // JOB_SUBMIT_H
