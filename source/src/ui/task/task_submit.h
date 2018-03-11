#ifndef TASK_SUBMIT_H
#define TASK_SUBMIT_H

#include <ui/task/task.fwd.h>
#include <ui/task/file.fwd.h>
#include <ui/task/util.h>

#include <QWidget>

namespace Ui {
class TaskSubmit;
}

namespace ui {
namespace task {

QWidget * create_viewer_for_file(QWidget *parent, FileSP const &);

class FileTableWatcher : public QObject
{
	Q_OBJECT

	bool eventFilter(QObject * obj, QEvent * event) override;

	Q_SIGNAL void backspace_pressed();

public:
	using QObject::QObject;
};

class TaskSubmit : public QWidget
{
    Q_OBJECT

public:
    explicit TaskSubmit(TaskSP const &, QWidget *parent = 0);
    ~TaskSubmit();


private Q_SLOTS:
	void update_ui_from_task();

	void on_name_textChanged(QString const &);
	//void on_name_textEdited(QString const &text);


	void on_version_textChanged(QString const &);

	void on_description_textChanged();

	void on_add_files_clicked();

	void on_add_input_structure_clicked();
	void on_add_native_structure_clicked();

	void on_add_job_clicked();
	void s_on_jobs_tabMoved(int from, int to);

	void backspace_pressed_on_files();


	void on_files_clicked(const QModelIndex &index);

	void s_on_queues_finished();

	void on_submit_clicked();


	//void on_app_activated(const QString &text);
	//void on_nstruct_valueChanged(int);
	//void on_flags_textChanged();
	//void on_flags_from_file_clicked();


private:
	Ui::TaskSubmit *ui;

	TaskSP task_;

	NetworkCall queues;

	QWidget *viewer_ = nullptr;
};


} // namespace task
} // namespace ui

#endif // TASK_SUBMIT_H
