#ifndef TASK_VIEW_H
#define TASK_VIEW_H

#include <ui/task/task.fwd.h>
#include <ui/task/file.fwd.h>

#include <ui/task/util.h>

#include <QWidget>

namespace Ui {
class TaskView;
}

namespace ui {
namespace task {

void open_file_viewer( std::pair<QString const, FileSP> const & file, TaskSP const &, QWidget *parent=nullptr);


class TaskView : public QWidget
{
    Q_OBJECT
public:
    explicit TaskView(TaskSP const &, QWidget *parent = 0);
	~TaskView();

Q_SIGNALS:

public Q_SLOTS:

private Q_SLOTS:
	// void on_description_textChanged(void);
	// void on_nstruct_valueChanged(int);

	// void on_input_set_from_file_clicked();
	// void on_flags_set_from_file_clicked();
	// void on_script_set_from_file_clicked();

	void update_ui_from_task();
	void update_ui_file_list_from_task();

	//void on_submit_clicked();

    void on_output_clicked(const QModelIndex &index);
	void on_output_doubleClicked(const QModelIndex &index);

	void create_output_context_menu(const QPoint &pos);

	void action_output_open();
	void action_output_save_as();

	void on_export_all_files_clicked();

	void update_syncing_progress();

	void file_changed(QString const &file_name);

private:
	void preview_file(std::pair<QString const, FileSP> const & name_and_file);

private:
	//QWidget * create_viewer_for_file(FileSP const &);

	Ui::TaskView *ui;

	TaskSP task_;

	QWidget *viewer_ = nullptr;

	QString previewed_file_name_;
};

} // namespace task
} // namespace ui

#endif // TASK_VIEW_H
