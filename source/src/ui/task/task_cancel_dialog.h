#ifndef TASK_CANCEL_DIALOG_H
#define TASK_CANCEL_DIALOG_H

#include <ui/task/task.fwd.h>
#include <ui/task/util.h>


#include <QDialog>

namespace Ui {
class TaskCancelDialog;
}

namespace ui {
namespace task {


class TaskCancelDialog : public QDialog
{
    Q_OBJECT

public:
    explicit TaskCancelDialog(QWidget *parent = nullptr);
    ~TaskCancelDialog();

	int run(TaskSP const &);

private Q_SLOTS:
	void network_call_finished();

private:
	TaskSP task_;

	NetworkCall network_call_;

    Ui::TaskCancelDialog *ui;
};

} // namespace task
} // namespace ui

#endif // TASK_CANCEL_DIALOG_H
