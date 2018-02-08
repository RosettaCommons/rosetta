#ifndef JOB_VIEW_H
#define JOB_VIEW_H

#include <ui/task/task.h>

#include <QWidget>

namespace Ui {
class JobView;
}

namespace ui {
namespace task {

class JobView : public QWidget
{
    Q_OBJECT

public:
    explicit JobView(Task::JobSP const &, TaskWP const &task, QWidget *parent = nullptr);
    ~JobView();

private:
	Task::JobSP job_;
	TaskWP task_;

    Ui::JobView *ui;
};

} // namespace task
} // namespace ui

#endif // JOB_VIEW_H
