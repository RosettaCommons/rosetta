#ifndef PROJECT_VIEW_H
#define PROJECT_VIEW_H

#include <QMainWindow>

namespace Ui {
class ProjectView;
}

namespace ui {
namespace task {


class ProjectView : public QMainWindow
{
    Q_OBJECT

public:
    explicit ProjectView(QWidget *parent = 0);
    ~ProjectView();

private:
    Ui::ProjectView *ui;
};

} // namespace task
} // namespace ui


#endif // PROJECT_VIEW_H
