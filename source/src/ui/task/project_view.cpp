#include "project_view.h"
#include "ui_project_view.h"

#include <ui/task/project_model.h>


namespace ui {
namespace task {


ProjectView::ProjectView(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::ProjectView)
{
    ui->setupUi(this);

	ProjectModel *model = new ProjectModel();
        ui->project->setModel(model);

}

ProjectView::~ProjectView()
{
    delete ui;
}

} // namespace task
} // namespace ui
