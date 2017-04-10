#ifndef UI_TASK_PROJECT_H
#define UI_TASK_PROJECT_H

#include <ui/task/node.h>


namespace ui {
namespace task {

class Project : public Node
{
public:
	explicit Project(QUuid _project_id);

	std::string type() const override;


private:
	QUuid const project_id;
};

} // namespace task
} // namespace ui

#endif // UI_TASK_PROJECT_H
