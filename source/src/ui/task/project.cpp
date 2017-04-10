#include <ui/task/project.h>

namespace ui {
namespace task {


Project::Project(QUuid _project_id) : Node(QUuid(), nullptr), project_id(_project_id)
{
}


std::string Project::type() const
{
	return "project";
}


} // namespace task
} // namespace ui
