#include "task.h"

namespace ui {
namespace task {


Task::Task(QUuid _node_id, Node *_parent) : Node(_node_id, _parent)
{

}


std::string Task::type() const
{
	return "task";
}


} // namespace task
} // namespace ui
