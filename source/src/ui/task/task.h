#ifndef TASK_H
#define TASK_H

#include <ui/task/node.h>

namespace ui {
namespace task {

class Task : public Node
{
public:
	explicit Task(QUuid _node_id, Node *_parent);

	std::string type() const override;

};

} // namespace task
} // namespace ui


#endif // TASK_H
