#pragma once

#include <QPointer>

#include <memory>

namespace ui {
namespace task {

class TaskSyncer_NodeStrategy;

using TaskSyncer_NodeStrategyQP  = QPointer<TaskSyncer_NodeStrategy>;

using TaskSyncer_NodeStrategySP  = std::shared_ptr<TaskSyncer_NodeStrategy>;
using TaskSyncer_NodeStrategyCSP = std::shared_ptr<TaskSyncer_NodeStrategy const >;

} // namespace task
} // namespace ui
