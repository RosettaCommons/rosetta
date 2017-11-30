#ifndef TASK_TASK_FWD_H
#define TASK_TASK_FWD_H

#include <QPointer>
#include <memory>

namespace ui {
namespace task {

class Task;
using TaskSP  = std::shared_ptr< Task >;
using TaskCSP = std::shared_ptr< Task const >;
using TaskQP  = QPointer<Task>;

class File;
using FileSP  = std::shared_ptr< File >;
using FileCSP = std::shared_ptr< File const >;
using FileQP  = QPointer<File>;

} // namespace task
} // namespace ui


#endif // TASK_TASK_H
