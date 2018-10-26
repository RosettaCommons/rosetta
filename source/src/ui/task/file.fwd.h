#pragma once

#include <QPointer>
#include <memory>

namespace ui {
namespace task {

struct FileID;

class File;
using FileSP  = std::shared_ptr< File >;
using FileCSP = std::shared_ptr< File const >;
using FileQP  = QPointer<File>;

} // namespace task
} // namespace ui
