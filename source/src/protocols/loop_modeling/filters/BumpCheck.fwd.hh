#pragma once

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loop_modeling {
namespace filters {

class BumpCheck;

typedef utility::pointer::owning_ptr<BumpCheck> BumpCheckOP;
typedef utility::pointer::owning_ptr<BumpCheck const> BumpCheckCOP;

}
}
}

