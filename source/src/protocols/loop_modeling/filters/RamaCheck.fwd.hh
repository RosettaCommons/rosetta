#pragma once

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loop_modeling {
namespace filters {

class RamaCheck;

typedef utility::pointer::owning_ptr<RamaCheck> RamaCheckOP;
typedef utility::pointer::owning_ptr<RamaCheck const> RamaCheckCOP;

}
}
}

