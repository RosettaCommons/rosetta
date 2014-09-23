#ifndef INCLUDED_devel_optest_B2_fwd_hh
#define INCLUDED_devel_optest_B2_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace optest {

class B2;

typedef utility::pointer::shared_ptr< B2 > B2OP;
typedef utility::pointer::shared_ptr< B2 const > B2COP;

}
}

#endif
