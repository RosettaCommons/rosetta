#ifndef INCLUDED_devel_optest_B1_fwd_hh
#define INCLUDED_devel_optest_B1_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace optest {

class B1;

typedef utility::pointer::shared_ptr< B1 > B1OP;
typedef utility::pointer::shared_ptr< B1 const > B1COP;

}
}

#endif
