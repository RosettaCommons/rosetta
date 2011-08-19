#ifndef INCLUDED_devel_optest_A_fwd_hh
#define INCLUDED_devel_optest_A_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace optest {

class A;

typedef utility::pointer::owning_ptr< A > AOP;
typedef utility::pointer::owning_ptr< A const > ACOP;

}
}

#endif
