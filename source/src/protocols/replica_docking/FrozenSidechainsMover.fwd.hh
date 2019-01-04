
#ifndef INCLUDED_protocols_replica_docking_FrozenSidechainsMover_fwd_hh
#define INCLUDED_protocols_replica_docking_FrozenSidechainsMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace replica_docking{

class FrozenSidechainsMover;
typedef utility::pointer::owning_ptr< FrozenSidechainsMover > FrozenSidechainsMoverOP;
typedef utility::pointer::owning_ptr< FrozenSidechainsMover const > FrozenSidechainsMoverCOP;

}
}
#endif
