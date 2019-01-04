
#ifndef INCLUDED_protocols_replica_docking_FrozenSidechainsCreator_hh
#define INCLUDED_protocols_replica_docking_FrozenSidechainsCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace replica_docking {

class FrozenSidechainsMoverCreator : public protocols::moves::MoverCreator {
public:
  protocols::moves::MoverOP create_mover() const override;
  std::string keyname() const override;
  static std::string mover_name();

};

}
}
#endif
