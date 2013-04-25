#ifndef INCLUDED_protocols_docking_DockingInitialPerturbationCreator_hh
#define INCLUDED_protocols_docking_DockingInitialPerturbationCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace docking {

class DockingInitialPerturbationCreator : public protocols::moves::MoverCreator {

public:
  virtual protocols::moves::MoverOP create_mover() const;
  virtual std::string keyname() const;
  static std::string mover_name();

};

} // docking
} // protocols
#endif
