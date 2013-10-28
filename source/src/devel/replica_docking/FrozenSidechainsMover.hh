
#ifndef INCLUDED_devel_replica_docking_FrozenSidechainsMover_hh
#define INCLUDED_devel_replica_docking_FrozenSidechainsMover_hh

#include <devel/replica_docking/FrozenSidechainsMover.fwd.hh>
//#include <protocols/docking/types.hh>

//#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>

#include <protocols/moves/Mover.hh>
//#include <protocols/filters/Filter.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/scoring/InterfaceInfo.fwd.hh>
#include <protocols/simple_moves/ReturnSidechainMover.fwd.hh>


namespace devel {
namespace replica_docking {

class FrozenSidechainsMover : public protocols::moves::Mover {

public:

  FrozenSidechainsMover();

  //  virtual ~AddEncounterConstraintMover();

  virtual protocols::moves::MoverOP clone() const;

  virtual void apply( core::pose::Pose & pose );
  void compute_interface( core::pose::Pose & pose ) const;

  protocols::scoring::InterfaceInfo const & interface_from_pose( core::pose::Pose const & ) const;
  protocols::scoring::InterfaceInfo & nonconst_interface_from_pose( core::pose::Pose & ) const;

  virtual std::string get_name() const;

  void parse_my_tag(
		    utility::tag::TagCOP const tag,
		    basic::datacache::DataMap &,
		    protocols::filters::Filters_map const &,
		    protocols::moves::Movers_map const &,
		    core::pose::Pose const &
  );

private:

  utility::vector1<bool> allow_copy_chi_;
  protocols::simple_moves::ReturnSidechainMoverOP frozen_sidechains_;
  core::pose::Pose recover_sidechains_pose_;

};

}
}

#endif
