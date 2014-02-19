/// @brief inherit from protocols/rigid/RigidBodyPerturbNoCenterMover aimed to set it free from DockSetupMover, which was originally required to set RigidBodyPerturbNoCenter's jump. Intead of that, we store the RigidBodyInfo into basic::datacache::DataMap in DockSetupMover, then later RigidBodyInfo could be retrieved whenever it is needed (i.e. DockingInitialPerturbation) during docking. This is especially useful in TempWeightedReplica for docking in which we have multiple RigidBodyPerturbNoCenterMover for each replica (each correspond to a unique temperature), it is not possible to use DockSetupMover to set the jump for each RigidBodyPerturbNoCenterMover.

/// @author Zhe Zhang

#ifndef INCLUDED_devel_replica_docking_ThermodynamicRigidBodyMover_hh
#define INCLUDED_devel_replica_docking_ThermodynamicRigidBodyMover_hh

#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <protocols/canonical_sampling/MpiHamiltonianExchange.hh>
#include <protocols/canonical_sampling/MpiHamiltonianExchange.fwd.hh>

#include <devel/replica_docking/TempInterpolator.hh>
#include <devel/replica_docking/TempInterpolator.fwd.hh>

#include <protocols/rigid/RigidBodyMover.fwd.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/docking/RigidBodyInfo.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <protocols/docking/types.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

#include <utility/pointer/ReferenceCount.hh>


namespace devel {
namespace replica_docking {

class ThermodynamicRigidBodyPerturbNoCenterMover : public protocols::rigid::RigidBodyPerturbNoCenterMover {
  typedef RigidBodyPerturbNoCenterMover Parent;
  //  typedef std::map< std::string, devel::replica_docking::TempInterpolatorOP > Interpolators;
  //  typedef utility::vector1< protocols::canonical_sampling::ThermodynamicMoverOP > MoverOPs;

public:
  ThermodynamicRigidBodyPerturbNoCenterMover();
  ThermodynamicRigidBodyPerturbNoCenterMover( ThermodynamicRigidBodyPerturbNoCenterMover const & );

  virtual ~ThermodynamicRigidBodyPerturbNoCenterMover();

  ///@brief overload it to use random unit quaternion to unbiasedly sample rotation instead of Jump::gaussian_move
  virtual void apply( core::pose::Pose & pose );

  virtual std::string get_name() const;

  protocols::moves::MoverOP clone() const;

  virtual protocols::moves::MoverOP fresh_instance() const;

  virtual void parse_my_tag(
       utility::tag::TagCOP tag,
       basic::datacache::DataMap &,
       protocols::filters::Filters_map const &,
       protocols::moves::Movers_map const &,
       core::pose::Pose const &
  );

  virtual void
  initialize_simulation(
     core::pose::Pose& pose,
     protocols::canonical_sampling::MetropolisHastingsMover const& metropolis_hastings_mover,
     core::Size cycle
  );

  core::Real generate_rotation_angle( core::Real );
  numeric::xyzVector< core::Real > generate_uniform_rotation_axis();
private:

  protocols::docking::RigidBodyInfoOP rigid_body_info_;
  protocols::docking::DockJumps movable_jumps_;

};

} // namespace replica_docking
}

#endif
