/// @brief inherit from protocols/rigid/RigidBodyPerturbNoCenterMover aimed to set it free from DockSetupMover, which was originally required to set RigidBodyPerturbNoCenter's jump. Intead of that, we store the RigidBodyInfo into basic::datacache::DataMap in DockSetupMover, then later RigidBodyInfo could be retrieved whenever it is needed (i.e. DockingInitialPerturbation) during docking. This is especially useful in TempWeightedReplica for docking in which we have multiple RigidBodyPerturbNoCenterMover for each replica (each correspond to a unique temperature), it is not possible to use DockSetupMover to set the jump for each RigidBodyPerturbNoCenterMover.

/// @author Zhe Zhang

#ifndef INCLUDED_devel_replica_docking_UnbiasedRigidBodyMover_hh
#define INCLUDED_devel_replica_docking_UnbiasedRigidBodyMover_hh

#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <protocols/canonical_sampling/HamiltonianExchange.hh>
#include <protocols/canonical_sampling/HamiltonianExchange.fwd.hh>

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

class UnbiasedRigidBodyPerturbNoCenterMover : public protocols::rigid::RigidBodyPerturbNoCenterMover {
	typedef RigidBodyPerturbNoCenterMover Parent;
	//  typedef std::map< std::string, devel::replica_docking::TempInterpolatorOP > Interpolators;
	//  typedef utility::vector1< protocols::canonical_sampling::UnbiasedMoverOP > MoverOPs;

public:
	UnbiasedRigidBodyPerturbNoCenterMover();
	UnbiasedRigidBodyPerturbNoCenterMover( UnbiasedRigidBodyPerturbNoCenterMover const & );

	~UnbiasedRigidBodyPerturbNoCenterMover() override;

	/// @brief overload it to use random unit quaternion to unbiasedly sample rotation instead of Jump::gaussian_move
	void apply( core::pose::Pose & pose ) override;

	std::string get_name() const override;

	protocols::moves::MoverOP clone() const override;

	protocols::moves::MoverOP fresh_instance() const override;

	void initialize( core::pose::Pose const &);

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	//   virtual void
	//   initialize_simulation(
	//      core::pose::Pose& pose,
	//      protocols::canonical_sampling::MetropolisHastingsMover const& metropolis_hastings_mover,
	//      core::Size cycle
	//   );

private:

	protocols::docking::RigidBodyInfoOP rigid_body_info_;
	protocols::docking::DockJumps movable_jumps_;
	bool initialized_;
	bool restrict_; // if restrict the searching space from a reference structure
	bool max_move_;
	std::string ref_file_;
	core::Real max_trans_dist_;
	core::Real max_rot_angle_;
	numeric::xyzVector<core::Real> ref_T_;
	numeric::xyzMatrix<core::Real> ref_R_;

};

} // namespace replica_docking
}

#endif
