// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  inheried from RigidBodyPerturbNoCenterMover, 1. use Rigid_Body_Info to parse the movable_jump to
/// get independant from DockSetupMover; 2. use unbiased rotation sampling and gaussing translation instead
/// of gaussian_move() in which rotation is not unbiased.
/// @author Zhe Zhang


#include <devel/replica_docking/UnbiasedRigidBodyMover.hh>
#include <devel/replica_docking/UnbiasedRigidBodyMoverCreator.hh>

#include <protocols/docking/RigidBodyInfo.hh>
#include <basic/datacache/DataMap.hh>

#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Jump.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/AtomID.hh>

//#include <numeric/rotation.functions.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/trig.functions.hh>

#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

static basic::Tracer tr( "devel.replica_docking.UnbiasedRigidBodyMover" );
static numeric::random::RandomGenerator rigid_RG(24205385);

namespace devel {
namespace replica_docking {

std::string
UnbiasedRigidBodyPerturbNoCenterMoverCreator::keyname() const {
  return UnbiasedRigidBodyPerturbNoCenterMoverCreator::mover_name();
}

protocols::moves::MoverOP
UnbiasedRigidBodyPerturbNoCenterMoverCreator::create_mover() const {
  return new UnbiasedRigidBodyPerturbNoCenterMover;
}

std::string
UnbiasedRigidBodyPerturbNoCenterMoverCreator::mover_name() {
  return "UnbiasedRigidBodyPerturbNoCenter";
}

UnbiasedRigidBodyPerturbNoCenterMover::UnbiasedRigidBodyPerturbNoCenterMover() :
  rigid_body_info_( NULL )
{}

UnbiasedRigidBodyPerturbNoCenterMover::UnbiasedRigidBodyPerturbNoCenterMover( UnbiasedRigidBodyPerturbNoCenterMover const& other ) : RigidBodyPerturbNoCenterMover( other ) {
  ///copy value of every private variables
  rigid_body_info_ = other.rigid_body_info_;
  movable_jumps_ = other.movable_jumps_;
}

UnbiasedRigidBodyPerturbNoCenterMover::~UnbiasedRigidBodyPerturbNoCenterMover() {}

std::string
UnbiasedRigidBodyPerturbNoCenterMover::get_name() const
{
  return "UnbiasedRigidBodyPerturbNoCenter";
}

protocols::moves::MoverOP
UnbiasedRigidBodyPerturbNoCenterMover::clone() const
{
  return new UnbiasedRigidBodyPerturbNoCenterMover(*this);
}

protocols::moves::MoverOP
UnbiasedRigidBodyPerturbNoCenterMover::fresh_instance() const
{
  return new UnbiasedRigidBodyPerturbNoCenterMover;
}

void
UnbiasedRigidBodyPerturbNoCenterMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
) {


  if ( !data.has( "RigidBodyInfo", "docking_setup" ) ) {
    tr << "RigidBodyInfo not found in basic::datacache::DataMap" << std::endl;
    rigid_body_info_ = new protocols::docking::RigidBodyInfo;
    data.add( "RigidBodyInfo", "docking_setup", rigid_body_info_ );
    //	     throw utility::excn::EXCN_RosettaScriptsOption( "RigidBodyInfo not found in basic::datacache::DataMap, DockingInitialPerturbation can not be done, so exit here!" );
  } else {
    rigid_body_info_ = data.get< protocols::docking::RigidBodyInfo* >( "RigidBodyInfo", "docking_setup" );
    tr.Debug << "get RigidBodyInfo pointer from basic::datacache::DataMap" << std::endl;
  }

  Parent::parse_my_tag( tag, data, filters, movers, pose );
}


void UnbiasedRigidBodyPerturbNoCenterMover::initialize_simulation(
  core::pose::Pose & pose,
  protocols::canonical_sampling::MetropolisHastingsMover const & mhm,
  core::Size cycle
) {
  tr.Debug << "initialize_simulation: " << std::endl;
  if ( rigid_body_info_ ) {
    movable_jumps_ = rigid_body_info_->movable_jumps();
    tr.Debug << "finished reading movable_jumps_ from RigidBodyInfo" << std::endl;
    if ( movable_jumps_.empty() ) {
      utility_exit_with_message( "DockSetupMover has to be applied before this !" );
    }
  }

  runtime_assert( !movable_jumps_.empty() );
  for ( protocols::docking::DockJumps::const_iterator it=movable_jumps_.begin(); it != movable_jumps_.end(); ++it ) {
    Parent::add_jump( *it );
    tr.Debug << "movable_jump " << *it << " was added!" << std::endl;
  }
  Parent::initialize_simulation( pose, mhm, cycle );
}

void UnbiasedRigidBodyPerturbNoCenterMover::apply( core::pose::Pose& pose ) {
  // Parent::apply( pose );
  // overload to use random unit quaternion to generate uniformly distributed rotation.

  // set baseclass rb_jump_ randomly from list of movable jumps
  if ( movable_jumps_.size() > 1 ) {
    rb_jump_ = rigid_RG.random_element( movable_jumps_ );
  } else if ( movable_jumps_.size() == 1 ) {
    rb_jump_ = movable_jumps_[1];
    tr.Debug <<"set rb_jump_ movable_jumps_[1]" << std::endl;
  } else {
    rb_jump_ = 1;
    tr.Debug <<"set rb_jump_ 1 " << std::endl;
  }

  tr.Debug << "Set movable jump# " << rb_jump_ << std::endl;
  tr.Debug <<"rot_mag "<< rot_mag_ <<" ; trans_mag " << trans_mag_ << std::endl;

  core::kinematics::Jump flexible_jump = pose.jump( rb_jump_ );
  numeric::xyzMatrix< core::Real> const rot=flexible_jump.get_rotation();
  numeric::xyzVector< core::Real> const trans=flexible_jump.get_translation();

  numeric::xyzVector<core::Real> delta_trans/*(0.0245597,-0.181304,-0.276021);*/( trans_mag_*rigid_RG.gaussian(), trans_mag_*rigid_RG.gaussian(), trans_mag_*rigid_RG.gaussian() );
  core::Real theta = numeric::random::random_rotation_angle( rot_mag_, rigid_RG );
  numeric::xyzVector<core::Real> axis = numeric::random::random_point_on_unit_sphere< core::Real >( rigid_RG );
	tr.Info << "axis=["<<axis.x() <<" "<<axis.y()<<" "<<axis.z()<<"]"<<std::endl;
  numeric::xyzMatrix<core::Real> delta_rot = numeric::rotation_matrix_radians( axis, theta );
  flexible_jump.set_translation( delta_trans + trans );
  flexible_jump.set_rotation( delta_rot*rot );
  //  flexible_jump.gaussian_trans_random_rotation( dir_, trans_mag_, rot_mag_);
  pose.set_jump( rb_jump_, flexible_jump );


}

}
}
