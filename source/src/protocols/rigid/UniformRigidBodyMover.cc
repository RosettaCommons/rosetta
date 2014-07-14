// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rigid/UniformRigidBodyMover.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/rigid/UniformRigidBodyMover.hh>
#include <protocols/rigid/UniformRigidBodyMoverCreator.hh>

// Package headers

// Project headers

//Utility Headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random.functions.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace rigid {

int const NO_JUMP = 0;

static basic::Tracer tr("protocols.rigid.UniformRigidBodyMover", basic::t_info);
static numeric::random::RandomGenerator RG(434855);

using namespace core::environment;
using namespace protocols::environment;

// creator
std::string
UniformRigidBodyMoverCreator::keyname() const {
  return UniformRigidBodyMoverCreator::mover_name();
}

protocols::moves::MoverOP
UniformRigidBodyMoverCreator::create_mover() const {
  return new UniformRigidBodyMover;
}

std::string
UniformRigidBodyMoverCreator::mover_name() {
  return "UniformRigidBodyMover";
}

UniformRigidBodyMover::UniformRigidBodyMover():
  ThermodynamicMover(),
  target_jump_( NO_JUMP ),
  rotation_mag_( 3.0 ),
  translation_mag_( 8.0 )
{}

UniformRigidBodyMover::UniformRigidBodyMover( JumpNumber target_jump,
                                              core::Real rotation_mag,
                                              core::Real translation_mag ):
  ThermodynamicMover(),
  target_jump_( target_jump ),
  rotation_mag_( rotation_mag ),
  translation_mag_( translation_mag )
{}

void UniformRigidBodyMover::apply( core::pose::Pose& pose ){
  using namespace numeric;
  using namespace numeric::random;

  if( target_jump_ == NO_JUMP ){
    std::ostringstream ss;
    ss << "The target jump number of " << this->get_name()
       << " was " << NO_JUMP << ", which probably means this value wasn't set properly.  "
       << "If you're using RosettaScripts, try using the 'target_jump' option." << std::endl;
    throw utility::excn::EXCN_BadInput( ss.str() );
  }

  core::kinematics::Jump flexible_jump = pose.jump( target_jump_ );

  xyzMatrix< core::Real> const rot=flexible_jump.get_rotation();
  xyzVector< core::Real> const trans=flexible_jump.get_translation();

  xyzVector<core::Real> delta_trans = random_translation( translation_mag_, RG );

  core::Real theta = random_rotation_angle<core::Real>( rotation_mag_, RG );

  xyzVector<core::Real> axis = random_point_on_unit_sphere<core::Real>( RG );

  xyzMatrix<core::Real> delta_rot = rotation_matrix_radians( axis, theta );

  flexible_jump.set_translation( delta_trans + trans );
  flexible_jump.set_rotation( delta_rot*rot );

  pose.set_jump( (int) target_jump_, flexible_jump );
}

void UniformRigidBodyMover::jump_number( JumpNumber jnum ) {
  target_jump_ = jnum;
}

UniformRigidBodyMover::JumpNumber UniformRigidBodyMover::jump_number() const {
  return target_jump_;
}

void UniformRigidBodyMover::parse_my_tag( utility::tag::TagCOP tag,
                                          basic::datacache::DataMap&,
                                          protocols::filters::Filters_map const&,
                                          protocols::moves::Movers_map const&,
                                          core::pose::Pose const& ) {
  rotation_mag_ = tag->getOption< core::Real >( "rotation_magnitude", 3.0 );
  translation_mag_ = tag->getOption< core::Real >( "translation_magnitude", 8.0 );
  target_jump_ = tag->getOption< JumpNumber >("target_jump", NO_JUMP );
}

std::string UniformRigidBodyMover::get_name() const {
    return "UniformRigidBodyMover";
}

utility::vector1<core::id::TorsionID_Range>
UniformRigidBodyMover::torsion_id_ranges( core::pose::Pose & ) {
  return utility::vector1<core::id::TorsionID_Range>();
}

moves::MoverOP UniformRigidBodyMover::fresh_instance() const {
  return new UniformRigidBodyMover();
}

moves::MoverOP UniformRigidBodyMover::clone() const{
  return new UniformRigidBodyMover( *this );
}

} // rigid
} // protocols
