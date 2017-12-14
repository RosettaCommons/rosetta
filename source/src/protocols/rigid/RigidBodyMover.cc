// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file docking_initialization_protocols
/// @brief initialization protocols for docking
/// @details
///  This contains the functions that create initial positions for docking
///  You can either randomize partner 1 or partner 2, spin partner 2, or
///  perform a simple perturbation.
/// @author Monica Berrondo

// Unit headers
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RigidBodyMoverCreator.hh>

// Package headers
#include <protocols/rigid/RB_geometry.hh>
#include <core/pose/PDBInfo.hh>
// Rosetta Headers
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility>
#include <utility/tag/Tag.hh>

// Random number generator
#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/TorsionID_Range.hh>

// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/xyz.functions.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <string>
#include <core/id/TorsionID_Range.hh>
#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <numeric/xyz.functions.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/options/IntegerVectorOption.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


using basic::Error;
using basic::Warning;

namespace protocols {
namespace rigid {

using namespace core;

static basic::Tracer TR( "protocols.rigid.RigidBodyMover" );
static basic::Tracer TRBM( "protocols.moves.RigidBodyMover" );

// Large rotational perturbations produce a distribution of orientations
// that are neither uniform nor similar to the input orientation.
// What we probably wanted was a random orientation (possibly plus a translation).
const Real max_allowed_rot_mag ( 60.0 );


// default constructor
RigidBodyMover::RigidBodyMover() :
	protocols::canonical_sampling::ThermodynamicMover(),
	rb_jump_( 1 ), dir_( n2c ), rot_center_( 0.0 ), freeze_(false)
{
	Mover::type( "RigidBodyBase" );
}

// constructor with arguments
RigidBodyMover::RigidBodyMover(
	int const rb_jump_in,
	Direction dir_in
):
	protocols::canonical_sampling::ThermodynamicMover(),
	rb_jump_( rb_jump_in ), dir_( dir_in ), rot_center_( 0.0 ), freeze_(false)
{
	Mover::type( "RigidBodyBase" );
	if ( dir_ == random ) {
		dir_ = ( numeric::random::uniform() < 0.5 ? c2n : n2c );
	} else {
		runtime_assert( dir_ == n2c || dir_ == c2n );
	}
}

RigidBodyMover::RigidBodyMover( RigidBodyMover const & /*src*/ ) = default;

RigidBodyMover::~RigidBodyMover() = default;

std::string
RigidBodyMover::get_name() const {
	return "RigidBodyMover";
}

void
RigidBodyMover::show(std::ostream & output) const
{
	Mover::show(output);
	output << "Jump number:  " << rb_jump() << std::endl;
}

utility::vector1<core::id::TorsionID_Range>
RigidBodyMover::torsion_id_ranges(
	core::pose::Pose & //pose
)
{
	return utility::vector1<core::id::TorsionID_Range>();
}


RigidBodyPerturbMover::RigidBodyPerturbMover(
	int const rb_jump_in,
	core::Real const rot_mag_in,
	core::Real const trans_mag_in,
	Partner const partner_in,
	bool interface_in      //rot_center calculated at interface
):
	RigidBodyMover( rb_jump_in ),
	rot_mag_( rot_mag_in ),
	trans_mag_( trans_mag_in ),
	partner_( partner_in ),
	interface_(interface_in)
{
	movable_jumps_.push_back( rb_jump_in );
	TRBM.Trace << "rb_jump " << rb_jump_in << std::endl;
	TRBM.Trace << "rot_mag " << rot_mag_ << std::endl;
	TRBM.Trace << "trans_mag " << trans_mag_ << std::endl;
	Mover::type( "RigidBodyPerturb" );
}

RigidBodyPerturbMover::RigidBodyPerturbMover() :
	RigidBodyMover(),
	rot_mag_( 3.0 ),
	trans_mag_( 8.0 )
{
	Mover::type( "RigidBodyPerturb" );
}

RigidBodyPerturbMover::RigidBodyPerturbMover(
	int const rb_jump_in,
	core::Real const rot_mag_in,
	core::Real const trans_mag_in,
	Partner const partner_in,
	utility::vector1< bool > const & ok_for_centroid_calculation
):
	RigidBodyMover( rb_jump_in ),
	rot_mag_( rot_mag_in ),
	trans_mag_( trans_mag_in ),
	partner_( partner_in ),
	interface_( false ),
	ok_for_centroid_calculation_( ok_for_centroid_calculation )
{
	movable_jumps_.push_back( rb_jump_in );
	TRBM.Trace << "rb_jump " << rb_jump_in << std::endl;
	TRBM.Trace << "rot_mag " << rot_mag_ << std::endl;
	TRBM.Trace << "trans_mag " << trans_mag_ << std::endl;
	Mover::type( "RigidBodyPerturb" );
}

RigidBodyPerturbMover::RigidBodyPerturbMover(
	core::pose::Pose const & pose_in,
	core::kinematics::MoveMap const & mm,
	core::Real const rot_mag_in,
	core::Real const trans_mag_in,
	Partner const partner_in,
	bool interface_in //rot_center calculated at interface
):
	RigidBodyMover(),
	rot_mag_( rot_mag_in ),
	trans_mag_( trans_mag_in ),
	partner_( partner_in ),
	interface_(interface_in)
{
	TRBM.Trace << "rot_mag " << rot_mag_ << std::endl;
	TRBM.Trace << "trans_mag " << trans_mag_ << std::endl;
	Mover::type( "RigidBodyPerturb" );
	for ( Size i=1, i_end = pose_in.num_jump(); i<= i_end; ++i ) {
		if ( mm.get_jump(i) ) {
			movable_jumps_.push_back( i );
		}
	}

	if ( movable_jumps_.empty() ) {
		TRBM.Warning << "no movable jumps!" << std::endl;
		return;
	}
}

RigidBodyPerturbMover::RigidBodyPerturbMover(
	core::Real const rot_mag_in,
	core::Real const trans_mag_in,
	Partner const partner_in,
	bool interface_in      //rot_center calculated at interface
):
	RigidBodyMover(),
	rot_mag_( rot_mag_in ),
	trans_mag_( trans_mag_in ),
	partner_( partner_in ),
	interface_(interface_in)
{
	TRBM.Trace << "rb_jump " << rb_jump_ << std::endl;
	TRBM.Trace << "rot_mag " << rot_mag_ << std::endl;
	TRBM.Trace << "trans_mag " << trans_mag_ << std::endl;
	Mover::type( "RigidBodyPerturb" );
}

RigidBodyPerturbMover::RigidBodyPerturbMover( RigidBodyPerturbMover const & src ) :
	//utility::pointer::ReferenceCount(),
	parent( src ),
	rot_mag_( src.rot_mag_ ),
	trans_mag_( src.trans_mag_ ),
	partner_( src.partner_ ),
	interface_( src.interface_ ),
	movable_jumps_( src.movable_jumps_ )
{}

RigidBodyPerturbMover::~RigidBodyPerturbMover() = default;

void
RigidBodyPerturbMover::apply( core::pose::Pose & pose )
{
	// Want to update our center of rotation every time we take a step.
	// baseclass jump is chosen at random from movable jumps every apply
	rb_jump_ = ( movable_jumps_.size() > 1 ) ?
		numeric::random::rg().random_element( movable_jumps_ )
		: movable_jumps_[1];

	TR.Trace << "Set movable jump #" << rb_jump_ << std::endl;

	// Set center of rotation unless we are frozen.
	if ( !freeze_ ) {
		core::Vector dummy_up, dummy_down;
		if ( interface_ ) {
			protocols::geometry::centroids_by_jump_int(pose, rb_jump_, dummy_up, dummy_down);
		} else {
			protocols::geometry::centroids_by_jump(pose, rb_jump_, dummy_up, dummy_down, ok_for_centroid_calculation_ );
		}
		rot_center_ = ( partner_ == partner_downstream ) ? dummy_down : dummy_up;
	}

	core::kinematics::Jump flexible_jump = pose.jump( rb_jump_ );
	core::kinematics::Stub downstream_stub = pose.conformation().downstream_jump_stub( rb_jump_ );
	flexible_jump.set_rb_center( dir_, downstream_stub, rot_center_ );

	// this should probably be changed so that instead of doing this here, it calls the randomize apply
	if ( rot_mag_ >= max_allowed_rot_mag ) {
		TR.Warning << "Large Gaussian rotational perturbations don't make sense!  Bad choices with -dock_pert?  Use -randomize[12] instead." << std::endl;
	}

	if ( !freeze_ ) {
		rb_delta_ = flexible_jump.gaussian_move( dir_, trans_mag_, rot_mag_ );
	} else {
		flexible_jump.set_rb_deltas(dir_, rb_delta_);
		flexible_jump.fold_in_rb_deltas();
	}
	pose.set_jump( rb_jump_, flexible_jump );
} // RigidBodyPerturbMover::apply()

std::string
RigidBodyPerturbMover::get_name() const {
	return "RigidBodyPerturbMover";
}

void
RigidBodyPerturbMover::show(std::ostream & output) const
{
	RigidBodyMover::show(output);
	output << "Magnitude of translational movement (deg): " << get_trans_mag() << std::endl <<
		"Magnitude of rotational movement (deg):    " << get_rot_mag() << std::endl;
}

core::Distance
RigidBodyPerturbMover::get_trans_mag() const {
	return trans_mag_;
}

core::Angle
RigidBodyPerturbMover::get_rot_mag() const {
	return rot_mag_;
}

void
RigidBodyPerturbMover::rot_center( core::Vector const & /*rot_center_in*/ )
{
	utility_exit_with_message("Rotation point is automatically determined ONLY");
}

std::ostream
&operator<< ( std::ostream &os, RigidBodyPerturbMover const &mover )
{
	mover.show(os);
	return os;
}

//Implementation of RigidBodyPerturbRandomJumpMover, which takes a random jump and calls RigidBodyPerturbMover
RigidBodyPerturbRandomJumpMover::RigidBodyPerturbRandomJumpMover() : rot_mag_in_(3.0), trans_mag_in_(8.0), num_jump_(0){}
RigidBodyPerturbRandomJumpMover::~RigidBodyPerturbRandomJumpMover()= default;

RigidBodyPerturbRandomJumpMover::RigidBodyPerturbRandomJumpMover(
	core::Real const& rot_mag_in,
	core::Real const& trans_mag_in,
	core::Size const& num_jump_in) :
	rot_mag_in_(rot_mag_in),
	trans_mag_in_(trans_mag_in),
	num_jump_(num_jump_in)
{}

void
RigidBodyPerturbRandomJumpMover::apply(core::pose::Pose& pose)
{
	auto random_jump_num = static_cast<core::Size>(numeric::random::rg().random_range(1,num_jump_));
	RigidBodyPerturbMover RBMover(random_jump_num,rot_mag_in_,trans_mag_in_);
	RBMover.apply(pose);
}

std::string
RigidBodyPerturbRandomJumpMover::get_name() const {
	return "RigidBodyPerturbRandomJumpMover";
}

void
RigidBodyPerturbNoCenterMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap&,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
	rot_mag_ = tag->getOption< core::Real >( "rot_mag", 0.1 );
	trans_mag_ = tag->getOption< core::Real >( "trans_mag", 0.4 );
}

// XRW TEMP std::string
// XRW TEMP RigidBodyPerturbNoCenterMoverCreator::keyname() const {
// XRW TEMP  return RigidBodyPerturbNoCenterMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RigidBodyPerturbNoCenterMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new RigidBodyPerturbNoCenterMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RigidBodyPerturbNoCenterMover::mover_name() {
// XRW TEMP  return "RigidBodyPerturbNoCenter";
// XRW TEMP }

moves::MoverOP
RigidBodyPerturbNoCenterMover::clone() const {
	return moves::MoverOP( new RigidBodyPerturbNoCenterMover(*this) );
}

RigidBodyPerturbNoCenterMover::RigidBodyPerturbNoCenterMover() :
	Parent(),
	rot_mag_( 3.0 ),
	trans_mag_( 8.0 )
{
	moves::Mover::type( "RigidBodyPerturbNoCenter" );
}

RigidBodyPerturbNoCenterMover::RigidBodyPerturbNoCenterMover(
	int const rb_jump_in,
	core::Real const rot_mag_in,
	core::Real const trans_mag_in
) : RigidBodyMover(),
	rot_mag_( rot_mag_in ),
	trans_mag_( trans_mag_in )
{
	movable_jumps_.push_back( rb_jump_in );
	moves::Mover::type( "RigidBodyPerturbNoCenter" );
}

/// @details constructor for the rbm that doesn't set a center
RigidBodyPerturbNoCenterMover::RigidBodyPerturbNoCenterMover(
	core::pose::Pose const & pose_in,
	kinematics::MoveMap const & mm,
	Real const rot_mag_in,
	Real const trans_mag_in,
	Direction dir_in
):
	RigidBodyMover(),
	rot_mag_( rot_mag_in ),
	trans_mag_( trans_mag_in )
{

	TRBM.Debug << "rot_mag " << rot_mag_in << std::endl;
	TRBM.Debug << "trans_mag " << trans_mag_in << std::endl;
	moves::Mover::type( "RigidBodyPerturbNoCenter" );
	for ( Size i=1, i_end = pose_in.num_jump(); i<= i_end; ++i ) {
		if ( mm.get_jump(i) ) {
			if ( std::find(movable_jumps_.begin(), movable_jumps_.end(), i) == movable_jumps_.end() ) { // if jump is not already in the list
				movable_jumps_.push_back( i );
			}
		}
	}

	if ( movable_jumps_.empty() ) {
		TRBM.Warning << "no movable jumps!" << std::endl;
		return;
	}
	dir_ = dir_in;
}

RigidBodyPerturbNoCenterMover::RigidBodyPerturbNoCenterMover( RigidBodyPerturbNoCenterMover const & /*src*/ ) = default;

RigidBodyPerturbNoCenterMover::~RigidBodyPerturbNoCenterMover() = default;

void
RigidBodyPerturbNoCenterMover::add_jump( core::Size jump_id ) {
	movable_jumps_.push_back( jump_id );
}

void
RigidBodyPerturbNoCenterMover::clear_jumps() {
	movable_jumps_.clear();
}

void
RigidBodyPerturbNoCenterMover::apply( core::pose::Pose & pose )
{
	// set baseclass rb_jump_ randomly from list of movable jumps
	if ( movable_jumps_.size() > 1 ) {
		rb_jump_ = numeric::random::rg().random_element( movable_jumps_ );
	} else if ( movable_jumps_.size() == 1 ) {
		rb_jump_ = movable_jumps_[1];
	} else {
		rb_jump_ = 1;
	}

	TR.Debug << "Set movable jump# " << rb_jump_ << std::endl;
	core::kinematics::Jump flexible_jump = pose.jump( rb_jump_ );
	flexible_jump.gaussian_move( dir_, trans_mag_, rot_mag_ );
	pose.set_jump( rb_jump_, flexible_jump );
}

// XRW TEMP std::string
// XRW TEMP RigidBodyPerturbNoCenterMover::get_name() const {
// XRW TEMP  return "RigidBodyPerturbNoCenterMover";
// XRW TEMP }

std::string RigidBodyPerturbNoCenterMover::get_name() const {
	return mover_name();
}

std::string RigidBodyPerturbNoCenterMover::mover_name() {
	return "RigidBodyPerturbNoCenter";
}

void RigidBodyPerturbNoCenterMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute::attribute_w_default("rot_mag", xsct_real, "XRW TO DO", "0.1")
		+ XMLSchemaAttribute::attribute_w_default("trans_mag", xsct_real, "XRW TO DO", "0.4");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string RigidBodyPerturbNoCenterMoverCreator::keyname() const {
	return RigidBodyPerturbNoCenterMover::mover_name();
}

protocols::moves::MoverOP
RigidBodyPerturbNoCenterMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RigidBodyPerturbNoCenterMover );
}

void RigidBodyPerturbNoCenterMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RigidBodyPerturbNoCenterMover::provide_xml_schema( xsd );
}



RigidBodyRandomizeMover::RigidBodyRandomizeMover() :
	parent(),
	partner_( partner_downstream ),
	phi_angle_(360),
	psi_angle_(360),
	update_center_after_move_(true)
{
	moves::Mover::type( "RigidBodyRandomize" );
}

// constructor with arguments
RigidBodyRandomizeMover::RigidBodyRandomizeMover(
	core::pose::Pose const & pose_in,
	int const rb_jump_in,
	Partner const partner_in,
	int phi_angle,
	int psi_angle,
	bool update_center_after_move
):
	RigidBodyMover( rb_jump_in ),
	partner_( partner_in ),
	phi_angle_(phi_angle),
	psi_angle_(psi_angle),
	update_center_after_move_(update_center_after_move)
{
	moves::Mover::type( "RigidBodyRandomize" );
	core::Vector upstream_dummy, downstream_dummy;
	protocols::geometry::centroids_by_jump(pose_in, rb_jump_in, upstream_dummy, downstream_dummy );
	if ( partner_in == partner_downstream ) rot_center_ = downstream_dummy;
	else rot_center_ = upstream_dummy;
}

RigidBodyRandomizeMover::RigidBodyRandomizeMover( RigidBodyRandomizeMover const & src ) :
	//utility::pointer::ReferenceCount(),
	parent( src ),
	partner_( src.partner_ ),
	phi_angle_( src.phi_angle_ ),
	psi_angle_( src.psi_angle_ )
{}

RigidBodyRandomizeMover::~RigidBodyRandomizeMover() = default;

void
RigidBodyRandomizeMover::apply( core::pose::Pose & pose )
{
	core::kinematics::Jump flexible_jump = pose.jump( rb_jump_ );
	TRBM << "Randomize: " << "Jump (before): " << flexible_jump << std::endl;
	core::kinematics::Stub upstream_stub = pose.conformation().upstream_jump_stub( rb_jump_ );
	core::kinematics::Stub downstream_stub = pose.conformation().downstream_jump_stub( rb_jump_ );
	TRBM << "Randomize: " << "Rot. Center (before): "
		<< rot_center_.x() << " "
		<< rot_center_.y() << " "
		<< rot_center_.z() << std::endl;
	// comments for set_rb_center() explain which stub to use when!
	flexible_jump.set_rb_center( dir_, downstream_stub, rot_center_ );
	if ( !freeze_ ) rotation_matrix_ = protocols::geometry::random_reorientation_matrix(phi_angle_, psi_angle_);
	flexible_jump.rotation_by_matrix( upstream_stub, rot_center_,  rotation_matrix_);
	TRBM << "Randomize: " << "Jump (after):  " << flexible_jump << std::endl;
	pose.set_jump( rb_jump_, flexible_jump );

	if ( update_center_after_move_ ) { // update rot_center_ // TODO fix this so we don't update center, because that ruins our freezing ability
		core::Vector dummy_up, dummy_down;
		protocols::geometry::centroids_by_jump(pose, rb_jump_, dummy_up, dummy_down);
		rot_center_ = ( partner_ == 2 ) ? dummy_down : dummy_up;
	}

	TRBM << "Randomize: " << "Rot. Center (after):  "
		<< rot_center_.x() << " "
		<< rot_center_.y() << " "
		<< rot_center_.z() << std::endl;
	TRBM << "Randomize: "  << "---" << std::endl;
}

std::string
RigidBodyRandomizeMover::get_name() const {
	return "RigidBodyRandomizeMover";
}

void
RigidBodyRandomizeMover::show(std::ostream & output) const
{
	RigidBodyMover::show(output);
	output << "\nPhi angle:   " << get_phi() <<
		"\nPsi angle:   " << get_psi() << std::endl;
}

core::Size
RigidBodyRandomizeMover::get_phi() const {
	return phi_angle_;
}

core::Size
RigidBodyRandomizeMover::get_psi() const {
	return psi_angle_;
}

std::ostream
&operator<< ( std::ostream &os, RigidBodyRandomizeMover const &randommover )
{
	randommover.show(os);
	return os;
}


RigidBodySpinMover::RigidBodySpinMover() : parent(), spin_axis_( 0.0 ), default_spin_mag_( true )
{
	moves::Mover::type( "RigidBodySpin" );
}

/// @brief constructor with arguments
///       spin axis is initialized to 0 and then calculated during apply()
RigidBodySpinMover::RigidBodySpinMover(
	int const rb_jump_in
):
	RigidBodyMover( rb_jump_in ), spin_axis_( 0.0 ), update_spin_axis_( true ), default_spin_mag_( true )
{
	moves::Mover::type( "RigidBodySpin" );
}

RigidBodySpinMover::RigidBodySpinMover( RigidBodySpinMover const & /*src*/ ) = default;

RigidBodySpinMover::~RigidBodySpinMover() = default;

void
RigidBodySpinMover::spin_axis ( core::Vector spin_axis_in )
{
	spin_axis_ = spin_axis_in;
	update_spin_axis_ = false;
}

void
RigidBodySpinMover::rot_center ( core::Vector const & rot_center_in )
{
	rot_center_ = rot_center_in;
	update_spin_axis_ = false;
}

void
RigidBodySpinMover::spin_mag( core::Real const & spin_mag )
{
	spin_mag_ = spin_mag;
	default_spin_mag_ = false;
}

void
RigidBodySpinMover::apply( core::pose::Pose & pose )
{
	core::kinematics::Jump flexible_jump = pose.jump( rb_jump_ );
	TRBM << "Spin: " << "Jump (before): " << flexible_jump << std::endl;
	core::kinematics::Stub upstream_stub = pose.conformation().upstream_jump_stub( rb_jump_ );
	core::kinematics::Stub downstream_stub = pose.conformation().downstream_jump_stub( rb_jump_ );

	core::Vector dummy_up, dummy_down;
	if ( update_spin_axis_ ) {
		protocols::geometry::centroids_by_jump(pose, rb_jump_, dummy_up, dummy_down);
		rot_center_ = dummy_down;
		spin_axis_ = dummy_up - rot_center_;
	}

	TRBM << "Spin: " << "Rot (before: "
		<< rot_center_.x() << " "
		<< rot_center_.y() << " "
		<< rot_center_.z() << std::endl;
	// comments for set_rb_center() explain which stub to use when!
	flexible_jump.set_rb_center( dir_, downstream_stub, rot_center_ );
	core::Real spin_angle;
	if ( default_spin_mag_ ) {
		spin_angle = 360.0f*numeric::random::rg().uniform();
	} else {
		spin_angle = spin_mag_ * numeric::random::rg().gaussian();
	}
	flexible_jump.rotation_by_axis( upstream_stub, spin_axis_, rot_center_, spin_angle );
	//flexible_jump.rotation_by_axis( upstream_stub, spin_axis_, rot_center_, 360.0f*numeric::random::rg().uniform() );
	TRBM << "Spin: " << "Jump (after):  " << flexible_jump << std::endl;
	pose.set_jump( rb_jump_, flexible_jump );
	protocols::geometry::centroids_by_jump(pose, rb_jump_, dummy_up, dummy_down);
	rot_center_ = dummy_down;
	TRBM << "Spin: " << "Rot  (after): "
		<< rot_center_.x() << " "
		<< rot_center_.y() << " "
		<< rot_center_.z() << std::endl;
	TRBM << "Spin: " << "---" << std::endl;
}

std::string
RigidBodySpinMover::get_name() const {
	return "RigidBodySpinMover";
}

////////////////////////////////////////////////////////////////////////////////

RigidBodyDeterministicSpinMover::RigidBodyDeterministicSpinMover() : parent()
{
	moves::Mover::type( "RigidBodyDeterministicSpin" );
	angle_magnitude_ = 0.0;

}

/// @brief constructor with arguments
///       takes a complete set of arguments needed for apply
RigidBodyDeterministicSpinMover::RigidBodyDeterministicSpinMover( int const rb_jump_in, core::Vector spin_axis, core::Vector rot_center, float angle_magnitude ):
	parent( rb_jump_in )
{
	moves::Mover::type( "RigidBodyDeterministicSpin" );
	spin_axis_ = spin_axis;
	rot_center_ = rot_center;
	angle_magnitude_ = angle_magnitude;
	update_spin_axis_ = false;
}

RigidBodyDeterministicSpinMover::RigidBodyDeterministicSpinMover( RigidBodyDeterministicSpinMover const & ) = default;

RigidBodyDeterministicSpinMover::~RigidBodyDeterministicSpinMover() = default;

void
RigidBodyDeterministicSpinMover::angle_magnitude( float angle_magnitude )
{
	angle_magnitude_ = angle_magnitude;
}

void
RigidBodyDeterministicSpinMover::apply( core::pose::Pose & pose )
{
	core::kinematics::Jump flexible_jump = pose.jump( rb_jump_ );
	TRBM << "Spin: " << "Jump (before): " << flexible_jump << std::endl;
	core::kinematics::Stub upstream_stub = pose.conformation().upstream_jump_stub( rb_jump_ );
	core::kinematics::Stub downstream_stub = pose.conformation().downstream_jump_stub( rb_jump_ );

	core::Vector dummy_up, dummy_down;
	if ( update_spin_axis_ ) {
		protocols::geometry::centroids_by_jump(pose, rb_jump_, dummy_up, dummy_down);
		rot_center_ = dummy_down;
		spin_axis_ = dummy_up - rot_center_;
	}

	TRBM << "Spin: " << "Rot (before: "
		<< rot_center_.x() << " "
		<< rot_center_.y() << " "
		<< rot_center_.z() << std::endl;
	// comments for set_rb_center() explain which stub to use when!
	flexible_jump.set_rb_center( dir_, downstream_stub, rot_center_ );
	flexible_jump.rotation_by_axis( upstream_stub, spin_axis_, rot_center_, angle_magnitude_ );
	TRBM << "Spin: " << "Jump (after):  " << flexible_jump << std::endl;
	pose.set_jump( rb_jump_, flexible_jump );
	protocols::geometry::centroids_by_jump(pose, rb_jump_, dummy_up, dummy_down);
	rot_center_ = dummy_down;
	TRBM << "Spin: " << "Rot  (after): "
		<< rot_center_.x() << " "
		<< rot_center_.y() << " "
		<< rot_center_.z() << std::endl;
	TRBM << "Spin: " << "---" << std::endl;
}

std::string
RigidBodyDeterministicSpinMover::get_name() const {
	return "RigidBodyDeterministicSpinMover";
}

////////////////////////////////////////////////////////////////////////////////

RigidBodyTiltMover::RigidBodyTiltMover():
	RigidBodyMover(),
	tilt1_mag_( 0.0 ),
	tilt2_mag_( 0.0 ),
	tilt1_center_( 0 ),
	tilt2_center_( 0 ),
	spin_axis_( 0.0 )
{
	moves::Mover::type( "RigidBodyTilt" );
}

/// @brief constructor with arguments
///       spin axis is initialized to 0 and then calculated during apply()
RigidBodyTiltMover::RigidBodyTiltMover(
	int const rb_jump_in,
	core::Real const tilt1_mag_in, //max tilt of partner1 in degrees
	core::Real const tilt2_mag_in, //max tilt of partner2 in degrees
	core::Size const tilt1_center_in, //resID for residue tilt of partner1 is centered at
	core::Size const tilt2_center_in //resID for residue tilt of partner2 is centered at
):
	RigidBodyMover( rb_jump_in ),
	tilt1_mag_(tilt1_mag_in),
	tilt2_mag_(tilt2_mag_in),
	tilt1_center_(tilt1_center_in),
	tilt2_center_(tilt2_center_in),
	spin_axis_( 0.0 )
{
	moves::Mover::type( "RigidBodyTilt" );
}

RigidBodyTiltMover::RigidBodyTiltMover( RigidBodyTiltMover const & ) = default;

RigidBodyTiltMover::~RigidBodyTiltMover() = default;

void
RigidBodyTiltMover::spin_axis ( core::Vector const & spin_axis_in )
{
	spin_axis_ = spin_axis_in;
}

void
RigidBodyTiltMover::apply( core::pose::Pose & pose )
{
	// retrieve stubs based on moveable jump
	core::kinematics::Stub upstream_stub = pose.conformation().upstream_jump_stub( rb_jump_ );
	core::kinematics::Stub downstream_stub = pose.conformation().downstream_jump_stub( rb_jump_ );

	// Find partner centers.
	core::Vector partner1_center, partner2_center;
	protocols::geometry::centroids_by_jump(pose, rb_jump_, partner1_center, partner2_center);

	// update spin axis, if needed
	core::Vector spin_axis( spin_axis_ );
	if ( spin_axis.is_zero() ) {
		TRBM <<"Updating spin axis"<<std::endl;
		spin_axis = partner1_center - partner2_center;
	}

	// tilt partner 1
	if ( tilt1_mag_ > 0 ) {
		core::Vector tilt_center( find_tilt_center(pose, tilt1_center_, partner1_center ) );
		tilt( pose, "1", tilt1_mag_, tilt_center, spin_axis, partner1_center, c2n, upstream_stub, downstream_stub);
	}

	// tilt partner 2
	if ( tilt2_mag_ > 0 ) {
		core::Vector tilt_center( find_tilt_center(pose, tilt2_center_, partner2_center ) );
		tilt( pose, "2", tilt2_mag_, tilt_center, spin_axis, partner2_center, dir_, upstream_stub, downstream_stub);
	}
}

void RigidBodyTiltMover::tilt(
	core::pose::Pose & pose,
	std::string const & which,
	core::Real tilt_mag,
	core::Vector const & tilt_center,
	core::Vector const & spin_axis,
	core::Vector const & partner_center,
	Direction dir,
	core::kinematics::Stub const & upstream_stub,
	core::kinematics::Stub const & downstream_stub
) const {
	// calculate a vector perpendicular to the spin axis
	core::Vector rot_vector1 = (core::Vector(-(spin_axis.y()+spin_axis.z())/spin_axis.x(),1,1));
	// calculate third vector that is perpendicular to spin axis and rot_vector1
	core::Vector rot_vector2 = spin_axis.cross(rot_vector1);

	core::Real tilt_angle1 = tilt_mag-2*tilt_mag*numeric::random::rg().uniform();
	core::Real tilt_angle2 = tilt_mag-2*tilt_mag*numeric::random::rg().uniform();
	TRBM <<"Tilt: tilt of partner "<<which<<" - " <<tilt_angle1 <<" degrees (axis 1) and "<<tilt_angle2<< " degrees (axis 2) (limit: "<<tilt_mag<< " degrees)" <<std::endl;

	// Get jump and apply rotation
	core::kinematics::Jump flexible_jump = pose.jump( rb_jump_ );
	flexible_jump.set_rb_center( dir, downstream_stub, partner_center );
	flexible_jump.rotation_by_axis( upstream_stub, rot_vector1, tilt_center, tilt_angle1);
	flexible_jump.rotation_by_axis( upstream_stub, rot_vector2, tilt_center, tilt_angle2);
	pose.set_jump( rb_jump_, flexible_jump );
}

core::Vector RigidBodyTiltMover::find_tilt_center(
	core::pose::Pose const & pose,
	core::Size tilt_center_res,
	core::Vector const & partner_center
) const {
	if ( tilt_center_res == 0 ) {
		TRBM.Debug << "Tilt center on center of mass." << std::endl;
		return partner_center;
	} else {
		TRBM.Debug << "Tilt center on residue " << tilt_center_res << std::endl;
		return pose.residue(tilt_center_res).atom("CA").xyz();
	}
}


// XRW TEMP std::string
// XRW TEMP RigidBodyTiltMover::get_name() const {
// XRW TEMP  return "RigidBodyTiltMover";
// XRW TEMP }

void
RigidBodyTiltMover::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose )
{
	rb_jump( tag->getOption< int >( "jump", 1 ) );
	core::Vector axis( tag->getOption< core::Real >( "x", 0.0 ), tag->getOption< core::Real >( "y", 0.0 ), tag->getOption< core::Real >( "z", 0.0 ));
	spin_axis( axis );
	tilt1_mag( tag->getOption< core::Real >( "tilt1_mag", 0.0 )  );
	tilt2_mag( tag->getOption< core::Real >( "tilt2_mag", 0.0 )  );
	if ( tag->hasOption( "tilt1_center" ) ) {
		core::Size tilt1_cent( core::pose::parse_resnum( tag->getOption< std::string >( "tilt1_center" ), pose ) );
		TRBM.Debug << "Tilt 1 center " << tilt1_cent << std::endl;
		if ( tilt1_cent == 0 ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "In RigidBodyTiltMover: can't understand value passed to tilt1_center.");
		}
		tilt1_center( tilt1_cent );
	}
	if ( tag->hasOption( "tilt2_center" ) ) {
		core::Size tilt2_cent( core::pose::parse_resnum( tag->getOption< std::string >( "tilt2_center" ), pose ) );
		TRBM.Debug << "Tilt 2 center " << tilt2_cent << std::endl;
		if ( tilt2_cent == 0 ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "In RigidBodyTiltMover: can't understand value passed to tilt2_center.");
		}
		tilt2_center( tilt2_cent );
	}
}

moves::MoverOP
RigidBodyTiltMover::clone() const {
	return moves::MoverOP( new RigidBodyTiltMover(*this) );
}

moves::MoverOP
RigidBodyTiltMover::fresh_instance() const {
	return moves::MoverOP( new RigidBodyTiltMover );
}

// XRW TEMP std::string
// XRW TEMP RigidBodyTiltMoverCreator::keyname() const {
// XRW TEMP  return RigidBodyTiltMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RigidBodyTiltMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new RigidBodyTiltMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RigidBodyTiltMover::mover_name() {
// XRW TEMP  return "RigidBodyTiltMover";
// XRW TEMP }

std::string RigidBodyTiltMover::get_name() const {
	return mover_name();
}

std::string RigidBodyTiltMover::mover_name() {
	return "RigidBodyTiltMover";
}

void RigidBodyTiltMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// TO DO!
	using namespace utility::tag;
	AttributeList attlist; // TO DO: add attributes to this list
	attlist + XMLSchemaAttribute::attribute_w_default("jump", xs_integer, "XRW TO DO", "1")
		+ XMLSchemaAttribute::attribute_w_default("x", xsct_real, "XRW TO DO", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("y", xsct_real, "XRW TO DO", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("z", xsct_real, "XRW TO DO", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("tilt1_mag", xsct_real, "XRW TO DO", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("tilt2_mag", xsct_real, "XRW TO DO", "0.0")
		+ XMLSchemaAttribute("tilt1_center", xs_string, "XRW TO DO")
		+ XMLSchemaAttribute("tilt2_center", xs_string, "XRW TO DO");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string RigidBodyTiltMoverCreator::keyname() const {
	return RigidBodyTiltMover::mover_name();
}

protocols::moves::MoverOP
RigidBodyTiltMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RigidBodyTiltMover );
}

void RigidBodyTiltMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RigidBodyTiltMover::provide_xml_schema( xsd );
}


////////////////////////////////////////////////////////////////////////////////

RigidBodyTransMover::RigidBodyTransMover( bool vary_stepsize ) : RigidBodyMover()
{
	moves::Mover::type( "RigidBodyTrans" );
	vary_stepsize_ = vary_stepsize;
}

// constructor with arguments
RigidBodyTransMover::RigidBodyTransMover(
	core::pose::Pose const & pose_in,
	int const rb_jump_in,
	bool vary_stepsize
):
	RigidBodyMover( rb_jump_in )
{
	moves::Mover::type( "RigidBodyTrans" );
	step_size_ = 1.0;
	trans_axis_ = centroid_axis(pose_in);
	vary_stepsize_ = vary_stepsize;
}

RigidBodyTransMover::RigidBodyTransMover( core::Vector const & trans_axis, int const rb_jump_in, bool vary_stepsize  ) :
	RigidBodyMover( rb_jump_in ), trans_axis_(trans_axis)
{
	moves::Mover::type( "RigidBodyTrans" );
	step_size_ = 1.0;
	vary_stepsize_ = vary_stepsize;
}

RigidBodyTransMover::RigidBodyTransMover( RigidBodyTransMover const & /*src*/ ) = default;

RigidBodyTransMover::~RigidBodyTransMover() = default;

core::Vector
RigidBodyTransMover::centroid_axis(core::pose::Pose const & pose_in) const
{
	core::Vector upstream_dummy, downstream_dummy;
	protocols::geometry::centroids_by_jump(pose_in, rb_jump(), upstream_dummy, downstream_dummy );
	return downstream_dummy - upstream_dummy;
}

void
RigidBodyTransMover::apply( core::pose::Pose & pose )
{
	// get translation axis
	core::Vector axis( trans_axis_ );
	if ( axis.is_zero() ) {
		axis = centroid_axis(pose);
	}

	// set flexible jump in the pose
	core::kinematics::Jump flexible_jump = pose.jump( rb_jump_ );
	TRBM << "Translate: " << "Jump (before): " << flexible_jump << std::endl;

	// get upstream stub
	core::kinematics::Stub upstream_stub = pose.conformation().upstream_jump_stub( rb_jump_ );

	// variable stepsize means that the steps are larger when the partners are
	// far away and that they get smaller the closer the partners are;
	// should save time
	if ( vary_stepsize_ == true ) {
		core::Vector centroid_vector( centroid_axis( pose ) );
		core::Real distance( centroid_vector.length() );

		// if partners are close, take small stepsize of 1A
		if ( distance <= 5 ) {
			step_size_ = 1.0;
		} else {
			step_size_ = std::sqrt( distance );
		}
		TRBM << "distance between partners: " << distance << std::endl;
	}

	// do the actual translation
	flexible_jump.translation_along_axis( upstream_stub, axis, step_size_ );
	TRBM << "Translate: " << "Jump (after):  " << flexible_jump << std::endl;

	// set flexible jump in the pose
	pose.set_jump( rb_jump_, flexible_jump );

}// apply

/////////////////// ROSETTA SCRIPTS STUFF //////////////////////////////

// XRW TEMP std::string
// XRW TEMP RigidBodyTransMover::get_name() const {
// XRW TEMP  return "RigidBodyTransMover";
// XRW TEMP }

void
RigidBodyTransMover::parse_my_tag( utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	step_size( tag->getOption< core::Real >( "distance", 1.0 ) );
	rb_jump( tag->getOption< int >( "jump", 1 ) );
	core::Vector axis( tag->getOption< core::Real >( "x", 0.0 ), tag->getOption< core::Real >( "y", 0.0 ), tag->getOption< core::Real >( "z", 0.0 ));
	trans_axis( axis );
	vary_stepsize_ = tag->getOption< bool >( "vary_stepsize", false );
}

moves::MoverOP
RigidBodyTransMover::clone() const {
	return moves::MoverOP( new RigidBodyTransMover(*this) );
}

moves::MoverOP
RigidBodyTransMover::fresh_instance() const {
	return moves::MoverOP( new RigidBodyTransMover );
}

// XRW TEMP std::string
// XRW TEMP RigidBodyTransMoverCreator::keyname() const {
// XRW TEMP  return RigidBodyTransMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP RigidBodyTransMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new RigidBodyTransMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP RigidBodyTransMover::mover_name() {
// XRW TEMP  return "RigidBodyTransMover";
// XRW TEMP }

std::string RigidBodyTransMover::get_name() const {
	return mover_name();
}

std::string RigidBodyTransMover::mover_name() {
	return "RigidBodyTransMover";
}

void RigidBodyTransMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	// TO DO!
	using namespace utility::tag;
	AttributeList attlist; // TO DO: add attributes to this list
	attlist + XMLSchemaAttribute::attribute_w_default("distance", xsct_real, "XRW TO DO", "1.0")
		+ XMLSchemaAttribute::attribute_w_default("jump", xs_integer, "XRW TO DO", "1")
		+ XMLSchemaAttribute::attribute_w_default("x", xsct_real, "XRW TO DO", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("y", xsct_real, "XRW TO DO", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("z", xsct_real, "XRW TO DO", "0.0")
		+ XMLSchemaAttribute::attribute_w_default("vary_stepsize", xsct_rosetta_bool, "XRW TO DO", "false");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}

std::string RigidBodyTransMoverCreator::keyname() const {
	return RigidBodyTransMover::mover_name();
}

protocols::moves::MoverOP
RigidBodyTransMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RigidBodyTransMover );
}

void RigidBodyTransMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RigidBodyTransMover::provide_xml_schema( xsd );
}


////////////////////////////////////////////////////////////////////////////////

UniformSphereTransMover::UniformSphereTransMover() : parent(), step_size_(1), random_step_(0), trans_axis_()
{
	moves::Mover::type( "UniformSphereTrans" );
	reset_trans_axis();
}

// constructor with arguments
UniformSphereTransMover::UniformSphereTransMover(
	int const rb_jump_in,
	core::Real step_size_in
):
	parent( rb_jump_in ),
	step_size_( step_size_in ),
	random_step_(0),
	trans_axis_()

{
	moves::Mover::type( "UniformSphereTrans" );
	reset_trans_axis(); // start with a random trans_axis, freeze is valid without first calling apply
}

UniformSphereTransMover::UniformSphereTransMover( UniformSphereTransMover const & /*src*/ ) = default;

UniformSphereTransMover::~UniformSphereTransMover() = default;

void UniformSphereTransMover::reset_trans_axis(){
	do {
		trans_axis_.assign(
			step_size_*2*(numeric::random::rg().uniform()-0.5),
			step_size_*2*(numeric::random::rg().uniform()-0.5),
			step_size_*2*(numeric::random::rg().uniform()-0.5) );
		random_step_ = trans_axis_.length();
	} while( random_step_ > step_size_ );
	trans_axis_.normalize();
}

/// @details Sample points in a cube randomly, and discard ones that are outside the sphere.
/// This gives us *uniform* sampling of the space inside, whereas
/// choosing a random distance and a random direction samples more near the center.
void UniformSphereTransMover::apply( core::pose::Pose & pose )
{
	if ( ! freeze_ ) reset_trans_axis();

	RigidBodyTransMover mover( pose, rb_jump_);
	mover.trans_axis( trans_axis_ );
	mover.step_size( random_step_ );
	mover.apply( pose );
}

std::string
UniformSphereTransMover::get_name() const {
	return "UniformSphereTransMover";
}

RigidBodyDofRandomizeMover::RigidBodyDofRandomizeMover() :
	RigidBodyMover()
{
	moves::Mover::type( "RigidBodyDofRandomize" );
}

// @details rigid body randomization according to SymDof information. It randomizes all
// allowed dof dor a single jump
RigidBodyDofRandomizeMover::RigidBodyDofRandomizeMover(
	int const rb_jump_in,
	core::conformation::symmetry::SymDof dof
):
	RigidBodyMover( rb_jump_in )
{
	moves::Mover::type( "RigidBodyDofRandomize" );
	rb_jump_ = rb_jump_in;
	dof_ = dof;
}


RigidBodyDofRandomizeMover::RigidBodyDofRandomizeMover( RigidBodyDofRandomizeMover const & /*src*/ ) = default;

RigidBodyDofRandomizeMover::~RigidBodyDofRandomizeMover() = default;

void RigidBodyDofRandomizeMover::apply( core::pose::Pose & pose )
{
	using namespace core::conformation::symmetry;

	core::kinematics::Jump flexible_jump = pose.jump( rb_jump_ );
	TRBM.Debug << "Randomize: " << "Jump (before): " << flexible_jump << std::endl;

	// randomize the translational dofs
	Vector trans = flexible_jump.get_translation();
	for ( int i = X_DOF; i <= Z_DOF; ++i ) {
		// randomize if dof is allowed and we have specified a range to randomize the translation
		// for the specific dof
		if ( dof_.allow_dof(i) && ( dof_.has_range1(i) || dof_.has_range1_lower(i) ) ) {
			core::Real new_trans(0);
			if ( dof_.has_range1(i) ) {
				new_trans = numeric::random::rg().uniform()*(dof_.range1_upper(i) - dof_.range1_lower(i) ) + dof_.range1_lower(i);
			} else {
				new_trans = dof_.range1_lower(i);
			}
			Real scale = (dof_.jump_direction(i) == c2n) ? -1 : 1;
			if ( i == X_DOF ) trans[0] = scale*new_trans;
			if ( i == Y_DOF ) trans[1] = scale*new_trans;
			if ( i == Z_DOF ) trans[2] = scale*new_trans;
		}
	}
	flexible_jump.set_translation( trans );
	pose.set_jump( rb_jump_, flexible_jump );

	// Now apply rotations
	for ( int i = X_ANGLE_DOF; i <= Z_ANGLE_DOF; ++i ) {
		// If the user has set 360 degrees rotation for X_ANGLE, Y_ANGLE and Z_ANGLE the n we want
		// uniform randomization as well. Observe that randomizing x,y,z independent does not give
		// uniform randomization!
		if ( std::abs( dof_.range1_lower(X_ANGLE_DOF) - dof_.range1_upper(X_ANGLE_DOF) ) == 360 &&
				std::abs( dof_.range1_lower(Y_ANGLE_DOF) - dof_.range1_upper(Y_ANGLE_DOF) ) == 360 &&
				std::abs( dof_.range1_lower(Z_ANGLE_DOF) - dof_.range1_upper(Z_ANGLE_DOF) ) == 360 ) {
			if ( dof_.jump_direction(i) == c2n ) flexible_jump.reverse();
			numeric::xyzMatrix< Real > rot = protocols::geometry::random_reorientation_matrix()*
				flexible_jump.get_rotation();
			flexible_jump.set_rotation( rot );
			if ( dof_.jump_direction(i) == c2n ) flexible_jump.reverse();
			pose.set_jump( rb_jump_, flexible_jump );
			break;
		}
		// randomize each rotational dof independently
		if ( dof_.allow_dof(i) && dof_.has_range1(i) ) {
			numeric::xyzMatrix< Real > rot;
			core::Real angle = numeric::random::rg().uniform()*(dof_.range1_upper(i) - dof_.range1_lower(i) ) + dof_.range1_lower(i);

			if ( i == X_ANGLE_DOF ) rot = numeric::x_rotation_matrix_degrees(angle);
			if ( i == Y_ANGLE_DOF ) rot = numeric::y_rotation_matrix_degrees(angle);
			if ( i == Z_ANGLE_DOF ) rot = numeric::z_rotation_matrix_degrees(angle);
			if ( dof_.jump_direction(i) == c2n ) flexible_jump.reverse();
			rot *= flexible_jump.get_rotation();
			flexible_jump.set_rotation( rot );
			if ( dof_.jump_direction(i) == c2n ) flexible_jump.reverse();
			pose.set_jump( rb_jump_, flexible_jump );
		}
		if ( dof_.allow_dof(i) && !dof_.has_range1(i) && dof_.has_range1_lower(i) ) {
			numeric::xyzMatrix< Real > rot;

			if ( i == X_ANGLE_DOF ) rot = numeric::x_rotation_matrix_degrees( dof_.range1_lower(i) );
			if ( i == Y_ANGLE_DOF ) rot = numeric::y_rotation_matrix_degrees( dof_.range1_lower(i) );
			if ( i == Z_ANGLE_DOF ) rot = numeric::z_rotation_matrix_degrees( dof_.range1_lower(i) );
			if ( dof_.jump_direction(i) == c2n ) flexible_jump.reverse();
			rot *= flexible_jump.get_rotation();
			flexible_jump.set_rotation( rot );
			if ( dof_.jump_direction(i) == c2n ) flexible_jump.reverse();
			pose.set_jump( rb_jump_, flexible_jump );
		}
	}
	TRBM.Debug << "Randomize: " << "Jump (after):  " << flexible_jump << std::endl;
	TRBM.Debug << "Randomize: "  << "---" << std::endl;
}

std::string
RigidBodyDofRandomizeMover::get_name() const {
	return "RigidBodyDofRandomizeMover";
}

RigidBodyDofSeqRandomizeMover::RigidBodyDofSeqRandomizeMover() :
	parent()
{
	moves::Mover::type( "RigidBodyDofRandomize" );
}

// constructor with arguments
RigidBodyDofSeqRandomizeMover::RigidBodyDofSeqRandomizeMover(
	std::map< Size, core::conformation::symmetry::SymDof > const & dofs
):
	parent()
{
	moves::Mover::type( "RigidBodyDofSeqRandomize" );
	dofs_ = dofs;
}

RigidBodyDofSeqRandomizeMover::RigidBodyDofSeqRandomizeMover( RigidBodyDofSeqRandomizeMover const & /*src*/ ) = default;

RigidBodyDofSeqRandomizeMover::~RigidBodyDofSeqRandomizeMover() = default;


// @details go through and perturb all movable dofs in sequence
void RigidBodyDofSeqRandomizeMover::apply( core::pose::Pose & pose )
{
	using namespace core::conformation::symmetry;

	std::map< Size, SymDof >::iterator it;
	auto it_begin = dofs_.begin();
	auto it_end = dofs_.end();
	for ( it = it_begin; it != it_end; ++it ) {
		int jump_nbr ( (*it).first );
		SymDof dof ( (*it).second );
		RigidBodyDofRandomizeMover dofrandommover( jump_nbr, dof );
		dofrandommover.apply( pose );
	}
}

std::string
RigidBodyDofSeqRandomizeMover::get_name() const {
	return "RigidBodyDofSeqRandomizeMover";
}


// default constructor
RigidBodyDofTransMover::RigidBodyDofTransMover() : parent(),
	last_slide_good_(false),
	jump_dir_(n2c)
{
	moves::Mover::type( "RigidBodyDofTrans" );
}


/// @details Constructor for a rigid body translation mover
/// moves only along directions defined by a vector
/// of dofs (x,y or z). If more than two directions are
/// allowed the move along them as well. This probably
/// never makes sense. Perhaps more logical to select one
/// direction randomly then?
RigidBodyDofTransMover::RigidBodyDofTransMover(
	core::conformation::symmetry::SymDof dof,
	int const rb_jump_in,
	core::Real step_size
):
	RigidBodyMover( rb_jump_in ),
	last_slide_good_(false)
{
	dof_ = dof;
	moves::Mover::type( "RigidBodyDofTrans" );
	jump_dir_ = n2c;
	// This is fishy. We should not have different directions for the same jump
	// need to put in checks for that...
	if ( dof.jump_direction(1) == c2n || dof.jump_direction(2) == c2n
			|| dof.jump_direction(3) == c2n ) jump_dir_ = c2n;
	step_size_ = step_size;
	core::Vector zero(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
	trans_axis_ = zero;
	if ( dof.allow_dof(1) ) {
		zero += x;
	}
	if ( dof.allow_dof(2) ) {
		zero += y;
	}
	if ( dof.allow_dof(3) ) {
		zero += z;
	}
	trans_axis_ = zero;
}

// constructor with arguments
RigidBodyDofTransMover::RigidBodyDofTransMover(
	std::map< Size, core::conformation::symmetry::SymDof > dofs
):
	RigidBodyMover(),
	last_slide_good_(false)
{

	utility::vector1< int > trans_jumps;

	moves::Mover::type( "RigidBodyDofTrans" );
	jump_dir_ = n2c;
	step_size_ = 0.5;

	// Save jumps that are allowed to move and have a translation dof
	std::map< Size, core::conformation::symmetry::SymDof >::iterator it;
	auto it_begin = dofs.begin();
	auto it_end = dofs.end();
	for ( it = it_begin; it != it_end; ++it ) {
		core::conformation::symmetry::SymDof dof ( (*it).second );
		if ( dof.allow_dof(1) || dof.allow_dof(2) || dof.allow_dof(3) ) {
			int jump_nbr ( (*it).first );
			trans_jumps.push_back( jump_nbr );
		}
	}

	if ( trans_jumps.empty() ) {
		TRBM.Warning << "no movable jumps!" << std::endl;
		return;
	}
	rb_jump_ = numeric::random::rg().random_element( trans_jumps );
	auto jump_iterator =
		dofs.find( rb_jump_ );
	if ( jump_iterator == dofs.end() ) {
		TRBM.Warning << "jump dof not found!" << std::endl;
	} else {
		dof_ = (*jump_iterator).second ;
		// This is fishy. We should not have different directions for the same jump
		// need to put in checks for that...
		if ( dof_.jump_direction(1) == c2n || dof_.jump_direction(2) == c2n
				|| dof_.jump_direction(3) == c2n ) jump_dir_ = c2n;
		core::Vector zero(0,0,0), x(1,0,0), y(0,1,0), z(0,0,1);
		trans_axis_ = zero;
		if ( dof_.allow_dof(1) ) {
			zero += x;
		}
		if ( dof_.allow_dof(2) ) {
			zero += y;
		}
		if ( dof_.allow_dof(3) ) {
			zero += z;
		}
		trans_axis_ = zero;
	}
}

RigidBodyDofTransMover::RigidBodyDofTransMover( RigidBodyDofTransMover const & src ) : parent( src ),
	last_slide_good_(false),
	jump_dir_( src.jump_dir_ ),
	step_size_( src.step_size_ ),
	trans_axis_( src.trans_axis_ )
{}

RigidBodyDofTransMover::~RigidBodyDofTransMover() = default;


void RigidBodyDofTransMover::apply( core::pose::Pose & pose )
{
	last_slide_good_ = true;
	core::kinematics::Jump flexible_jump = pose.jump( rb_jump_ );
	int c2n(-1);
	TRBM.Debug << "Translate: " << "Jump (before): " << flexible_jump << std::endl;
	Vector trans_start ( flexible_jump.get_translation() );
	if ( jump_dir_ == c2n ) flexible_jump.reverse();

	// if range2_is_bound is set, make sure jump stays within bound
	core::Vector x_i = trans_start + step_size_*trans_axis_;
	for ( int ii=1; ii<=3; ++ii ) {
		if ( dof_.allow_dof(ii) && dof_.range2_is_bound(ii) && dof_.has_range2_lower(ii) ) {
			if ( x_i[ii] < dof_.range2_lower(ii) || x_i[ii] > dof_.range2_upper(ii) ) {
				x_i[ii] = trans_start[ii];
				last_slide_good_ = false;
			}
		}
	}
	flexible_jump.set_translation( x_i );
	if ( jump_dir_ == c2n ) flexible_jump.reverse();
	TRBM.Debug << "Translate: " << "Jump (after):  " << flexible_jump << std::endl;
	pose.set_jump( rb_jump_, flexible_jump );
}

std::string
RigidBodyDofTransMover::get_name() const {
	return "RigidBodyDofTransMover";
}

// default constructor
RigidBodyDofSeqTransMover::RigidBodyDofSeqTransMover() : RigidBodyMover()
{
	moves::Mover::type( "RigidBodyDofSeqTrans" );
}


// @details go through all movable dofs for which translatiions are allowed
// and apply a translation
// constructor with arguments
RigidBodyDofSeqTransMover::RigidBodyDofSeqTransMover(
	std::map< Size, core::conformation::symmetry::SymDof > dofs
):
	RigidBodyMover()
{

	utility::vector1< int > trans_jumps;

	moves::Mover::type( "RigidBodyDofSeqTrans" );
	step_size_ = 0.5;
	trans_axis_ = Vector(1,0,0);
	// Save jumps that are allowed to move and have a translation dof
	std::map< Size, core::conformation::symmetry::SymDof >::iterator it;
	auto it_begin = dofs.begin();
	auto it_end = dofs.end();
	for ( it = it_begin; it != it_end; ++it ) {
		core::conformation::symmetry::SymDof dof ( (*it).second );
		if ( dof.allow_dof(1) || dof.allow_dof(2) || dof.allow_dof(3) ) {
			int jump_nbr ( (*it).first );
			trans_jumps.push_back( jump_nbr );
		}
	}

	if ( trans_jumps.empty() ) {
		TRBM.Warning << "no movable jumps!" << std::endl;
		return;
	}
	dofs_ = dofs;
	rb_jumps_ = trans_jumps;
}


RigidBodyDofSeqTransMover::RigidBodyDofSeqTransMover( RigidBodyDofSeqTransMover const & /*src*/ ) = default;

RigidBodyDofSeqTransMover::~RigidBodyDofSeqTransMover() = default;

void RigidBodyDofSeqTransMover::apply( core::pose::Pose & pose )
{

	std::map< Size, core::conformation::symmetry::SymDof >::iterator jump_iterator;
	utility::vector1< int >::iterator start, end, it;
	start = rb_jumps_.begin();
	end = rb_jumps_.end();

	//random__shuffle(rb_jumps_.begin(), rb_jumps_.end() );
	numeric::random::random_permutation(rb_jumps_.begin(), rb_jumps_.end(), numeric::random::rg());

	for ( it = start; it != end; ++it ) {
		jump_iterator = dofs_.find( *it );
		if ( jump_iterator == dofs_.end() ) {
			TRBM.Warning << "jump dof not found!" << std::endl;
		} else {
			core::conformation::symmetry::SymDof dof( (*jump_iterator).second );
			RigidBodyDofTransMover dofmover( dof, *it, step_size_ );
			// Silly, just reverse the direction if this vector is reversed
			// Since we don't store the direction in this mover
			// trans_axis_ serves as a storage for the direction
			if ( trans_axis_(1) < 0 ) dofmover.trans_axis().negate();
			dofmover.apply( pose );
		}
	}
}


std::string
RigidBodyDofSeqTransMover::get_name() const {
	return "RigidBodyDofSeqTransMover";
}


RigidBodyDofRandomTransMover::RigidBodyDofRandomTransMover() : parent()
{
	moves::Mover::type( "RigidBodyDofRandomTrans" );
}

// @details go randomly set a random translation. Select all
// movable dofs but apply the randomization in random order

// constructor with arguments
RigidBodyDofRandomTransMover::RigidBodyDofRandomTransMover(
	std::map< Size, core::conformation::symmetry::SymDof > dofs
):
	parent()
{

	utility::vector1< int > trans_jumps;

	moves::Mover::type( "RigidBodyDofRandomTrans" );
	step_size_ = 0.5;
	trans_axis_ = Vector(1,0,0);
	// Save jumps that are allowed to move and have a translation dof
	std::map< Size, core::conformation::symmetry::SymDof >::iterator it;
	auto it_begin = dofs.begin();
	auto it_end = dofs.end();
	for ( it = it_begin; it != it_end; ++it ) {
		core::conformation::symmetry::SymDof dof ( (*it).second );
		if ( dof.allow_dof(1) || dof.allow_dof(2) || dof.allow_dof(3) ) {
			int jump_nbr ( (*it).first );
			trans_jumps.push_back( jump_nbr );
		}
	}

	if ( trans_jumps.empty() ) {
		TRBM.Warning << "no movable jumps!" << std::endl;
		return;
	}
	dofs_ = dofs;
	rb_jumps_ = trans_jumps;
}

RigidBodyDofRandomTransMover::RigidBodyDofRandomTransMover( RigidBodyDofRandomTransMover const & /*src*/ ) = default;

RigidBodyDofRandomTransMover::~RigidBodyDofRandomTransMover() = default;


void RigidBodyDofRandomTransMover::apply( core::pose::Pose & pose )
{

	std::map< Size, core::conformation::symmetry::SymDof >::iterator jump_iterator;
	utility::vector1< int >::iterator start, end;
	start = rb_jumps_.begin();
	end = rb_jumps_.end();

	//random__shuffle(rb_jumps_.begin(), rb_jumps_.end() );
	numeric::random::random_permutation(rb_jumps_.begin(), rb_jumps_.end(), numeric::random::rg());

	int jump_;
	if ( rb_jumps_.size() < 1 ) return;
	else jump_ = rb_jumps_[1];
	jump_iterator = dofs_.find( jump_ );
	if ( jump_iterator == dofs_.end() ) {
		TRBM.Warning << "jump dof not found!" << std::endl;
	} else {
		core::conformation::symmetry::SymDof dof( (*jump_iterator).second );
		RigidBodyDofTransMover dofmover( dof, (*jump_iterator).first, step_size_ );
		// Silly, just reverse the direction if this vector is reversed
		// Since we don't store the direction in this mover
		// trans_axis_ serves as a storage for the direction
		if ( trans_axis_(1) < 0 ) dofmover.trans_axis().negate();
		dofmover.apply( pose );
	}
}


std::string
RigidBodyDofRandomTransMover::get_name() const {
	return "RigidBodyDofRandomTransMover";
}

// @details apply a random transformation to a certain jump and use
// the jump step according to dof information if dof_range1 exists.
// Othervise use the trans and rot magnitudes from the constructor
RigidBodyDofPerturbMover::RigidBodyDofPerturbMover(
	int const rb_jump_in,
	core::conformation::symmetry::SymDof dof,
	core::Real const rot_mag_in,
	core::Real const trans_mag_in
):
	RigidBodyMover(),
	rot_mag_( rot_mag_in ),
	trans_mag_( trans_mag_in )
{
	//  TRBM.Debug << "rot_mag " << rot_mag_in << std::endl;
	//  TRBM.Debug << "trans_mag " << trans_mag_in << std::endl;
	moves::Mover::type( "RigidBodyDofPerturbMover" );

	rb_jump_ = rb_jump_in;
	dof_ = dof;
}


// @details apply a random transformation to a certain jump and use
// the jump step according to dof information if dof_range1 exists.
// Othervise use the trans and rot magnitudes from the constructor.
// The jump is selected randomly from the allowed dofs
RigidBodyDofPerturbMover::RigidBodyDofPerturbMover(
	std::map< Size, core::conformation::symmetry::SymDof > dofs,
	Real const rot_mag_in,
	Real const trans_mag_in
):
	RigidBodyMover(),
	rot_mag_( rot_mag_in ),
	trans_mag_( trans_mag_in )
{
	utility::vector1< int > moving_jumps;

	//  TRBM.Debug << "rot_mag " << rot_mag_in << std::endl;
	//  TRBM.Debug << "trans_mag " << trans_mag_in << std::endl;
	moves::Mover::type( "RigidBodyDofPerturbMover" );

	std::map< Size, core::conformation::symmetry::SymDof >::iterator it;
	auto it_begin = dofs.begin();
	auto it_end = dofs.end();
	for ( it = it_begin; it != it_end; ++it ) {
		moving_jumps.push_back( (*it).first );
	}

	if ( moving_jumps.empty() ) {
		TRBM.Warning << "no movable jumps!" << std::endl;
		return;
	}
	rb_jump_ = numeric::random::rg().random_element( moving_jumps );
	auto jump_iterator =
		dofs.find( rb_jump_ );
	if ( jump_iterator == dofs.end() ) {
		TRBM.Warning << "jump dof not found!" << std::endl;
	} else {
		core::conformation::symmetry::SymDof dof( (*jump_iterator).second );
		dof_ = dof;
	}
}

RigidBodyDofPerturbMover::RigidBodyDofPerturbMover( RigidBodyDofPerturbMover const & /*src*/ ) = default;

RigidBodyDofPerturbMover::~RigidBodyDofPerturbMover() = default;

void RigidBodyDofPerturbMover::apply( core::pose::Pose & pose )
{
	core::kinematics::Jump flexible_jump = pose.jump( rb_jump_ );

	int c2n(-1);

	for ( Size i = 1; i<= 3; ++i ) {
		if ( dof_.allow_dof(i) ) {
			// the dat in the dof takes precedence
			core::Real transmag;
			if ( dof_.has_range2_lower(i) && !dof_.range2_is_bound(i) ) transmag = dof_.range2_lower(i);
			else transmag = trans_mag_;
			if ( dof_.jump_direction(i) == c2n ) flexible_jump.reverse();
			flexible_jump.gaussian_move_single_rb( dir_, transmag, i );
			if ( dof_.jump_direction(i) == c2n ) flexible_jump.reverse();

			// if range2_is_bound is set, make sure jump stays within bound
			if ( dof_.has_range2_lower(i) && dof_.range2_is_bound(i) ) {
				core::Vector trans_i = flexible_jump.rt().get_translation();
				trans_i(i) = std::max( dof_.range2_lower(i), trans_i(i) );
				trans_i(i) = std::min( dof_.range2_upper(i), trans_i(i) );
				flexible_jump.set_translation( trans_i );
			}
		}
	}
	for ( Size i = 4; i<= 6; ++i ) {
		if ( dof_.allow_dof(i) ) {
			// the dat in the dof takes precedence
			core::Real rotmag;
			if ( dof_.has_range2_lower(i) && !dof_.range2_is_bound(i) ) rotmag = dof_.range2_lower(i);
			else rotmag = rot_mag_;
			if ( dof_.jump_direction(i) == c2n ) flexible_jump.reverse();
			flexible_jump.gaussian_move_single_rb( dir_, rotmag, i );
			if ( dof_.jump_direction(i) == c2n ) flexible_jump.reverse();
		}
	}
	pose.set_jump( rb_jump_, flexible_jump );
}

std::string
RigidBodyDofPerturbMover::get_name() const {
	return "RigidBodyDofPerturbMover";
}

// @details apply perturbations to all allowed dofs. Apply them in sequential order.
RigidBodyDofSeqPerturbMover::RigidBodyDofSeqPerturbMover(
	std::map< Size, core::conformation::symmetry::SymDof > dofs,
	Real const rot_mag_in,
	Real const trans_mag_in
):
	RigidBodyMover(),
	rot_mag_( rot_mag_in ),
	trans_mag_( trans_mag_in )
{
	utility::vector1< int > moving_jumps;

	//  TRBM.Debug << "rot_mag " << rot_mag_in << std::endl;
	//  TRBM.Debug << "trans_mag " << trans_mag_in << std::endl;
	moves::Mover::type( "RigidBodyDofSeqPerturbMover" );

	std::map< Size, core::conformation::symmetry::SymDof >::iterator it;
	auto it_begin = dofs.begin();
	auto it_end = dofs.end();
	for ( it = it_begin; it != it_end; ++it ) {
		moving_jumps.push_back( (*it).first );
	}

	if ( moving_jumps.empty() ) {
		TRBM.Warning << "no movable jumps!" << std::endl;
		return;
	}
	rb_jumps_ = moving_jumps;
	dofs_ = dofs;
}

RigidBodyDofSeqPerturbMover::RigidBodyDofSeqPerturbMover( RigidBodyDofSeqPerturbMover const & /*src*/ ) = default;

RigidBodyDofSeqPerturbMover::~RigidBodyDofSeqPerturbMover() = default;

void RigidBodyDofSeqPerturbMover::apply( core::pose::Pose & pose )
{
	std::map< Size, core::conformation::symmetry::SymDof >::iterator jump_iterator;
	utility::vector1< int >::iterator start, end, it;
	start = rb_jumps_.begin();
	end = rb_jumps_.end();
	// Shuffle the order by which we visit the jumps
	//random__shuffle(rb_jumps_.begin(), rb_jumps_.end() );
	numeric::random::random_permutation(rb_jumps_.begin(), rb_jumps_.end(), numeric::random::rg());
	// Iterate over all available translation jumps
	// and do a translation for its allowd translation dofs
	for ( it = start; it != end; ++it ) {
		jump_iterator = dofs_.find( *it );
		if ( jump_iterator == dofs_.end() ) {
			TRBM.Warning << "jump dof not found!" << std::endl;
		} else {
			core::conformation::symmetry::SymDof dof( (*jump_iterator).second );
			RigidBodyDofPerturbMover dofmover( *it, dof, rot_mag_, trans_mag_ );
			dofmover.apply( pose );
		}
	}
}

std::string
RigidBodyDofSeqPerturbMover::get_name() const {
	return "RigidBodyDofSeqPerturbMover";
}


}  // namespace rigid
}  // namespace protocols
