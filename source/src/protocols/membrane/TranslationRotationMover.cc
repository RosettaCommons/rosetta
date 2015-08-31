// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief      Does translation and rotation of a pose or the membrane
/// @details Translates along a vector and rotates a pose or membrane around
///    the axis perpendicular to both vectors. Works for both
///    fixed protein/flexible membrane or fixed membrane/flexible protein
///    Takes the jump number as input, default is the membrane jump
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_TranslationRotationMover_cc
#define INCLUDED_protocols_membrane_TranslationRotationMover_cc

// Unit Headers
#include <protocols/membrane/TranslationRotationMover.hh>
#include <protocols/membrane/TranslationRotationMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/LipidAccInfo.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/util.hh>

// Package Headers
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.membrane.TranslationRotationMover" );

namespace protocols {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::moves;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Translation vector can be defined in -mp:setup center
///   flag to translate the new pose to. The mover is a general mover
///   but used mainly on membrane proteins, that's why we use this flag here
///   The default jump is going to be the membrane jump, but you can also
///   specify your own
TranslationMover::TranslationMover() :
	protocols::moves::Mover(),
	translation_vector_(0, 0, 0),
	jumpnum_( 0 )
{
	register_options();
	init_from_cmd();
}

/// @brief Custom Constructor
/// @details User can specify a translation vector
TranslationMover::TranslationMover(
	Vector translation_vector
) :
	protocols::moves::Mover(),
	translation_vector_( translation_vector ),
	jumpnum_(0)
{
	register_options();
	init_from_cmd();
}

/// @brief Custom Constructor
/// @details User can specify a translation vector and jumpnum
TranslationMover::TranslationMover(
	Vector translation_vector,
	Size jumpnum
) :
	protocols::moves::Mover(),
	translation_vector_( translation_vector ),
	jumpnum_( jumpnum )
{
	register_options();
	init_from_cmd();
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
TranslationMover::TranslationMover( TranslationMover const & src ) :
	protocols::moves::Mover( src ),
	translation_vector_( src.translation_vector_ ),
	jumpnum_( src.jumpnum_ )
{}

/// @brief Destructor
TranslationMover::~TranslationMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
TranslationMover::clone() const {
	return ( protocols::moves::MoverOP( new TranslationMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
TranslationMover::fresh_instance() const {
	return protocols::moves::MoverOP( new TranslationMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
TranslationMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	// Read in membrane center
	if ( tag->hasOption( "center" ) ) {
		std::string center = tag->getOption< std::string >( "center" );
		utility::vector1< std::string > str_cen = utility::string_split_multi_delim( center, ":,'`~+*&|;." );

		if ( str_cen.size() != 3 ) {
			utility_exit_with_message( "Cannot read in xyz center vector from string - incorrect length!" );
		} else {
			translation_vector_.x() = std::atof( str_cen[1].c_str() );
			translation_vector_.y() = std::atof( str_cen[2].c_str() );
			translation_vector_.z() = std::atof( str_cen[3].c_str() );
		}
	}
}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
TranslationMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new TranslationMover() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
TranslationMoverCreator::keyname() const {
	return TranslationMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
TranslationMoverCreator::mover_name() {
	return "TranslationMover";
}


/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover
std::string
TranslationMover::get_name() const {
	return "TranslationMover";
}

/// @brief Translate the pose along the defined vector
void
TranslationMover::apply( Pose & pose ) {

	using namespace core::conformation::membrane;
	using namespace protocols::membrane::geometry;

	// checking input
	core::Vector zero( 0, 0, 0 );
	if ( translation_vector_ == zero ) {
		TR << "WARNING: Old and new centers are identical. Skipping rotation!" << std::endl;
	} else {

		TR << "Translating the pose" << std::endl;

		// starting foldtree
		TR << "Starting foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
		pose.fold_tree().show( TR );
		core::kinematics::FoldTree orig_ft = pose.fold_tree();

		// if pose is membrane pose and jump is undefined, use membrane jump as default
		if ( jumpnum_ == 0 && pose.conformation().is_membrane() ) {
			jumpnum_ = pose.conformation().membrane_info()->membrane_jump();
		} else if ( jumpnum_ == 0 && ! pose.conformation().is_membrane() ) {
			jumpnum_ = 2;
		}

		// get upstream stub ( =MEM for fixed membrane or pose for fixed pose )
		core::kinematics::Stub up_stub = pose.conformation().upstream_jump_stub( jumpnum_ );

		// get length of vector
		Real length = translation_vector_.length();

		// create jump, translate it along axis
		core::kinematics::Jump flexible_jump = pose.jump( jumpnum_ );
		flexible_jump.translation_along_axis( up_stub, translation_vector_, length );
		pose.set_jump( jumpnum_, flexible_jump );

		// reset foldtree and show final one
		pose.fold_tree( orig_ft );
		TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
		pose.fold_tree().show( TR );
	}

}// apply

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
/// @details Register mover-relevant options with JD2 - includes
/// mp:setup options: center
void
TranslationMover::register_options() {

	using namespace basic::options;

	option.add_relevant( OptionKeys::mp::setup::center );

}

/// @brief Initialize Mover options from the commandline
/// @details Initialize mover settings from the commandline
/// using the mp:setup center flag
void
TranslationMover::init_from_cmd() {

	using namespace basic::options;

	// Read in Center Parameter
	// TODO: Is this correct??? Shouldn't the translation vector be new_center minus current_center???
	if ( option[ OptionKeys::mp::setup::center ].user() ) {
		translation_vector_.x() = option[ OptionKeys::mp::setup::center ]()[1];
		translation_vector_.y() = option[ OptionKeys::mp::setup::center ]()[2];
		translation_vector_.z() = option[ OptionKeys::mp::setup::center ]()[3];
	}
}// init from cmd

/// @brief Rotates the pose such that a vector in the old orientation will be
///  overlayed in the new orientation. Requires two vectors (old and new)
///  and a point on the new vector around which the rotation takes place.
///  The mover is a general mover but used mainly on membrane proteins.
///  For membrane proteins, the two vectors will be the old and new normal
///  and the point will be the new center. The default jump is going to be
///  the membrane jump, but you can also specify your own. The rotation
///  happens around the axis perpendicular to both vectors with an angle
///  enclosed by both vectors.
RotationMover::RotationMover() :
	protocols::moves::Mover(),
	old_normal_( 0, 0, 1 ),
	new_normal_( 0, 1, 0 ),
	rot_center_( 0, 0, 0 ),
	jumpnum_(0)
{
	register_options();
	init_from_cmd();
}

/// @brief Custom Constructor
/// @details User can specify an old normal, a new normal, and a new center
///   around which the rotation takes place
RotationMover::RotationMover(
	Vector old_normal,
	Vector new_normal,
	Vector rot_center
) :
	protocols::moves::Mover(),
	old_normal_( old_normal ),
	new_normal_( new_normal ),
	rot_center_( rot_center ),
	jumpnum_(0)
{
	register_options();
	init_from_cmd();
}

/// @brief Custom Constructor
/// @details User can specify an old normal, a new normal, and a new center
///   around which the rotation takes place on this particular jump;
///   operation happens on the downstream stub
RotationMover::RotationMover(
	Vector old_normal,
	Vector new_normal,
	Vector rot_center,
	Size jumpnum
) :
	protocols::moves::Mover(),
	old_normal_( old_normal ),
	new_normal_( new_normal ),
	rot_center_( rot_center ),
	jumpnum_( jumpnum )
{
	register_options();
	init_from_cmd();
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
RotationMover::RotationMover( RotationMover const & src ) :
	protocols::moves::Mover( src ),
	old_normal_( src.old_normal_ ),
	new_normal_( src.new_normal_ ),
	rot_center_( src.rot_center_ ),
	jumpnum_( src.jumpnum_ )
{}

/// @brief Destructor
RotationMover::~RotationMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
RotationMover::clone() const {
	return ( protocols::moves::MoverOP( new RotationMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
RotationMover::fresh_instance() const {
	return protocols::moves::MoverOP( new RotationMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
RotationMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	// read center and normal tag
	read_center_normal_from_tag( rot_center_, new_normal_, tag );

}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
RotationMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RotationMover() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
RotationMoverCreator::keyname() const {
	return RotationMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
RotationMoverCreator::mover_name() {
	return "RotationMover";
}


/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover
std::string
RotationMover::get_name() const {
	return "RotationMover";
}

/// @brief Translate the pose along the defined vector
void
RotationMover::apply( Pose & pose ) {

	using namespace core::conformation::membrane;
	using namespace protocols::membrane::geometry;

	// checking input
	core::Vector diff = new_normal_ - old_normal_;
	core::Vector zero( 0, 0, 0 );
	if ( diff == zero ) {
		TR << "WARNING: Old and new normal are identical. Skipping rotation!" << std::endl;
	} else {
		TR << "Rotating the pose..." << std::endl;

		// starting foldtree
		TR << "Starting foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
		pose.fold_tree().show( TR );
		core::kinematics::FoldTree orig_ft = pose.fold_tree();

		// normalize the new normal
		new_normal_.normalize();

		// if pose membrane pose, use membrane jump as default
		if ( jumpnum_ == 0 && pose.conformation().is_membrane() ) {
			jumpnum_ = pose.conformation().membrane_info()->membrane_jump();
		} else if ( jumpnum_ == 0 && ! pose.conformation().is_membrane() ) {
			jumpnum_ = 2;
		}

		// get upstream stub ( =MEM )
		core::kinematics::Stub up_stub = pose.conformation().upstream_jump_stub( jumpnum_ );

		// get rotation axis = perpendicular of two normals
		// order of vectors is crucial!
		Vector rot_axis = cross( old_normal_, new_normal_ );

		// get rotation angle
		Real angle = numeric::conversions::degrees( angle_of( new_normal_, old_normal_ ) );

		// do the actual rotation
		rigid::RigidBodyDeterministicSpinMover spinmover = rigid::RigidBodyDeterministicSpinMover( jumpnum_, rot_axis, rot_center_, angle );
		spinmover.apply( pose );

		// reset foldtree and show final one
		pose.fold_tree( orig_ft );
		TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
		pose.fold_tree().show( TR );
	}

}// apply

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
/// @details Register mover-relevant options with JD2 - includes
/// mp:setup options: center
void
RotationMover::register_options() {

	using namespace basic::options;

	option.add_relevant( OptionKeys::mp::setup::center );
	option.add_relevant( OptionKeys::mp::setup::normal );

}

/// @brief Initialize Mover options from the commandline
/// @details Initialize mover settings from the commandline
/// using the mp:setup center flag
void
RotationMover::init_from_cmd() {

	// read center and normal from cmd
	read_center_normal_from_cmd( rot_center_, new_normal_ );

}// init from cmd

/// @brief TranslationRotation vector can be defined in -mp:setup center
///   flag to translate the new pose to. The mover is a general mover
///   but used mainly on membrane proteins, that's why we use this flag here
///   The default jump is going to be the membrane jump, but you can also
///   specify your own
TranslationRotationMover::TranslationRotationMover() :
	protocols::moves::Mover(),
	old_center_( 0, 0, 0 ),
	old_normal_( 0, 0, 1 ),
	new_center_( 1,  0, 0 ),
	new_normal_(  0, 1, 0 ),
	jumpnum_( 2 )
{
	register_options();
	init_from_cmd();
}

/// @brief Custom Constructor
/// @details User can specify a TranslationRotation vector
TranslationRotationMover::TranslationRotationMover(
	Vector old_center,
	Vector old_normal,
	Vector new_center,
	Vector new_normal
) :
	protocols::moves::Mover(),
	old_center_( old_center ),
	old_normal_( old_normal ),
	new_center_( new_center ),
	new_normal_( new_normal ),
	jumpnum_( 2 )
{
	register_options();
	init_from_cmd();
}

/// @brief Custom Constructor
/// @details User can specify a TranslationRotation vector and jumpnum
TranslationRotationMover::TranslationRotationMover(
	Vector old_center,
	Vector old_normal,
	Vector new_center,
	Vector new_normal,
	Size jumpnum
) :
	protocols::moves::Mover(),
	old_center_( old_center ),
	old_normal_( old_normal ),
	new_center_( new_center ),
	new_normal_( new_normal ),
	jumpnum_( jumpnum )
{
	register_options();
	init_from_cmd();
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
TranslationRotationMover::TranslationRotationMover( TranslationRotationMover const & src ) :
	protocols::moves::Mover(),
	old_center_( src.old_center_ ),
	old_normal_( src.old_normal_ ),
	new_center_( src.new_center_ ),
	new_normal_( src.new_normal_ ),
	jumpnum_( src.jumpnum_ )
{}

/// @brief Destructor
TranslationRotationMover::~TranslationRotationMover() {}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
TranslationRotationMover::clone() const {
	return ( protocols::moves::MoverOP( new TranslationRotationMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
TranslationRotationMover::fresh_instance() const {
	return protocols::moves::MoverOP( new TranslationRotationMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
TranslationRotationMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	// read center and normal tag
	read_center_normal_from_tag( new_center_, new_normal_, tag );

}

/// @brief Create a new copy of this mover
protocols::moves::MoverOP
TranslationRotationMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new TranslationRotationMover() );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
TranslationRotationMoverCreator::keyname() const {
	return TranslationRotationMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
TranslationRotationMoverCreator::mover_name() {
	return "TranslationRotationMover";
}


/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover
std::string
TranslationRotationMover::get_name() const {
	return "TranslationRotationMover";
}

/// @brief Translate the pose along the defined vector
void
TranslationRotationMover::apply( Pose & pose ) {

	using namespace core::conformation::membrane;
	using namespace protocols::membrane::geometry;

	TR << "Translating and rotating the pose" << std::endl;

	// starting foldtree
	TR << "Starting foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );
	core::kinematics::FoldTree orig_ft = pose.fold_tree();

	// normalize new normal
	new_normal_.normalize();

	// if pose membrane pose, use membrane jump as default
	if ( jumpnum_ == 0 && pose.conformation().is_membrane() ) {
		jumpnum_ = pose.conformation().membrane_info()->membrane_jump();
	} else if ( jumpnum_ == 0 && ! pose.conformation().is_membrane() ) {
		jumpnum_ = 2;
	}

	// translate pose
	Vector translation_vector = new_center_ - old_center_;
	if ( translation_vector.length() != 0 ) {
		TranslationMoverOP translate( new TranslationMover( translation_vector, jumpnum_ ) );
		translate->apply( pose );
	}

	// rotate pose
	RotationMoverOP rotate( new RotationMover( old_normal_, new_normal_, new_center_, jumpnum_ ) );
	rotate->apply( pose );

	// reset foldtree and show final one
	pose.fold_tree( orig_ft );
	TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

}// apply

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
/// @details Register mover-relevant options with JD2 - includes
/// mp:setup options: center
void
TranslationRotationMover::register_options() {

	using namespace basic::options;

	option.add_relevant( OptionKeys::mp::setup::center );
	option.add_relevant( OptionKeys::mp::setup::normal );

}

/// @brief Initialize Mover options from the commandline
/// @details Initialize mover settings from the commandline
/// using the mp:setup center flag
void
TranslationRotationMover::init_from_cmd() {

	// read center and normal from cmd
	read_center_normal_from_cmd( new_center_, new_normal_ );

}// init from cmd


} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_TranslationRotationMover_cc
