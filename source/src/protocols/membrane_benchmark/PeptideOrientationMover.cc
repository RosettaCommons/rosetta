// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane_benchmark/PeptideOrientationMover.cc
/// @brief Orient a peptide to a specific orientation relative to the membrane
/// @author Rebecca Alford (rfalford12@gmail.com)

// Unit headers
#include <protocols/membrane_benchmark/PeptideOrientationMover.hh>
#include <protocols/membrane_benchmark/PeptideOrientationMoverCreator.hh>

// Pakcage headers
#include <protocols/membrane/TranslationRotationMover.hh>
#include <protocols/membrane/util.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

#include <numeric/conversions.hh>
#include <numeric/constants.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/down_cast.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

// C++ Headers
#include <cstdlib>
#include <ostream>
#include <fstream>

static basic::Tracer TR( "protocols.membrane_benchmark.PeptideOrientationMover" );

namespace protocols {
namespace membrane_benchmark {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PeptideOrientationMover::PeptideOrientationMover():
	protocols::moves::Mover( PeptideOrientationMover::mover_name() ),
	sfxn_(),
	sfxn_weights_(),
	rotation_type_( "YZ" ),
	zcoord_( 0.0 ),
	angle_( 0.0 ),
	interface_( false )
{}

/// @brief Non-defualt cnstructor for landscape sampling
PeptideOrientationMover::PeptideOrientationMover(
	std::string sfxn_weights,
	std::string rotation_type,
	core::Real zcoord,
	core::Real angle,
	bool interface
) : protocols::moves::Mover( PeptideOrientationMover::mover_name() ),
	sfxn_( core::scoring::ScoreFunctionFactory::create_score_function( sfxn_weights ) ),
	sfxn_weights_( sfxn_weights ),
	rotation_type_( rotation_type ),
	zcoord_( zcoord ),
	angle_( angle ),
	interface_( interface )
{}

/// @brief Copy constructor
PeptideOrientationMover::PeptideOrientationMover( PeptideOrientationMover const & src ):
	protocols::moves::Mover( src ),
	sfxn_( src.sfxn_ ),
	sfxn_weights_( src.sfxn_weights_ ),
	rotation_type_( src.rotation_type_ ),
	zcoord_( src.zcoord_ ),
	angle_( src.angle_ ),
	interface_( src.interface_ )
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PeptideOrientationMover::~PeptideOrientationMover() {}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
PeptideOrientationMover::apply( core::pose::Pose & pose ) {

	using namespace numeric;
	using namespace core;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace protocols::membrane;
	using namespace protocols::membrane::geometry;

	// Check that the pose is a membrane protein before continuing
	if ( !pose.conformation().is_membrane() ) {
		utility_exit_with_message( "Pose is not a membrane protein, and I cannot score the membrane energy landscape of a non membrane protein!");
	}

	// Abort if the membrane is not fixed. Leads to significant rounding errors
	if ( is_membrane_fixed( pose ) ) {
		TR << "Membrane is fixed and the pose is moveable" << std::endl;
	} else {
		utility_exit_with_message( "The pose is moveable and the membrane is fixed. Exiting..." );
	}

	TR << "Initialize energy landscape sampling with rotation type " << rotation_type_ << std::endl;

	// Get the membrane jump
	Size membrane_jump( pose.conformation().membrane_info()->membrane_jump() );

	// Apply an initial translation of the membrane center
	Vector initial_move( 0, 0, zcoord_ );
	TranslationMoverOP initialize( new TranslationMover( initial_move, membrane_jump ) );
	initialize->apply( pose );
	Vector init_center( pose.conformation().membrane_info()->membrane_center( pose.conformation() ) );
	TR << "Shifting the protein to a new center position of (" << initial_move.x() << "," << initial_move.y() << "," << initial_move.z() << ")" << std::endl;

	// Setup the axis of rotation (formerly the membrane normal
	core::Vector axis(0,0,0);
	if ( rotation_type_ == "YZ" || rotation_type_ == "XZ" ) {
		axis.z() = 1;
	} else if ( rotation_type_ == "XY" ) {
		axis.x() = 1;
	} else {
		utility_exit_with_message( "Unknown rotation type " + rotation_type_ );
	}

	// If testing interface sampling, rotate the pose by 90 degrees and reset the rotation axis
	if ( interface_ ) {

		TR << "Interface Mode: Rotate the pose by 90 degrees in the YZ plane so it lies horizontal to the membrane plane" << std::endl;

		// this step gets the peptide onto its side
		core::Vector rot_center = pose_tm_com( pose );
		core::Real delta_degrees( 90 );
		core::Real rotation_delta( numeric::conversions::radians( delta_degrees ) );
		Vector new_normal( 0, axis.y()*cos(rotation_delta) + axis.z()*sin(rotation_delta), -axis.y()*sin(rotation_delta) + axis.z()*cos(rotation_delta) );
		RotationMoverOP rotate_onto_side( new RotationMover( axis, new_normal, rot_center, membrane_jump ) );
		rotate_onto_side->apply( pose );

	}

	// Define the rotation center from the transmembrane center of mass
	Vector rot_center = pose_tm_com( pose );

	// Verify the membrane information
	pose.conformation().membrane_info()->membrane_center( pose.conformation() ).show();

	// Make a copy of the pose to avoid rounding errors
	using namespace core::conformation::membrane;
	core::pose::PoseOP pose_copy = pose.clone();

	// Calculate the rotation based on the pose's current position and degree of rotation
	RotationMoverOP rotate( get_rotation( angle_, axis, rot_center, membrane_jump, rotation_type_ ) );
	rotate->apply( *pose_copy );

	pose_copy->dump_pdb( "oriented_pose.pdb" );

} // apply

protocols::membrane::RotationMoverOP
PeptideOrientationMover::get_rotation( core::Real normal_angle, core::Vector axis, core::Vector rot_center, core::Size membrane_jump, std::string rotation_type ) {

	using namespace protocols::membrane;
	using namespace numeric::conversions;

	// Pick the rotation delta from the normal angle
	core::Real delta( radians( normal_angle ) );

	// Calculate a new normal based on an XY or XZ rotation
	core::Vector updated_axis( 0, 0, 0 );
	if ( rotation_type == "XZ" ) {
		updated_axis.x() = axis.x() * cos(delta) + axis.z() * sin(delta);
		updated_axis.y() = 0;
		updated_axis.z() = -axis.x() * sin(delta) + axis.z() * cos(delta);
	} else if ( rotation_type == "YZ" ) {
		updated_axis.x() = 0;
		updated_axis.y() = axis.y() * cos(delta) + axis.z() * sin(delta);
		updated_axis.z() = -axis.y() * sin(delta) + axis.z() * cos(delta);
	} else if ( rotation_type == "XY" ) {
		updated_axis.x() = axis.x() * cos(delta) + axis.y() * sin(delta);
		updated_axis.y() = -axis.x() * sin(delta) + axis.y() * cos(delta);
		updated_axis.z() = 0;
	}

	RotationMoverOP rotation( new RotationMover( axis, updated_axis, rot_center, membrane_jump ) );
	return rotation;

}


////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
PeptideOrientationMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
PeptideOrientationMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

	// Get scoring function and set the weights
	using namespace protocols::rosetta_scripts;
	using namespace core::scoring;

	ScoreFunctionOP candidate_sfxn( parse_score_function( tag, "scorefxn", data ) );
	if ( candidate_sfxn == nullptr ) {
		utility_exit_with_message("No energy function specified!" );
	} else {
		sfxn_ = candidate_sfxn;
		sfxn_weights_ = sfxn_->get_name();
		TR << "Setting up a new scoring function with weights " << sfxn_->get_name() << std::endl;
	}

	// Specify rotation type
	if ( tag->hasOption( "rotation_type" ) ) {
		rotation_type_ = tag->getOption< std::string >( "rotation_type" );
	}

	// Should I treat this like an interface landscape problem?
	if ( tag->hasOption( "interface" ) ) {
		interface_ = tag->getOption< bool >( "interface" );
	}

	// Specify the angle and z coordinate
	if ( tag->hasOption( "zcoord" ) ) {
		zcoord_ = tag->getOption< core::Real >( "zcoord" );
	}

	if ( tag->hasOption( "angle" ) ) {
		angle_ = tag->getOption< core::Real >( "angle" );
	}

}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PeptideOrientationMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new PeptideOrientationMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
PeptideOrientationMover::clone() const
{
	return protocols::moves::MoverOP( new PeptideOrientationMover( *this ) );
}

std::string PeptideOrientationMover::get_name() const {
	return mover_name();
}

std::string PeptideOrientationMover::mover_name() {
	return "PeptideOrientationMover";
}

void PeptideOrientationMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_score_function( attlist, "scorefxn" );
	attlist
		+ XMLSchemaAttribute( "sfxn_weights", xs_string, "Energy function weights file" )
		+ XMLSchemaAttribute( "rotation_type", xs_string, "Rotation axis: XY, XZ, or YZ" )
		+ XMLSchemaAttribute( "zcoord", xsct_real, "Z coordinate" )
		+ XMLSchemaAttribute( "angle", xsct_real, "Angle" )
		+ XMLSchemaAttribute( "interface", xsct_rosetta_bool, "Should I treat this like an interface landscape scanning problem?");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Sample the membrane energy function landscape using a given energy function", attlist );

}

/////////////// Creator ///////////////

protocols::moves::MoverOP
PeptideOrientationMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new PeptideOrientationMover );
}

std::string
PeptideOrientationMoverCreator::keyname() const
{
	return PeptideOrientationMover::mover_name();
}

void PeptideOrientationMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	PeptideOrientationMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

std::ostream &
operator<<( std::ostream & os, PeptideOrientationMover const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //membrane_benchmark
