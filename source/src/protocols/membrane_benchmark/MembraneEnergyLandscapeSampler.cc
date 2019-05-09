// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane_benchmark/MembraneEnergyLandscapeSampler.cc
/// @brief Sample a protein energy landscape as a function of orientation in the membrane
/// @details Orientation is sampled as a function of two degrees of freedom: helix tilt
/// relative to the membrane normal and distance between the membrane center and peptide residue
/// center of mass.
/// @author Rebecca F. Alford (rfalford12@gmail.com)

// Unit headers
#include <protocols/membrane_benchmark/MembraneEnergyLandscapeSampler.hh>
#include <protocols/membrane_benchmark/MembraneEnergyLandscapeSamplerCreator.hh>

// Pakcage headers
#include <protocols/membrane/TransformIntoMembraneMover.hh>
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
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>

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

static basic::Tracer TR( "protocols.membrane_benchmark.MembraneEnergyLandscapeSampler" );

namespace protocols {
namespace membrane_benchmark {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
MembraneEnergyLandscapeSampler::MembraneEnergyLandscapeSampler():
	protocols::moves::Mover( MembraneEnergyLandscapeSampler::mover_name() ),
	sfxn_(),
	sfxn_weights_(),
	rotation_type_( "YZ" ),
	interface_( false )
{}

/// @brief Non-defualt cnstructor for landscape sampling
MembraneEnergyLandscapeSampler::MembraneEnergyLandscapeSampler(
	std::string sfxn_weights,
	std::string rotation_type,
	bool interface
) : protocols::moves::Mover( MembraneEnergyLandscapeSampler::mover_name() ),
	sfxn_( core::scoring::ScoreFunctionFactory::create_score_function( sfxn_weights ) ),
	sfxn_weights_( sfxn_weights ),
	rotation_type_( rotation_type ),
	interface_( interface )
{}

/// @brief Copy constructor
MembraneEnergyLandscapeSampler::MembraneEnergyLandscapeSampler( MembraneEnergyLandscapeSampler const & src ):
	protocols::moves::Mover( src ),
	sfxn_( src.sfxn_ ),
	sfxn_weights_( src.sfxn_weights_ ),
	rotation_type_( src.rotation_type_ ),
	interface_( src.interface_ )
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
MembraneEnergyLandscapeSampler::~MembraneEnergyLandscapeSampler() {}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Apply the mover
void
MembraneEnergyLandscapeSampler::apply( core::pose::Pose & pose ) {

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

	// Perform an initial transformation of the pose into the membrane
	TransformIntoMembraneMoverOP transform_into_memb( new TransformIntoMembraneMover() );
	transform_into_memb->apply( pose );

	TR << "Initialize energy landscape sampling with rotation type " << rotation_type_ << std::endl;

	// Setup energy function based on user specified weights file
	TR << "Creating an energy function with weights " << sfxn_weights_ << std::endl;

	// Setup output filename and Header
	TR << "Configuring output file with 2D lanscape data" << std::endl;
	utility::vector1< std::string > temp( utility::string_split( pose.pdb_info()->name(), '/') );
	std::string tempstr = temp[ temp.size() ].substr(0, temp[ temp.size() ].size()-4 );
	std::string filename( tempstr +  "_" + sfxn_->get_name() + "_landscape.dat" );
	utility::io::ozstream output( filename );

	// Write header, including score types available in the given energy function
	output << "zcoord angle total_score ";
	ScoreTypes const nonzero_types( sfxn_->get_nonzero_weighted_scoretypes() );
	for ( core::Size ii = 1; ii <= nonzero_types.size(); ++ii ) {
		output << name_from_score_type( nonzero_types[ii] ) << " ";
	}
	output << std::endl;

	// Create a pack rotamers mover
	using namespace protocols::minimization_packing;
	using namespace core::pack::task;
	PackerTaskOP pack_task( TaskFactory::create_packer_task( pose ) );
	PackRotamersMoverOP pack_mover( new PackRotamersMover( sfxn_, pack_task ) );

	// Calculate the protein length to determine bounds
	TR << "Computing helix length to determine landscape boundaries" << std::endl;
	core::Real start_z( pose.residue( 1 ).xyz( 2 ).z() );
	core::Real end_z( pose.residue( pose.total_residue() ).xyz( 2 ).z() );
	core::Real length( std::abs( start_z - end_z ) );
	core::Real limit( 0 );
	if ( interface_ ) {
		limit = 70;
	} else {
		limit = 15 + 3*length;
	}

	// get the membrane jump
	Size membrane_jump( pose.conformation().membrane_info()->membrane_jump() );

	// Apply an initial translation of the membrane center
	Vector initial_move( 0, 0, -limit );
	TranslationMoverOP initialize( new TranslationMover( initial_move, membrane_jump ) );
	initialize->apply( pose );
	Vector init_center( pose.conformation().membrane_info()->membrane_center( pose.conformation() ) );
	TR << "Shifting the protein to a new center position of (" << initial_move.x() << "," << initial_move.y() << "," << initial_move.z() << ")" << std::endl;

	// Set up deltas for rotations & translations and spin axes
	// We densely sample the grid to more easily inspect for discontinuities
	core::Real rotation_delta( 1 ); // Degrees
	core::Vector translation_delta( 0, 0, 1 ); // Angstroms
	TR << "Setting a translation delta of: " << translation_delta.z() << std::endl;
	TR << "Setting a rotation delta of " << rotation_delta << std::endl;

	// Set up a constant translation mover along the z axis
	TranslationMoverOP translate_memb( new TranslationMover( translation_delta, membrane_jump ) );

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

	// Translate along the z axis in very small steps
	for ( core::Real z_coord = -limit; z_coord <= limit; z_coord += translation_delta.z() ) {

		// Define the rotation center from the transmembrane center of mass
		Vector rot_center = pose_tm_com( pose );

		// Rotate from zero to 360 degrees in very small steps
		for ( core::Real normal_angle = 0; normal_angle <= 360; normal_angle += rotation_delta ) {

			// Verify the membrane information
			pose.conformation().membrane_info()->membrane_center( pose.conformation() ).show();

			// Make a copy of the pose to avoid rounding errors
			using namespace core::conformation::membrane;
			core::pose::PoseOP pose_copy = pose.clone();

			// Calculate the rotation based on the pose's current position and degree of rotation
			RotationMoverOP rotate( get_rotation( normal_angle, axis, rot_center, membrane_jump, rotation_type_ ) );
			rotate->apply( *pose_copy );

			// Repack the sidechains of the current configuration
			//pack_mover->apply( *pose_copy );

			// Write data to output file
			output << z_coord << " " << normal_angle;
			output << " " << sfxn_->score( *pose_copy );
			for ( core::Size ii = 1; ii <= nonzero_types.size(); ++ii ) {
				output << " " << sfxn_->score_by_scoretype( *pose_copy, nonzero_types[ii] );
			} // finish writing scores
			output << std::endl;

		} // end rotation loop

		// increment the translation
		translate_memb->apply( pose );

	} // end translation loop

} // apply

protocols::membrane::RotationMoverOP
MembraneEnergyLandscapeSampler::get_rotation( core::Real normal_angle, core::Vector axis, core::Vector rot_center, core::Size membrane_jump, std::string rotation_type ) {

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
MembraneEnergyLandscapeSampler::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
MembraneEnergyLandscapeSampler::parse_my_tag(
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

}

////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
MembraneEnergyLandscapeSampler::fresh_instance() const
{
	return protocols::moves::MoverOP( new MembraneEnergyLandscapeSampler );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
MembraneEnergyLandscapeSampler::clone() const
{
	return protocols::moves::MoverOP( new MembraneEnergyLandscapeSampler( *this ) );
}

std::string MembraneEnergyLandscapeSampler::get_name() const {
	return mover_name();
}

std::string MembraneEnergyLandscapeSampler::mover_name() {
	return "MembraneEnergyLandscapeSampler";
}

void MembraneEnergyLandscapeSampler::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_score_function( attlist, "scorefxn" );
	attlist
		+ XMLSchemaAttribute( "sfxn_weights", xs_string, "Energy function weights file" )
		+ XMLSchemaAttribute( "rotation_type", xs_string, "Rotation axis: XY, XZ, or YZ" )
		+ XMLSchemaAttribute( "interface", xsct_rosetta_bool, "Should I treat this like an interface landscape scanning problem?");
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Sample the membrane energy function landscape using a given energy function", attlist );

}

/////////////// Creator ///////////////

protocols::moves::MoverOP
MembraneEnergyLandscapeSamplerCreator::create_mover() const
{
	return protocols::moves::MoverOP( new MembraneEnergyLandscapeSampler );
}

std::string
MembraneEnergyLandscapeSamplerCreator::keyname() const
{
	return MembraneEnergyLandscapeSampler::mover_name();
}

void MembraneEnergyLandscapeSamplerCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MembraneEnergyLandscapeSampler::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////

std::ostream &
operator<<( std::ostream & os, MembraneEnergyLandscapeSampler const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //membrane_benchmark
