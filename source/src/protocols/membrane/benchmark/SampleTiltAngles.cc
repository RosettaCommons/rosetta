// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/membrane/benchmark/SampleTiltAngles.cc
/// @brief Calculates the energy at all possible tilt angles (0->180 degrees)
/// @author Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <protocols/membrane/benchmark/SampleTiltAngles.hh>
#include <protocols/membrane/benchmark/SampleTiltAnglesCreator.hh>

// Package Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/util.hh>

#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

// Project Headers
#include <protocols/rigid/RigidBodyMover.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>

// Utility Headers
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

#include <numeric/xyzVector.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>
#include <ostream>
#include <fstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.membrane.benchmark.SampleTiltAngles" );

namespace protocols {
namespace membrane {
namespace benchmark {

SampleTiltAngles::SampleTiltAngles():
	protocols::moves::Mover( "SampleTiltAngles" ),
	prefix_( "helix" ),
	ref_sfxn1_(),
	ref_sfxn2_(),
	ref_sfxn3_()
{
	using namespace core::scoring;

	ref_sfxn1_ = get_score_function();
	ref_sfxn2_ = ScoreFunctionFactory::create_score_function( "mpframework_fa_2007" );
	ref_sfxn3_ = ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012" );

}

SampleTiltAngles::SampleTiltAngles(
	std::string prefix,
	core::scoring::ScoreFunctionOP ref_sfxn1,
	core::scoring::ScoreFunctionOP ref_sfxn2,
	core::scoring::ScoreFunctionOP ref_sfxn3
):
	protocols::moves::Mover( "SampleTiltAngles" ),
	prefix_( prefix ),
	ref_sfxn1_( ref_sfxn1 ),
	ref_sfxn2_( ref_sfxn2 ),
	ref_sfxn3_( ref_sfxn3 )
{}

SampleTiltAngles::~SampleTiltAngles(){}

SampleTiltAngles::SampleTiltAngles( SampleTiltAngles const & src ):
	protocols::moves::Mover( src ),
	prefix_( src.prefix_ ),
	ref_sfxn1_( src.ref_sfxn1_ ),
	ref_sfxn2_( src.ref_sfxn2_ ),
	ref_sfxn3_( src.ref_sfxn3_ )
{}

void
SampleTiltAngles::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

	using namespace core::scoring;

	// Parse prefix tag
	if ( tag->hasOption( "prefix" ) ) {
		prefix_ = tag->getOption< std::string >( "prefix" );
	}

	// Parse for scoring functions
	// ref1
	if ( tag->hasOption( "ref1" ) ) {
		std::string const ref1( tag->getOption< std::string >( "ref1" ) );
		ref_sfxn1_ = data.get< ScoreFunction * >( "scorefxns", ref1 )->clone();
	}

	// ref2
	if ( tag->hasOption( "ref2" ) ) {
		std::string const ref2( tag->getOption< std::string >( "ref2" ) );
		ref_sfxn2_ = data.get< ScoreFunction * >( "scorefxns", ref2 )->clone();
	}

	// ref3
	if ( tag->hasOption( "ref3" ) ) {
		std::string const ref3( tag->getOption< std::string >( "ref3" ) );
		ref_sfxn3_ = data.get< ScoreFunction * >( "scorefxns", ref3 )->clone();
	}




}

protocols::moves::MoverOP
SampleTiltAngles::clone() const{
	return protocols::moves::MoverOP( new SampleTiltAngles( *this ) );
}

moves::MoverOP
SampleTiltAngles::fresh_instance() const
{
	return protocols::moves::MoverOP( new SampleTiltAngles );
}

// XRW TEMP std::string
// XRW TEMP SampleTiltAngles::get_name() const {
// XRW TEMP  return "SampleTiltAngles";
// XRW TEMP }

void
SampleTiltAngles::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
	output << "Ouptut Prefix: " << prefix_ << std::endl;
	output << "Reference Energy Function 1: " << ref_sfxn1_->get_name() << std::endl;
	output << "Reference Energy Function 2: " << ref_sfxn2_->get_name() << std::endl;
	output << "Reference Energy Function 3: " << ref_sfxn3_->get_name() << std::endl;

}


void
SampleTiltAngles::apply( core::pose::Pose& pose ){

	using namespace core;
	using namespace core::scoring;
	using namespace protocols::membrane::geometry;
	using namespace protocols::rigid;
	using namespace numeric;

	// Determine the membrane jump number (pre-determined in the framework)
	Size membrane_jump( pose.conformation().membrane_info()->membrane_jump() );
	Vector center( pose.conformation().membrane_info()->membrane_center( pose.conformation() ) );
	TR << "Setting up an initial membrane orientation by aligning the peptide to the current membrane normal" << std::endl;
	// Calculate the current helix axes based on the COM of the transmembrane
	// span start and end posiitons
	numeric::xyzVector< core::Real > helix_axis = calc_helix_axis( pose, 1 /* only for single span poses */ );
	Vector spin_axis(1, 0, 0);

	// Rotate the helix to align with the membrane normal (representing a tilt angle of 0 degrees)
	EmbeddingDefOP embedding = compute_structure_based_embedding( pose );
	embedding->set_normal( helix_axis );

	// Transform the helix into a membrane coordinate frame such that the helix axis is parallel
	// to the membrane normal.
	TransformIntoMembraneMoverOP transform_into_memb( new TransformIntoMembraneMover( embedding ) );
	transform_into_memb->apply( pose );

	// Intermediate checks!
	TR << calc_helix_tilt_angle( pose, 1 /* only working with 1 helix poses for now */ ) << std::endl;

	core::Real ref1_score( ref_sfxn1_->score( pose ) );
	core::Real ref2_score( ref_sfxn2_->score( pose ) );
	core::Real ref3_score( ref_sfxn3_->score( pose ) );

	// Setup score & angle arrays
	utility::vector1< core::Real > angles;
	utility::vector1< core::Real > ref1_scores;
	utility::vector1< core::Real > ref2_scores;
	utility::vector1< core::Real > ref3_scores;

	// Add the first scores and angles
	angles.push_back( 0 ); // angle should be zero or close to it
	ref1_scores.push_back( ref1_score );
	ref2_scores.push_back( ref2_score );
	ref3_scores.push_back( ref3_score );

	// Reset the deterministic step size to 1
	// Rotate the helix to align with the membrane normal (representing a tilt angle of 0 degrees)
	core::Real step_size( 1.0 );
	RigidBodyDeterministicSpinMoverOP deterministic_z( new RigidBodyDeterministicSpinMover( membrane_jump, spin_axis, center, step_size ) );

	// For each tilt angle (sample every degree for a smooth surface), score the protein and add to arrays
	for ( core::Size angle = 1; angle <= 180; ++angle ) {

		deterministic_z->apply( pose );

		// Intermediate checks!
		TR << calc_helix_tilt_angle( pose, 1 /* only working with 1 helix poses for now */ ) << std::endl;

		// Get the current state of the pose
		ref1_score = ref_sfxn1_->score( pose );
		ref2_score = ref_sfxn2_->score( pose );
		ref3_score = ref_sfxn3_->score( pose );

		// Add the first scores and angles
		angles.push_back( angle ); // angle should be zero or close to it
		ref1_scores.push_back( ref1_score );
		ref2_scores.push_back( ref2_score );
		ref3_scores.push_back( ref3_score );
	}

	// Write scores to output
	write_score_to_outfiles( angles, ref1_scores, ref2_scores, ref3_scores );
	TR << "Sample Tilt Anlge Protocol Complete!" << std::endl;
}

void
SampleTiltAngles::write_score_to_outfiles(
	utility::vector1< core::Real > angles,
	utility::vector1< core::Real > ref1_scores,
	utility::vector1< core::Real > ref2_scores,
	utility::vector1< core::Real > ref3_scores
) {

	using namespace std;

	// Open three output files for writing - one for each energy function type
	std::ofstream ref1_out( (prefix_ + "_talaris.out").c_str() );
	std::ofstream ref2_out( (prefix_ + "_membrane07.out").c_str() );
	std::ofstream ref3_out( (prefix_ + "_membrane12.out").c_str() );

	// Write the score/angle pairs to an output file
	ref1_out << "TiltAngle Score" << std::endl;
	ref2_out << "TiltAngle Score" << std::endl;
	ref3_out << "TiltAngle Score" << std::endl;

	// Check that angles and score lists are of equal length
	debug_assert( angles.size() == ref1_scores.size() );
	debug_assert( angles.size() == ref2_scores.size() );
	debug_assert( angles.size() == ref3_scores.size() );

	for ( core::Size i = 1; i <= angles.size(); ++i ) {
		ref1_out << angles[i] << " " << ref1_scores[i] << std::endl;
		ref2_out << angles[i] << " " << ref2_scores[i] << std::endl;
		ref3_out << angles[i] << " " << ref3_scores[i] << std::endl;
	}

	ref1_out.close();
	ref2_out.close();
	ref3_out.close();
}

/////////////// Creator ///////////////

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SampleTiltAnglesCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SampleTiltAngles );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SampleTiltAnglesCreator::keyname() const {
// XRW TEMP  return SampleTiltAngles::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SampleTiltAngles::mover_name(){
// XRW TEMP  return "SampleTiltAngles";
// XRW TEMP }

std::string SampleTiltAngles::get_name() const {
	return mover_name();
}

std::string SampleTiltAngles::mover_name() {
	return "SampleTiltAngles";
}

void SampleTiltAngles::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute("prefix", xs_string, "Prefix for output files.")
		+ XMLSchemaAttribute("ref1", xs_string, "Name attribute of previously defined score function to use as first reference score function.")
		+ XMLSchemaAttribute("ref2", xs_string, "Name attribute of previously defined score function to use as second reference score function.")
		+ XMLSchemaAttribute("ref3", xs_string, "Name attribute of previously defined score function to use as third reference score function.");

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Sample different tilt angles of the pose within the membrane.", attlist );
}

std::string SampleTiltAnglesCreator::keyname() const {
	return SampleTiltAngles::mover_name();
}

protocols::moves::MoverOP
SampleTiltAnglesCreator::create_mover() const {
	return protocols::moves::MoverOP( new SampleTiltAngles );
}

void SampleTiltAnglesCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SampleTiltAngles::provide_xml_schema( xsd );
}


} //protocols
} //membrane
} //benchmark


