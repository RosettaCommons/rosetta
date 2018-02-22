// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farnesyl/SampleFarnesylMover.cc
/// @brief Modifies a free cysteine residue with a branch of 3 DMA residues (the terpene monomer) to create farnesyl-cysteine
/// @author Andy Watkins (andy.watkins2@gmail.com)

// Unit headers
#include <protocols/farnesyl/SampleFarnesylMover.hh>
#include <protocols/farnesyl/SampleFarnesylMoverCreator.hh>

#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerComb.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerOneTorsion.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tools/make_vector1.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.farnesyl.SampleFarnesylMover" );

namespace protocols {
namespace farnesyl {

using namespace core;
using namespace core::id;
using namespace core::pose;
using namespace utility::tools;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
SampleFarnesylMover::SampleFarnesylMover():
	protocols::moves::Mover( SampleFarnesylMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
SampleFarnesylMover::SampleFarnesylMover( SampleFarnesylMover const & src ):
	protocols::moves::Mover( src ),
	enumerate_( src.enumerate_ )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
SampleFarnesylMover::~SampleFarnesylMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Sample a particular farnesyl
/// @details Draw torsion combinations from a small library of options. Someone else (TM)
/// can turn this rigorous (with explicit management of their relative probabilities, and
/// perhaps registering these conformations in a useful way, and perhaps expanded and/or
/// made dependent on CYS more explicitly. Or the InstallFarnesylMover could give you a
/// custom library, optionally for this *particular* farnesyl. Plus, this should all be
/// read in from a DB file). What a useful task list for someone else (TM)!
/// Actually -- the sampling process we've been using is pretty fast. I am going to actually
/// use that every time. Why not! The fixed library subset paradigm can be used later. So now
/// we do all this sampling and then take the best conformation forward.
void
SampleFarnesylMover::sample_farnesyl( core::pose::Pose & pose, Size const cys_idx, Size const dma_one_idx, Size const dma_two_idx, Size const dma_three_idx ) {
	if ( !sfxn_ ) sfxn_ = core::scoring::get_score_function();

	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_bb( cys_idx, true );
	movemap->set_chi( cys_idx, true );
	movemap->set_bb( dma_one_idx, true );
	movemap->set_chi( dma_one_idx, true );
	movemap->set_bb( dma_two_idx, true );
	movemap->set_chi( dma_two_idx, true );
	movemap->set_bb( dma_three_idx, true );
	movemap->set_chi( dma_three_idx, true );

	movemap->set( TorsionID( dma_one_idx, BB, 2 ), false );
	movemap->set( TorsionID( dma_two_idx, BB, 2 ), false );
	movemap->set( TorsionID( dma_three_idx, BB, 2 ), false );
	movemap->show( TR );
	protocols::minimization_packing::MinMover minmover( movemap, sfxn_, "lbfgs_armijo_nonmonotone", 0.001, true );
	//minmover.set_movemap( movemap );

	pose.set_torsion( TorsionID( dma_one_idx, BB, 2 ), 180 );
	pose.set_torsion( TorsionID( dma_two_idx, BB, 2 ), 180 );
	pose.set_torsion( TorsionID( dma_three_idx, BB, 2 ), 180 );

	auto sp3_rots = make_vector1( 60, 180, 300 );
	using namespace protocols::stepwise::sampler;
	StepWiseSamplerComb comb;
	comb.set_random( !enumerate_ );
	//comb.add_external_loop_rotamer( protocols::stepwise::sampler::StepWiseSamplerOneTorsion( TorsionID( cys_idx, CHI, 1 ), sp3_rots );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( cys_idx, CHI, 2 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_one_idx, BB, 4 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_one_idx, BB, 3 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_one_idx, BB, 1 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_two_idx, BB, 4 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_two_idx, BB, 3 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_two_idx, BB, 1 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_two_idx, BB, 4 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_two_idx, BB, 3 ), sp3_rots ) ) );
	comb.init();

	TR << "Initted comb. " << std::endl;
	// These arbitrary score cutoffs are a major weakness. They are intended to save minimization work.
	// There are innumerable ways to do this better. This is based entirely on 1kzo (example). Should be
	// based on prior knowledge and settable, probably? Or dynamically figured out.
	Real first_cutoff = 5000;
	Pose best_pose;
	Real best_score = first_cutoff + 100; // greater than first_cutoff so guaranteed irrelevant
	Size ii = 0;
	Size const limit = 500;
	while ( ( comb.not_end() && enumerate_ ) || ( ii < limit && !enumerate_ ) ) {
		// We need a local pose so we don't just keep minning the same pose.
		// In theory, we could just copy before minning if this ends up super slow.

		if ( !enumerate_ ) {
			TR << "Sample " << ii << " of " << limit << std::endl;
		}

		comb.apply( pose );
		++comb; ++ii;

		Real const score = ( *sfxn_ )( pose );
		if ( score > first_cutoff + 10000 ) continue; // generally bump check esque. maybe add a screener?
		TR << "Score for sample " << ii << ": " << score << std::endl;

		if ( score > first_cutoff ) continue;
		minmover.apply( pose );

		Real const min_score = ( *sfxn_ )( pose );
		TR << "Minned core for sample " << ii << ": " << min_score << std::endl;
		if ( min_score < best_score ) {
			best_pose = pose;
			best_score = min_score;
		}
	}

	// Could just transfer relevant dofs if slow.
	pose = best_pose;
}

/// @brief Apply the mover
/// @details Loop through every residue in the pose. If it's not cysteine -- move on.
/// if it has no conjugation -- move on. If it doesn't have 3 pendant DMAs -- move on.
/// Then... you've found farnesyl cysteine; go to the sampling function.
void
SampleFarnesylMover::apply( core::pose::Pose & pose ){
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( pose.residue_type( ii ).aa() != core::chemical::aa_cys ) continue;
		if ( !pose.residue( ii ).has_variant_type( core::chemical::SC_BRANCH_POINT ) ) continue;

		// TODO: Maybe this 'walk' can be performed more smoothly using the FoldTree.
		TR << " SG conn id" << pose.residue_type( ii ).residue_connection_id_for_atom(
			pose.residue_type( ii ).atom_index( "SG" ) ) << std::endl;
		TR << " SG conn res" << pose.residue( ii ).connected_residue_at_resconn( pose.residue_type( ii ).residue_connection_id_for_atom(
			pose.residue_type( ii ).atom_index( "SG" ) ) ) << std::endl;
		Size jj = pose.residue( ii ).connected_residue_at_resconn(
			pose.residue_type( ii ).residue_connection_id_for_atom(
			pose.residue_type( ii ).atom_index( "SG" ) ) );
		if ( jj == 0 || jj > pose.size() || pose.residue_type( jj ).name3() != "DMA" ) continue;

		// The next one is connected by jj's LOWER, which should always be connection 1
		// (even if there come to be permissable connection-modifying variants, they will
		// create bigger connection numbers. And if 1 behaves differently we are
		// ok excluding this 'weird farnesyl')
		Size kk = pose.residue( jj ).connected_residue_at_resconn( 1 );
		if ( kk == 0 || kk > pose.size() || pose.residue_type( kk ).name3() != "DMA" ) continue;

		// Same situation as above.
		Size ll = pose.residue( kk ).connected_residue_at_resconn( 1 );
		if ( ll == 0 || ll > pose.size() || pose.residue_type( ll ).name3() != "DMA" ) continue;

		TR << "Determined seqposes: " << ii << ", " << jj << ", " << kk << ", " << ll << std::endl;
		sample_farnesyl( pose, ii, jj, kk, ll );
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
SampleFarnesylMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
SampleFarnesylMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	set_score_function( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	set_enumerate( tag->getOption< bool >( "enumerate", true ) );
}

void SampleFarnesylMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.
	rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist
		+ XMLSchemaAttribute::required_attribute(
		"enumerate", xsct_rosetta_bool,
		"Enumerate possible farnesyl conformations (there may be a ton!) or just sample a subset from them." );


	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "DOCUMENTATION STRING", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
SampleFarnesylMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new SampleFarnesylMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
SampleFarnesylMover::clone() const
{
	return protocols::moves::MoverOP( new SampleFarnesylMover( *this ) );
}

std::string SampleFarnesylMover::get_name() const {
	return mover_name();
}

std::string SampleFarnesylMover::mover_name() {
	return "SampleFarnesylMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
SampleFarnesylMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new SampleFarnesylMover );
}

std::string
SampleFarnesylMoverCreator::keyname() const
{
	return SampleFarnesylMover::mover_name();
}

void SampleFarnesylMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SampleFarnesylMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, SampleFarnesylMover const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //farnesyl
