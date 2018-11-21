// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farnesyl/InstallFarnesylMover.cc
/// @brief Modifies a free cysteine residue with a branch of 3 DMA residues (the terpene monomer) to create farnesyl-cysteine
/// @author Andy Watkins (andy.watkins2@gmail.com)

// Unit headers
#include <protocols/farnesyl/InstallFarnesylMover.hh>
#include <protocols/farnesyl/InstallFarnesylMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// These are only used in one case -- where we want to
// do some amount of sampling the second that we install.
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerComb.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerOneTorsion.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility/tools/make_vector1.hh>
#include <core/id/TorsionID.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/rosetta_scripts/util.hh>

static basic::Tracer TR( "protocols.farnesyl.InstallFarnesylMover" );

using namespace core;
using namespace core::id;
using namespace utility::tools;

namespace protocols {
namespace farnesyl {

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
InstallFarnesylMover::InstallFarnesylMover():
	protocols::moves::Mover( InstallFarnesylMover::mover_name() )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
InstallFarnesylMover::InstallFarnesylMover( InstallFarnesylMover const & src ):
	protocols::moves::Mover( src ),
	selector_( src.selector_->clone() ),
	sample_per_residue_( src.sample_per_residue_ )
{

}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
InstallFarnesylMover::~InstallFarnesylMover(){}

////////////////////////////////////////////////////////////////////////////////
/// Mover Methods ///
/////////////////////

void
InstallFarnesylMover::sample_first( core::pose::Pose & pose, Size const cys_idx ) {
	auto sfxn_ = core::scoring::get_score_function();

	Size const dma_one_idx = pose.size();
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_bb( cys_idx, true );
	movemap->set_chi( cys_idx, true );
	movemap->set_bb( dma_one_idx, true );
	movemap->set_chi( dma_one_idx, true );

	movemap->set( TorsionID( dma_one_idx, BB, 2 ), false );
	movemap->show( TR );
	protocols::minimization_packing::MinMover minmover( movemap, sfxn_, "lbfgs_armijo_nonmonotone", 0.001, true );
	//minmover.set_movemap( movemap );

	pose.set_torsion( TorsionID( dma_one_idx, BB, 2 ), 180 );

	// Sample cys chis very precisely
	auto chi_rots = make_vector1( -165, -150, -135, -120, -105, -90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180 );
	// Sample sp3 rotamers but also +/- 10 degrees. This is better than sampling every n degrees.
	auto sp3_rots = make_vector1( 50, 60, 70, 170, 180, 190, 290, 300, 310 );
	using namespace protocols::stepwise::sampler;
	StepWiseSamplerComb comb;
	comb.set_random( false );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( cys_idx, CHI, 1 ), chi_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( cys_idx, CHI, 2 ), chi_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_one_idx, BB, 4 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_one_idx, BB, 3 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_one_idx, BB, 1 ), sp3_rots ) ) );
	comb.init();

	TR << "Initted comb. " << std::endl;
	// These arbitrary score cutoffs are a major weakness. They are intended to save minimization work.
	// There are innumerable ways to do this better. This is based entirely on 1kzo (example). Should be
	// based on prior knowledge and settable, probably? Or dynamically figured out.
	Real first_cutoff = 5000;
	Pose best_pose;
	Real best_score = first_cutoff + 100; // greater than first_cutoff so guaranteed irrelevant
	Size ii = 0;
	while ( comb.not_end() ) {
		// We need a local pose so we don't just keep minning the same pose.
		// In theory, we could just copy before minning if this ends up super slow.
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
}

void
InstallFarnesylMover::sample_second( core::pose::Pose & pose ) {
	auto sfxn_ = core::scoring::get_score_function();

	Size const dma_one_idx = pose.size() - 1;
	Size const dma_two_idx = pose.size();
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_bb( dma_one_idx, true );
	movemap->set_chi( dma_one_idx, true );
	movemap->set_bb( dma_two_idx, true );
	movemap->set_chi( dma_two_idx, true );

	movemap->set( TorsionID( dma_one_idx, BB, 2 ), false );
	movemap->set( TorsionID( dma_two_idx, BB, 2 ), false );
	movemap->show( TR );
	protocols::minimization_packing::MinMover minmover( movemap, sfxn_, "lbfgs_armijo_nonmonotone", 0.001, true );
	//minmover.set_movemap( movemap );

	pose.set_torsion( TorsionID( dma_one_idx, BB, 2 ), 180 );
	pose.set_torsion( TorsionID( dma_two_idx, BB, 2 ), 180 );

	// Sample sp3 rotamers but also +/- 10 degrees. This is better than sampling every n degrees.
	auto sp3_rots = make_vector1( 50, 60, 70, 170, 180, 190, 290, 300, 310 );
	using namespace protocols::stepwise::sampler;
	StepWiseSamplerComb comb;
	comb.set_random( false );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_one_idx, BB, 4 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_one_idx, BB, 3 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_one_idx, BB, 1 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_two_idx, BB, 4 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_two_idx, BB, 3 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_two_idx, BB, 1 ), sp3_rots ) ) );
	comb.init();

	TR << "Initted comb. " << std::endl;
	// These arbitrary score cutoffs are a major weakness. They are intended to save minimization work.
	// There are innumerable ways to do this better. This is based entirely on 1kzo (example). Should be
	// based on prior knowledge and settable, probably? Or dynamically figured out.
	Real first_cutoff = 5000;
	Pose best_pose;
	Real best_score = first_cutoff + 100; // greater than first_cutoff so guaranteed irrelevant
	Size ii = 0;
	while ( comb.not_end() ) {
		// We need a local pose so we don't just keep minning the same pose.
		// In theory, we could just copy before minning if this ends up super slow.
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
}

void
InstallFarnesylMover::sample_third( core::pose::Pose & pose ) {
	auto sfxn_ = core::scoring::get_score_function();

	Size const dma_two_idx = pose.size() - 1;
	Size const dma_three_idx = pose.size();
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_bb( dma_two_idx, true );
	movemap->set_chi( dma_two_idx, true );
	movemap->set_bb( dma_three_idx, true );
	movemap->set_chi( dma_three_idx, true );

	movemap->set( TorsionID( dma_two_idx, BB, 2 ), false );
	movemap->set( TorsionID( dma_three_idx, BB, 2 ), false );
	movemap->show( TR );
	protocols::minimization_packing::MinMover minmover( movemap, sfxn_, "lbfgs_armijo_nonmonotone", 0.001, true );
	//minmover.set_movemap( movemap );

	pose.set_torsion( TorsionID( dma_two_idx, BB, 2 ), 180 );
	pose.set_torsion( TorsionID( dma_three_idx, BB, 2 ), 180 );

	// Sample cys chis very precisely
	auto chi_rots = make_vector1( -165, -150, -135, -120, -105, -90, -75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75, 90, 105, 120, 135, 150, 165, 180 );
	// Sample sp3 rotamers but also +/- 10 degrees. This is better than sampling every n degrees.
	auto sp3_rots = make_vector1( 50, 60, 70, 170, 180, 190, 290, 300, 310 );
	using namespace protocols::stepwise::sampler;
	StepWiseSamplerComb comb;
	comb.set_random( false );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_two_idx, BB, 4 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_two_idx, BB, 3 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_two_idx, BB, 1 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_three_idx, BB, 4 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_three_idx, BB, 3 ), sp3_rots ) ) );
	comb.add_external_loop_rotamer( StepWiseSamplerOneTorsionOP( new StepWiseSamplerOneTorsion( TorsionID( dma_three_idx, BB, 1 ), sp3_rots ) ) );
	comb.init();

	TR << "Initted comb. " << std::endl;
	// These arbitrary score cutoffs are a major weakness. They are intended to save minimization work.
	// There are innumerable ways to do this better. This is based entirely on 1kzo (example). Should be
	// based on prior knowledge and settable, probably? Or dynamically figured out.
	Real first_cutoff = 5000;
	Pose best_pose;
	Real best_score = first_cutoff + 100; // greater than first_cutoff so guaranteed irrelevant
	Size ii = 0;
	while ( comb.not_end() ) {
		// We need a local pose so we don't just keep minning the same pose.
		// In theory, we could just copy before minning if this ends up super slow.
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
}

/// @brief Apply the mover
void
InstallFarnesylMover::apply( core::pose::Pose & pose ) {

	using namespace core::id;

	auto cys_subset = selector_->apply( pose );
	utility::vector1< Size > cys_pos;
	for ( Size ii = 1; ii <= cys_subset.size(); ++ii ) {
		if ( cys_subset[ ii ] ) cys_pos.push_back( ii );
	}

	for ( Size const cys_seqpos : cys_pos ) {
		runtime_assert( pose.residue_type( cys_seqpos ).aa() == core::chemical::aa_cys );

		auto chm = core::chemical::ChemicalManager::get_instance();
		auto rts = chm->residue_type_set( "fa_standard" );
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::SC_BRANCH_POINT, cys_seqpos );

		auto prenyl_type = core::chemical::ResidueTypeFinder( *rts ).name3( "DMA" ).get_representative_type();
		auto prenyl = core::conformation::Residue( prenyl_type, true );

		pose.append_residue_by_atoms( prenyl, true, "C4", cys_seqpos, "SG" );
		pose.conformation().set_bond_length( AtomID( pose.residue_type( cys_seqpos ).atom_index( "SG" ), cys_seqpos ), AtomID( pose.residue_type( pose.size() ).atom_index( "C4" ), pose.size() ), 1.795 );
		pose.conformation().set_bond_angle( AtomID( pose.residue_type( cys_seqpos ).atom_index( "CB" ), cys_seqpos ), AtomID( pose.residue_type( cys_seqpos ).atom_index( "SG" ), cys_seqpos ), AtomID(pose.residue_type( pose.size() ).atom_index( "C4" ), pose.size() ), 1.7610372 ); // 100.9 degrees in radians
		// Sample
		if ( sample_per_residue_ ) sample_first( pose, cys_seqpos );

		// What is the new seqpos? In theory we should walk the foldtree, but realistically
		// it is pose.size()
		pose.append_residue_by_atoms( prenyl, true, "C4", pose.size(), "C1" );
		// Sample
		if ( sample_per_residue_ ) sample_second( pose );

		pose.append_residue_by_atoms( prenyl, true, "C4", pose.size(), "C1" );
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::LOWER_TERMINUS_VARIANT, pose.size() );

		// Sample
		if ( sample_per_residue_ ) sample_third( pose );
	}
}

////////////////////////////////////////////////////////////////////////////////
/// @brief Show the contents of the Mover
void
InstallFarnesylMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

////////////////////////////////////////////////////////////////////////////////
/// Rosetta Scripts Support ///
///////////////////////////////

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
InstallFarnesylMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	core::select::residue_selector::ResidueSelectorCOP selector = protocols::rosetta_scripts::parse_residue_selector( tag, datamap );
	if ( selector ) {
		TR << "Setting selector " << selector->get_name() << std::endl;
		set_selector( selector );
	}
	sample_per_residue_ = tag->getOption< bool >( "sample_per_residue", false );
}

void InstallFarnesylMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	core::select::residue_selector::attributes_for_parse_residue_selector_default_option_name( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "sample_per_residue", xsct_rosetta_bool, "Sample each residue as it's added", "false" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "DOCUMENTATION STRING", attlist );
}


////////////////////////////////////////////////////////////////////////////////
/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
InstallFarnesylMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new InstallFarnesylMover );
}

/// @brief required in the context of the parser/scripting scheme
protocols::moves::MoverOP
InstallFarnesylMover::clone() const
{
	return protocols::moves::MoverOP( new InstallFarnesylMover( *this ) );
}

std::string InstallFarnesylMover::get_name() const {
	return mover_name();
}

std::string InstallFarnesylMover::mover_name() {
	return "InstallFarnesylMover";
}



/////////////// Creator ///////////////

protocols::moves::MoverOP
InstallFarnesylMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new InstallFarnesylMover );
}

std::string
InstallFarnesylMoverCreator::keyname() const
{
	return InstallFarnesylMover::mover_name();
}

void InstallFarnesylMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	InstallFarnesylMover::provide_xml_schema( xsd );
}

////////////////////////////////////////////////////////////////////////////////
/// private methods ///
///////////////////////


std::ostream &
operator<<( std::ostream & os, InstallFarnesylMover const & mover )
{
	mover.show(os);
	return os;
}

} //protocols
} //farnesyl
