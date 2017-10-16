// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DeleteMover
/// @brief Adds an RNA residue from a chain terminus.
/// @details
/// @author Rhiju Das

#include <protocols/stepwise/monte_carlo/mover/AddMover.hh>
#include <protocols/stepwise/monte_carlo/mover/AddMoverCreator.hh>
#include <protocols/stepwise/monte_carlo/rna/RNA_TorsionMover.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/chemical/rna/util.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#include <map>

//Req'd on WIN32
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>

using namespace core;
using namespace core::pose;
using namespace core::pose::full_model_info;


using core::Real;
using utility::make_tag_with_dashes;
using namespace protocols::stepwise::modeler;

//////////////////////////////////////////////////////////////////////////
// Adds one residue to a 5' or 3' chain terminus, and appropriately
//  updates the pose full_model_info object.
//
// Now updated to add whole chunks (if stored in another pose), and also
//  to add nucleotides 'by jump' (a.k.a. the skip-bulge move), and also
//  inter-chain docking.
//
// Could probably be cleaned up -- append vs. prepend are pretty similar,
//  and jump-addition moves are actually identical. Best bet might
//  be to create a pose (even for single nucleotides) and have a single
//  addition function.
//
//////////////////////////////////////////////////////////////////////////

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.monte_carlo.mover.AddMover" );

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {


//////////////////////////////////////////////////////////////////////////
//constructor!
AddMover::AddMover( scoring::ScoreFunctionCOP scorefxn ):
	scorefxn_( scorefxn ),
	rna_torsion_mover_( rna::RNA_TorsionMoverOP( new rna::RNA_TorsionMover ) )
{}

//////////////////////////////////////////////////////////////////////////
//destructor
AddMover::~AddMover()
{}

//////////////////////////////////////////////////////////////////////////
void
AddMover::apply( core::pose::Pose & pose  )
{
	// This syntax is only used by RoesttaScripts and assumes you have
	// already set swa_move_.
	//std::cout << "not defined yet" << std::endl;

	// Not sure -- haven't even checked -- how to best differentiate the
	// two poses.
	apply( pose, swa_move_ );
}

protocols::moves::MoverOP
AddMover::clone() const {
	return AddMoverOP( new AddMover( *this ) );
}

void
AddMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose ) {

	// Oh man, we have access to the pose here! hooray.

	// We have to define the swa_move here, in addition to possibly setting the members
	// that are by default set in the ctor
	scorefxn_  = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
	//scorefxn_ = tag->getOption<bool>("ncbb_design_first", ncbb_design_first_);
	presample_added_residue_ = tag->getOption< bool >( "presample_added_residue", presample_added_residue_ );
	presample_by_swa_ = tag->getOption< bool >( "presample_by_swa", presample_by_swa_ );
	minimize_single_res_ = tag->getOption< bool >( "minimize_single_res", minimize_single_res_ );
	start_added_residue_in_aform_ = tag->getOption< bool >( "start_added_residue_in_aform", start_added_residue_in_aform_ );
	internal_cycles_ = tag->getOption< Size >( "internal_cycles", internal_cycles_ );
	sample_range_small_ = tag->getOption< Real >( "sample_range_small", sample_range_small_ );
	sample_range_large_ = tag->getOption< Real >( "sample_range_large", sample_range_large_ );
	kT_ = tag->getOption< Real >( "kT", kT_ );

	// required.
	std::string move_str = tag->getOption< std::string >( "swa_move" );

	// swa_move_ -- formatted as ADD $move_element $attachments
	// for example, ADD 4 BOND_TO_PREVIOUS 3
	swa_move_ = StepWiseMove( utility::string_split( move_str ), const_full_model_info( pose ).full_model_parameters() );
}

std::string AddMover::get_name() const {
	return mover_name();
}

std::string AddMover::mover_name() {
	return "AddMover";
}


void AddMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "presample_added_residue", xsct_rosetta_bool, "Presample the residue being added", "true" );
	attlist + XMLSchemaAttribute::attribute_w_default( "presample_by_swa", xsct_rosetta_bool, "Presample the residue being added using SWA", "false" );
	attlist + XMLSchemaAttribute::attribute_w_default( "minimize_single_res", xsct_rosetta_bool, "Run in single residue minimization mode", "false" );
	attlist + XMLSchemaAttribute::attribute_w_default( "start_added_residue_in_aform", xsct_rosetta_bool, "Just start the added residue in A-form (for RNA)", "false" );
	attlist + XMLSchemaAttribute::attribute_w_default( "internal_cycles", xsct_non_negative_integer, "Number of presample monte carlo cycles", "50" );
	attlist + XMLSchemaAttribute::attribute_w_default( "sample_range_small", xsct_real, "For monte carlo presampling, small-size move size (degrees)", "5.0" );
	attlist + XMLSchemaAttribute::attribute_w_default( "sample_range_large", xsct_real, "For monte carlo presampling, large-size move size (degrees)", "40.0" );
	attlist + XMLSchemaAttribute::attribute_w_default( "kT", xsct_real, "Simulation kT", "0.5" );
	attlist + XMLSchemaAttribute::required_attribute( "swa_move", xs_string, "String describing the SWA move to perform" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}


//////////////////////////////////////////////////////////////////////
void
AddMover::apply( core::pose::Pose & viewer_pose, Size const res_to_add_in_full_model_numbering, Size const res_to_build_off_in_full_model_numbering )
{
	int offset = static_cast<int>(res_to_add_in_full_model_numbering) - static_cast<int>(res_to_build_off_in_full_model_numbering);
	runtime_assert( offset != 0 );
	AttachmentType attachment_type = (std::abs( offset ) == 1) ?
		(offset > 0 ? BOND_TO_PREVIOUS : BOND_TO_NEXT) :
		(offset > 0 ? JUMP_TO_PREV_IN_CHAIN : JUMP_TO_NEXT_IN_CHAIN );
	if ( !check_same_chain( viewer_pose, res_to_add_in_full_model_numbering, res_to_build_off_in_full_model_numbering ) ) attachment_type = JUMP_DOCK;
	StepWiseMove swa_move( res_to_add_in_full_model_numbering,
		Attachment( res_to_build_off_in_full_model_numbering, attachment_type ), ADD );
	apply( viewer_pose, swa_move );
}

//////////////////////////////////////////////////////////////////////
void
AddMover::apply(
#ifdef GL_GRAPHICS
	core::pose::Pose & viewer_pose,
#else
	core::pose::Pose & pose,
#endif
	StepWiseMove const & swa_move
) {
	using namespace core::pose;
	using namespace core::pose::full_model_info;

#ifdef GL_GRAPHICS
	Pose pose = viewer_pose; // hard copy -- try to avoid graphics problems when variants are added.
#endif
	runtime_assert( pose.size() >= 1 );

	// will record which new dofs added.
	suite_num_ = 0;
	nucleoside_num_ = 0;

	swa_move_ = swa_move;
	res_to_add_in_full_model_numbering_   = get_add_res( swa_move, pose );
	res_to_build_off_in_full_model_numbering_ = swa_move.attached_res();
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );

	runtime_assert( res_list.has_value( res_to_build_off_in_full_model_numbering_ ) );
	runtime_assert( !res_list.has_value( res_to_add_in_full_model_numbering_ ) );
	if ( swa_move.attachment_type() == JUMP_DOCK ) { runtime_assert( !check_same_chain( pose ) );
	} else { runtime_assert( check_same_chain( pose ) ); }

	// need to encapsulate the following residue addition & domain addition functions...
	if ( res_to_add_in_full_model_numbering_ > res_to_build_off_in_full_model_numbering_ ) {
		do_append(  pose );
	} else {
		do_prepend( pose );
	}
	if ( pose.size() == 2 && pose.fold_tree().num_jump() == 0 ) pose.fold_tree( kinematics::FoldTree( 2 ) ); //special 'from scratch' case.

	fix_up_jump_atoms_and_residue_type_variants( pose );

#ifdef GL_GRAPHICS
	viewer_pose = pose;
#else
	# define viewer_pose pose
#endif
	setup_initial_torsions( viewer_pose );

	///////////////////////////////////
	// Premodeler added residue
	///////////////////////////////////
	if ( presample_added_residue_ ) {
		if ( presample_by_swa_ ) {
			Size swa_sample_res = nucleoside_num_;
			if ( swa_sample_res == 0 ) swa_sample_res = suite_num_; // nucleoside_num_ = 0 in domain addition.
			sample_by_swa( viewer_pose, swa_sample_res );
		} else {
			sample_by_monte_carlo_internal( viewer_pose );
		}
	}
#ifndef GL_GRAPHICS
	# undef viewer_pose
#endif
}


////////////////////////////////////////////////////////////////////////
void
AddMover::do_append( core::pose::Pose & pose ){

	FullModelInfo & full_model_info = nonconst_full_model_info( pose );
	std::string const & full_sequence  = full_model_info.full_sequence();
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
	Size const res_to_build_off = res_list.index( res_to_build_off_in_full_model_numbering_ );

	Size const offset = res_to_add_in_full_model_numbering_ - res_to_build_off_in_full_model_numbering_;
	// addition to strand ending (append)
	if ( swa_move_.attachment_type() != JUMP_DOCK ) {
		runtime_assert( res_to_build_off == pose.size() ||
			res_list[ res_to_build_off ] < res_list[ res_to_build_off + 1 ] - 1 );
		runtime_assert( res_list[ res_to_build_off ] < full_sequence.size() );
	}

	Size const other_pose_idx = full_model_info.get_idx_for_other_pose_with_residue( res_to_add_in_full_model_numbering_ );
	if ( other_pose_idx ) { // addition of a domain (a whole sister pose)
		append_other_pose( pose, offset, other_pose_idx );
	} else { // single residue addition -- can this just be combine with above?
		append_residue( pose, offset );
	}
}

////////////////////////////////////////////////////////////////////////
void
AddMover::do_prepend( core::pose::Pose & pose ){

	FullModelInfo & full_model_info = nonconst_full_model_info( pose );
	runtime_assert( res_to_add_in_full_model_numbering_ < res_to_build_off_in_full_model_numbering_ );
	Size const offset = res_to_build_off_in_full_model_numbering_ - res_to_add_in_full_model_numbering_;
	Size const other_pose_idx = full_model_info.get_idx_for_other_pose_with_residue( res_to_add_in_full_model_numbering_ );

	//  TR << "About to add onto " << res_to_build_off_in_full_model_numbering_ << " the following residue (in full model numbering) " << res_to_add_in_full_model_numbering_ << " which may be part of other pose " << other_pose_idx << std::endl;

	if ( other_pose_idx ) { // addition of a domain (a whole sister pose)
		prepend_other_pose( pose, offset, other_pose_idx );
	} else {  // single residue addition -- can this just be combine with above?
		prepend_residue( pose, offset );
	}
}

////////////////////////////////////////////////////////////////////////
void
AddMover::append_other_pose( pose::Pose & pose, Size const offset,
	Size const other_pose_idx ){
	runtime_assert( other_pose_idx > 0);
	FullModelInfo & full_model_info = nonconst_full_model_info( pose );
	Pose & other_pose = *(full_model_info.other_pose_list()[ other_pose_idx ]);
	if ( swa_move_.attachment_type() == BOND_TO_PREVIOUS ) {
		runtime_assert( offset == 1 );
		merge_in_other_pose_by_bond( pose, other_pose, res_to_build_off_in_full_model_numbering_ );
	} else {
		runtime_assert( swa_move_.is_jump() );
		merge_in_other_pose_by_jump( pose, other_pose,
			res_to_build_off_in_full_model_numbering_, res_to_add_in_full_model_numbering_ );
	}
	nonconst_full_model_info( pose ).remove_other_pose_at_idx( other_pose_idx );

	suite_num_ = 0;
	nucleoside_num_ = get_res_list_from_full_model_info( pose ).index( res_to_add_in_full_model_numbering_ );
}

////////////////////////////////////////////////////////////////////////
void
AddMover::prepend_other_pose( pose::Pose & pose, Size const offset,
	Size const other_pose_idx ){
	runtime_assert( other_pose_idx > 0);
	FullModelInfo & full_model_info = nonconst_full_model_info( pose );
	Pose & other_pose = *(full_model_info.other_pose_list()[ other_pose_idx ]);
	if ( swa_move_.attachment_type() == BOND_TO_NEXT ) {
		runtime_assert( offset == 1 );
		merge_in_other_pose_by_bond( pose, other_pose, res_to_build_off_in_full_model_numbering_ - 1 /*merge_res*/ );
	} else {
		runtime_assert( swa_move_.is_jump() );
		merge_in_other_pose_by_jump( pose, other_pose,
			res_to_add_in_full_model_numbering_, res_to_build_off_in_full_model_numbering_  );
	}
	nonconst_full_model_info( pose ).remove_other_pose_at_idx( other_pose_idx );
	suite_num_ = 0;
	nucleoside_num_ = get_res_list_from_full_model_info( pose ).index( res_to_add_in_full_model_numbering_ );
}

////////////////////////////////////////////////////////////////////////
void
AddMover::append_residue( pose::Pose & pose, Size const offset ){

	using namespace core::chemical;
	using namespace core::chemical::rna;
	using namespace protocols::stepwise::modeler::rna;

	//FullModelInfo & full_model_info = nonconst_full_model_info( pose );
	//std::string const & full_sequence  = full_model_info.full_sequence();
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
	Size const res_to_build_off = res_list.index( res_to_build_off_in_full_model_numbering_ );
	Size res_to_add = res_to_build_off + 1;

	core::conformation::ResidueOP new_rsd = create_residue_to_add( pose );

	Size actual_offset = offset;
	if ( swa_move_.attachment_type() == BOND_TO_PREVIOUS ) {
		runtime_assert( offset == 1 );
		remove_variant_type_from_pose_residue( pose, core::chemical::UPPER_TERMINUS_VARIANT, res_to_build_off ); // got to be safe.
		remove_variant_type_from_pose_residue( pose, core::chemical::THREE_PRIME_PHOSPHATE, res_to_build_off ); // got to be safe.
		remove_variant_type_from_pose_residue( pose, core::chemical::C_METHYLAMIDATION, res_to_build_off ); // got to be safe.
		//   fix_protein_jump_atom( pose, res_to_build_off, " N  " );
		pose.append_polymer_residue_after_seqpos( *new_rsd, res_to_build_off, true /*build ideal geometry*/ );
		suite_num_ = res_to_add - 1;
	} else {
		runtime_assert( swa_move_.is_jump() );

		// following is getting really clunky. May want to instead use merge_two_poses code.
		Size takeoff_res( 0 ), k( 0 );
		for ( k = res_to_add_in_full_model_numbering_; k >= res_to_build_off_in_full_model_numbering_; k-- ) {
			if ( res_list.has_value( k ) ) { takeoff_res = res_list.index( k ); break; }
		}
		actual_offset = res_to_add_in_full_model_numbering_ - k;
		res_to_add = takeoff_res + 1;

		if ( swa_move_.attachment_type() == JUMP_DOCK ) {
			takeoff_res = res_to_build_off;
		}

		pose.insert_residue_by_jump( *new_rsd, res_to_add, takeoff_res,
			default_jump_atom( pose.residue_type( takeoff_res ) ),
			default_jump_atom( new_rsd->type() ) );

		kinematics::FoldTree f = pose.fold_tree();
		Size const jump_nr = f.jump_nr( takeoff_res, res_to_add );
		f.slide_jump( jump_nr, res_to_build_off, res_to_add );
		f.set_jump_atoms( jump_nr, default_jump_atom( pose.residue_type( res_to_build_off ) ),
			default_jump_atom( pose.residue_type( res_to_add ) ) );
		pose.fold_tree( f );

		if ( pose.residue_type( res_to_add ).is_RNA() ) {
			add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, res_to_add );
			if ( res_to_add_in_full_model_numbering_ < res_list[ res_to_build_off+1 ]-1 ) {
				add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_RIBOSE, res_to_add );
			}
		}
		suite_num_ = 0;
	}

	reorder_full_model_info_after_append( pose, res_to_add, actual_offset );

	nucleoside_num_ = res_to_add;
}


////////////////////////////////////////////////////////////////////////
void
AddMover::prepend_residue( pose::Pose & pose, Size const offset ){

	using namespace core::chemical;
	using namespace core::chemical::rna;
	using namespace protocols::stepwise::modeler::rna;

	//FullModelInfo & full_model_info = nonconst_full_model_info( pose );
	//std::string const & full_sequence  = full_model_info.full_sequence();
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info( pose );
	Size const res_to_build_off = res_list.index( res_to_build_off_in_full_model_numbering_ );
	Size res_to_add = res_to_build_off;

	runtime_assert( res_list[ res_to_add ] > 1 );

	core::conformation::ResidueOP new_rsd = create_residue_to_add( pose );

	Size actual_offset = offset;
	if ( swa_move_.attachment_type() == BOND_TO_NEXT ) {
		runtime_assert( offset == 1 );
		remove_variant_type_from_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE,res_to_build_off ); // got to be safe.
		remove_variant_type_from_pose_residue( pose, core::chemical::FIVE_PRIME_PHOSPHATE, res_to_build_off ); // got to be safe.
		remove_variant_type_from_pose_residue( pose, core::chemical::LOWER_TERMINUS_VARIANT,   res_to_build_off ); // got to be safe.
		remove_variant_type_from_pose_residue( pose, core::chemical::N_ACETYLATION,res_to_build_off ); // got to be safe.
		pose.prepend_polymer_residue_before_seqpos( *new_rsd, res_to_add, true /*build ideal geometry*/ );
		suite_num_ = res_to_add;
	} else {
		runtime_assert( swa_move_.is_jump() );

		// following is getting really clunky. May want to instead use merge_two_poses code.
		Size takeoff_res( 0 ), k( 0 );
		for ( k = res_to_add_in_full_model_numbering_; k <= res_to_build_off_in_full_model_numbering_; k++ ) {
			if ( res_list.has_value( k ) ) { takeoff_res = res_list.index( k ); break; }
		}
		actual_offset = k - res_to_add_in_full_model_numbering_;
		res_to_add = takeoff_res;

		// not sure how to match logic in append_residue yet ..
		//if ( swa_move_.attachment_type() == JUMP_DOCK ) {
		// takeoff_res = res_to_build_off;
		//}

		pose.insert_residue_by_jump( *new_rsd, res_to_add, takeoff_res,
			default_jump_atom( pose.residue_type( takeoff_res ) ),
			default_jump_atom( new_rsd->type() ) );

		kinematics::FoldTree f = pose.fold_tree();
		Size const jump_nr = f.jump_nr( res_to_add, takeoff_res+1 );
		f.slide_jump( jump_nr, res_to_add, res_to_build_off+1 );
		f.set_jump_atoms( jump_nr, default_jump_atom( pose.residue_type( res_to_add ) ),
			default_jump_atom( pose.residue_type( res_to_build_off+1 ) ) );
		pose.fold_tree( f );

		if ( res_to_add_in_full_model_numbering_ > res_list[ res_to_build_off-1 ]+1 ) {
			if ( pose.residue_type( res_to_add ).is_RNA() ) {
				add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_RIBOSE, res_to_add );
				if ( !pose.residue_type( res_to_add+1).has_variant_type( core::chemical::FIVE_PRIME_PHOSPHATE ) ) {
					add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_PHOSPHATE, res_to_add+1 );
				}
			}
		}
		suite_num_ = 0;
	}
	reorder_full_model_info_after_prepend( pose, res_to_add, actual_offset );

	// initialize with a random torsion... ( how about an A-form + perturbation ... or go to a 'reasonable' rotamer)
	nucleoside_num_ = res_to_add;
}

////////////////////////////////////////////////////////////////////
void
AddMover::sample_by_swa( pose::Pose & pose, Size const res_to_add ) const{
	using namespace core::pose::full_model_info;
	runtime_assert( stepwise_modeler_ != 0 );
	stepwise_modeler_->set_moving_res_and_reset( res_to_add );
	if ( !minimize_single_res_ ) stepwise_modeler_->set_working_minimize_res( get_moving_res_from_full_model_info( pose ) );

	stepwise_modeler_->apply( pose );
}


////////////////////////////////////////////////////////////////////
void
AddMover::sample_by_monte_carlo_internal( pose::Pose & pose ) const {

	using namespace protocols::moves;

	TR.Debug << "premodeler added residue! " << nucleoside_num_ << " over " << internal_cycles_ << " cycles " << std::endl;
	MonteCarloOP monte_carlo_internal( new MonteCarlo( pose, *scorefxn_, kT_ ) );

	std::string move_type( "" );
	for ( Size count_internal = 1; count_internal <= internal_cycles_; ++count_internal ) {

		rna_torsion_mover_->sample_near_suite_torsion( pose, suite_num_, sample_range_large_);
		if ( nucleoside_num_ > 0 ) rna_torsion_mover_->sample_near_nucleoside_torsion( pose, nucleoside_num_, sample_range_large_);
		monte_carlo_internal->boltzmann( pose, move_type );

		rna_torsion_mover_->sample_near_suite_torsion( pose, suite_num_, sample_range_small_);
		if ( nucleoside_num_ > 0 ) rna_torsion_mover_->sample_near_nucleoside_torsion( pose, nucleoside_num_, sample_range_small_);
		monte_carlo_internal->boltzmann( pose, move_type );
		//std::cout << "During premodeler: " << (*scorefxn_)( pose );
	} // monte carlo cycles
}

///////////////////////////////////////////////////////////////////
void
AddMover::set_stepwise_modeler( protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler ){
	stepwise_modeler_ = stepwise_modeler;
}

///////////////////////////////////////////////////////////////////
bool
AddMover::check_same_chain( pose::Pose const & pose, Size const res_to_add_in_full_model_numbering, Size const res_to_build_off_in_full_model_numbering ){
	return  ( get_chain_for_full_model_resnum( res_to_add_in_full_model_numbering, pose ) ==
		get_chain_for_full_model_resnum( res_to_build_off_in_full_model_numbering, pose ) );
}

///////////////////////////////////////////////////////////////
bool
AddMover::check_same_chain( pose::Pose const & pose ) {
	return check_same_chain( pose, res_to_add_in_full_model_numbering_, res_to_build_off_in_full_model_numbering_ );
}

///////////////////////////////////////////////////////////////
void
AddMover::setup_initial_torsions( pose::Pose & pose ) {

	if ( suite_num_ > 0 && pose.residue_type( suite_num_ ).is_RNA() ) {
		if ( start_added_residue_in_aform_ ) {
			rna_torsion_mover_->apply_suite_torsion_Aform( pose, suite_num_ );
			if ( nucleoside_num_ > 0 ) rna_torsion_mover_->apply_nucleoside_torsion_Aform( pose, nucleoside_num_ );
		} else {
			rna_torsion_mover_->apply_random_suite_torsion( pose, suite_num_ );
			if ( nucleoside_num_ > 0 ) rna_torsion_mover_->apply_random_nucleoside_torsion( pose, nucleoside_num_ );
		}

		// rna_torsion_mover_->sample_near_suite_torsion( pose, suite_num_, sample_range_large_);
		// if ( nucleoside_num_ > 0 ) rna_torsion_mover_->sample_near_nucleoside_torsion( pose, nucleoside_num_, sample_range_large_);
	} else if ( suite_num_ == 0 && swa_move_.is_jump() ) {
		setup_initial_jump( pose );
	}
}

// similar to get_remodel_res in ResampleMover but only looks inside full_model_info.
Size
AddMover::get_add_res( StepWiseMove const & swa_move, pose::Pose const & pose ) const {

	utility::vector1< Size > const & res_list = const_full_model_info( pose ).res_list();
	utility::vector1< Size > const & cutpoint_open_in_full_model = const_full_model_info( pose ).cutpoint_open_in_full_model();
	std::string const & full_sequence = const_full_model_info( pose ).full_sequence();

	Size const & attached_res =  swa_move.attached_res();
	runtime_assert( res_list.has_value( attached_res ) );

	MoveElement const & move_element = swa_move.move_element();
	for ( Size n = 1; n <= move_element.size(); n++ ) {
		Size const i = move_element[ n ];
		runtime_assert(  !res_list.has_value( i ) );
	}

	AttachmentType const & attachment_type = swa_move.attachment_type();
	Size add_res( 0 );
	if ( attachment_type == BOND_TO_PREVIOUS ) {
		add_res = attached_res + 1;
		runtime_assert( !cutpoint_open_in_full_model.has_value( attached_res ) );
	} else if ( attachment_type == BOND_TO_NEXT ) {
		add_res = attached_res - 1;
		runtime_assert( !cutpoint_open_in_full_model.has_value( add_res ) );
	} else if ( attachment_type == JUMP_TO_PREV_IN_CHAIN ) {
		add_res = attached_res + 1;
		while ( add_res < full_sequence.size() && !move_element.has_value( add_res ) ) add_res++;
	} else if ( attachment_type == JUMP_TO_NEXT_IN_CHAIN ) {
		add_res = attached_res - 1;
		while ( add_res > 1 && !move_element.has_value( add_res ) ) add_res--;
	} else {
		runtime_assert( attachment_type == JUMP_DOCK );
		runtime_assert( move_element.size() == 1 );
		add_res = move_element[ 1 ];
	}
	runtime_assert( move_element.has_value( add_res ) );

	return add_res;
}

void
AddMover::setup_initial_jump( pose::Pose & pose ) {

	using namespace core::kinematics;
	using namespace core::pose::full_model_info;
	using namespace core::chemical::rna;

	// get the res to add in local numbering
	FullModelInfo const & full_model_info = nonconst_full_model_info( pose );
	Size const & res_to_add = full_model_info.full_to_sub( res_to_add_in_full_model_numbering_ );

	// get parent res of res added, make sure they are connected by a jump ... is this necessary?
	bool connected_by_jump;
	FoldTree ft( pose.fold_tree() );
	Size const & parent_res = ft.get_parent_residue( res_to_add, connected_by_jump );
	int offset = static_cast<int>(parent_res) - static_cast<int>(res_to_add);
	if ( std::abs(offset) > 1 ) {
		runtime_assert( swa_move_.attachment_type() == JUMP_DOCK );
	}
	runtime_assert( connected_by_jump );

	// create a new jump, and initialize its translation
	Jump new_jump = Jump();
	numeric::xyzVector< Real > t( 0.0, 0.0, 7.0 );
	if ( swa_move_.attachment_type() == JUMP_DOCK  ) {
		t = 2.0;
	}
	new_jump.set_translation( t );

	// get jump number, and add the new jump into pose
	Size const & jump_nr = ft.get_jump_that_builds_residue( res_to_add );
	pose.set_jump( jump_nr, new_jump );
}


///////////////////////////////////////////////////////////////////////
/// @brief: figure out the res type of the residue to add and create it
/// @details: also handles the question of 'what if I have chemical modifications'
/// and 'what if I have protonation sampling turned on'. Probably should factor
/// out these concerns.
core::conformation::ResidueOP
AddMover::create_residue_to_add( pose::Pose const & pose ) {

	// figure out new residue type from sequence
	std::string const & full_sequence  = const_full_model_info( pose ).full_sequence();
	std::map< Size, std::string > const & nc_res_map = const_full_model_info( pose ).full_model_parameters()->non_standard_residue_map();

	if ( TR.Debug.visible() ) {
		for ( auto const & nc_res_pair : nc_res_map ) {
			TR.Debug << nc_res_pair.first << " " << nc_res_pair.second << std::endl;
		}
	}

	std::string newrestype3;
	if ( nc_res_map.find( res_to_add_in_full_model_numbering_ ) != nc_res_map.end() ) {
		newrestype3 = nc_res_map.at( res_to_add_in_full_model_numbering_ );
	} else {
		// Otherwise, we need to look at the fasta. Get the relevant character.
		char newrestype = full_sequence[ res_to_add_in_full_model_numbering_ - 1 ];
		if ( newrestype == 'n' ) {
			if ( ! designing_with_noncanonicals_ ) {
				// output parameter is also newrestype
				modeler::rna::choose_random_if_unspecified_nucleotide( newrestype );
				newrestype3 = chemical::name_from_aa( chemical::aa_from_oneletter_code( newrestype ) );
			} else {
				newrestype3 = modeler::rna::choose_randomly_from_allowed_at_position( pose, res_to_add_in_full_model_numbering_ );
			}
		} else {
			newrestype3 = chemical::name_from_aa( chemical::aa_from_oneletter_code( newrestype ) );
		}
	}

	TR.Debug << "Going to add a " << newrestype3 << " for " << res_to_add_in_full_model_numbering_ << "." << std::endl;
	TR.Debug << "Note: NSRM was " << std::endl;
	for ( auto const & elem : nc_res_map ) {
		TR.Debug << elem.first << " : " << elem.second << std::endl;
	}

	// figure out residue type from name3
	chemical::ResidueTypeSetCOP rsd_set = pose.residue_type_set_for_pose();

	// Handle pH sampling: if sample_pH_ is true, then we need to have a 50/50
	// shot of creating a protonated adenosine (50/50 because score will take
	// care of the rest)
	chemical::ResidueTypeCOP rsd_type = nullptr;
	if ( newrestype3 == "RAD" && sample_pH_ ) {
		TR << "Going to add adenosine." << std::endl;
		// There is something irritating about enforcing the exact equality of
		// [0, 0.5) and [0.5, 1) in floating point space, but we don't care for
		// now, certainly.
		TR << "Selected ";
		if ( numeric::random::rg().uniform() < 0.5 ) {
			TR << " PROTONATED ";
			newrestype3 += ":Protonated_N1_Adenosine";
			rsd_type = chemical::ResidueTypeCOP( rsd_set->name_mapOP( newrestype3 ) );
			if ( !rsd_type ) {
				utility_exit_with_message( "We couldn't generate a residue type of name " + newrestype3 + "; apologies." );
			}
		} else {
			TR << " NORMAL ";
			rsd_type = chemical::ResidueTypeCOP( rsd_set->get_representative_type_base_name( newrestype3 ) );
		}
		TR << " variant." << std::endl;
	} else {
		rsd_type = chemical::ResidueTypeCOP( rsd_set->get_representative_type_base_name( newrestype3 ) );
	}

	// Might not be able to get by base name if adding a ligand -- particularly
	// from the CCD?
	if ( !rsd_type ) {
		rsd_type = chemical::ResidueTypeCOP( rsd_set->name_map( newrestype3 ).get_self_ptr() );
	}
	return conformation::ResidueFactory::create_residue( *rsd_type );
}

std::string AddMoverCreator::keyname() const {
	return AddMover::mover_name();
}

protocols::moves::MoverOP
AddMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddMover );
}

void AddMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddMover::provide_xml_schema( xsd );
}



} //mover
} //monte_carlo
} //stepwise
} //protocols
