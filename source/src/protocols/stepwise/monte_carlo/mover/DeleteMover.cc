// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file DeleteMover
/// @brief Deletes an RNA residue from a chain terminus.
/// @details
/// @author Rhiju Das

#include <protocols/stepwise/monte_carlo/mover/DeleteMover.hh>
#include <protocols/stepwise/monte_carlo/mover/DeleteMoverCreator.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMoveSelector.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/modeler/packer/util.hh>
#include <protocols/stepwise/modeler/rna/phosphate/MultiPhosphateSampler.hh>
#include <protocols/stepwise/modeler/rna/phosphate/util.hh>
#include <protocols/stepwise/modeler/StepWiseModeler.hh>
#include <protocols/stepwise/modeler/util.hh>

// libRosetta headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <protocols/scoring/VDW_CachedRepScreenInfo.hh>
#include <protocols/moves/mover_schemas.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/tools/make_vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>

#include <basic/Tracer.hh>

//Req'd on WIN32

using namespace core;
using namespace core::pose::full_model_info;
using namespace protocols::stepwise::modeler;
using core::Real;

//////////////////////////////////////////////////////////////////////////
// Removes one residue from a 5' or 3' chain terminus, and appropriately
//  updates the pose full_model_info object.
//
// Now updated to slice off more than one residue, or even an entire
//  fixed domain, and save that slice in an 'other_pose' (cached in the
//  main pose inside full_model_info).
//
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.stepwise.monte_carlo.mover.DeleteMover" );

namespace protocols {
namespace stepwise {
namespace monte_carlo {
namespace mover {


//////////////////////////////////////////////////////////////////////////
//constructor!
DeleteMover::DeleteMover( ):
	options_( utility::pointer::make_shared< protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptions >() )
{}

//////////////////////////////////////////////////////////////////////////
//destructor
DeleteMover::~DeleteMover() = default;

protocols::moves::MoverOP
DeleteMover::clone() const {
	return utility::pointer::make_shared< DeleteMover >( *this );
}

void DeleteMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & // data
) {

	minimize_after_delete_ = tag->getOption< bool >( "minimize_after_delete", minimize_after_delete_ );

	move_str_ = tag->getOption< std::string >( "swa_move" );
}

void DeleteMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	//rosetta_scripts::attributes_for_parse_score_function( attlist );
	attlist + XMLSchemaAttribute::attribute_w_default( "minimize_after_delete", xsct_rosetta_bool, "Minimize after delete", "true" );
	attlist + XMLSchemaAttribute::required_attribute( "swa_move", xs_string, "String describing the SWA move to perform" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "XRW TO DO", attlist );
}


//////////////////////////////////////////////////////////////////////////
void
DeleteMover::apply( core::pose::Pose & pose )
{
	//std::cout << "not defined" << std::endl;
	// We get here from RosettaScripts applications
	runtime_assert( ! move_str_.empty() );
	StepWiseMove swa_move( utility::string_split( move_str_ ), const_full_model_info( pose ).full_model_parameters() );
	apply( pose, swa_move.move_element() );
}


//////////////////////////////////////////////////////////////////////
void
DeleteMover::apply( core::pose::Pose & pose, core::Size const res_to_delete_in_full_model_numbering )
{
	apply( pose, utility::tools::make_vector1( res_to_delete_in_full_model_numbering ) );
}


//////////////////////////////////////////////////////////////////////
void
DeleteMover::apply( core::pose::Pose & pose, utility::vector1< core::Size > const & residues_to_delete_in_full_model_numbering )
{
	using namespace core::pose;

	// in case the deletion is supposed to occur on a pose on which we are not focused.
	core::Size const idx = const_full_model_info( pose ).get_idx_for_other_pose_with_residue( residues_to_delete_in_full_model_numbering[ 1 ] );
	if ( idx > 0 ) switch_focus_to_other_pose( pose, idx );

	// have at it!
	FullModelInfo & full_model_info = nonconst_full_model_info( pose );
	utility::vector1< core::Size > const residues_to_delete = full_model_info.full_to_sub( residues_to_delete_in_full_model_numbering );
	interface_res_ = full_model_info.sub_to_full( packer::figure_out_working_interface_res( pose, get_unique_connection_res( pose, residues_to_delete ), options_->pack_protein_side_chains() ) );

	// do the slice.
	PoseOP sliced_out_pose_op( new Pose );
	slice_out_pose( pose, *sliced_out_pose_op, residues_to_delete );

	// get rid of pieces that are single residues. no need to hold on to those.
	bool keep_remainder_pose( true ), keep_sliced_out_pose( true );
	remove_singletons_and_update_pose_focus( pose, sliced_out_pose_op, keep_remainder_pose, keep_sliced_out_pose );

	if ( keep_remainder_pose  ) fix_up_jump_atoms_and_residue_type_variants( pose );
	if ( keep_sliced_out_pose ) fix_up_jump_atoms_and_residue_type_variants( *sliced_out_pose_op );

	if ( options_->minimize_after_delete() ) {
		if ( keep_remainder_pose  ) minimize_after_delete( pose );
		if ( keep_sliced_out_pose ) minimize_after_delete( *sliced_out_pose_op );
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
DeleteMover::decide_to_keep_pose( pose::Pose const & pose ) const {
	if ( const_full_model_info( pose ).is_a_submotif() ) return const_full_model_info( pose ).is_a_submotif_seed();
	if ( check_for_input_domain( pose ) ) return true;
	if ( ( (options_->from_scratch_frequency() > 0.0) || options_->allow_split_off() ) && pose.size() > 1 )  return true;
	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
DeleteMover::remove_singletons_and_update_pose_focus( core::pose::Pose & pose,
	core::pose::PoseOP sliced_out_pose_op,
	bool & keep_remainder_pose,
	bool & keep_sliced_out_pose ) const {

	FullModelInfo & full_model_info = nonconst_full_model_info( pose );

	// make decision on which poses to keep.
	keep_remainder_pose  = decide_to_keep_pose( pose );
	keep_sliced_out_pose = decide_to_keep_pose( *sliced_out_pose_op );

	if ( keep_sliced_out_pose ) full_model_info.add_other_pose( sliced_out_pose_op );
	if ( keep_remainder_pose ) return;
	// remainder pose is a single nucleotide. no need to keep track of it anymore.

	runtime_assert( full_model_info.res_list().size() > 0 );
	core::Size const res_in_remainder_pose = full_model_info.res_list()[ 1 ];
	bool pose_is_alone( false );
	if ( keep_sliced_out_pose ) { // go to the sliced out pose.
		core::Size const sliced_out_pose_idx = full_model_info.get_idx_for_other_pose( *sliced_out_pose_op );
		switch_focus_to_other_pose( pose, sliced_out_pose_idx );
	} else if ( full_model_info.other_pose_list().size() > 0 ) { // switch focus randomly.
		switch_focus_among_poses_randomly( pose, nullptr, true /*force_switch*/ );
	} else {
		pose_is_alone = true;
	}

	if ( pose_is_alone ) {
		TR << TR.Red << "EMPTY POSE! Keeping the pose with number of residues " << pose.size() << TR.Reset << std::endl;

		core::pose::PoseCOP old_pose_cop = pose.clone();
		FullModelInfoOP new_full_model_info = nonconst_full_model_info( pose ).clone_info();
		// Rosetta's Pose object craps out if it is blank.
		// Still want to have blank poses represent fully random chain in stepwise monte carlo,
		// and would be a good state to have, or at least transit through.
		// Creating a single virtual residue, with res_list cleared, acts as a very
		// special case, and rest of the SWA_MonteCarlo code checks for res_list.size().
		core::chemical::ResidueTypeCOP rsd_type( core::pose::virtual_type_for_pose(pose) );
		core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *rsd_type ) );
		pose.clear();
		pose.append_residue_by_bond( *new_res );

		new_full_model_info->clear_res_list();
		new_full_model_info->update_submotif_info_list();
		new_full_model_info->clear_other_pose_list();
		set_full_model_info( pose, new_full_model_info );
		protocols::scoring::set_vdw_cached_rep_screen_info_from_pose( pose, *old_pose_cop );
	} else {
		FullModelInfo & new_full_model_info = nonconst_full_model_info( pose );
		core::Size const remainder_pose_idx = new_full_model_info.get_idx_for_other_pose_with_residue( res_in_remainder_pose );
		new_full_model_info.remove_other_pose_at_idx( remainder_pose_idx );
		keep_remainder_pose = true; // will tell outside functions to clean up pose variants, etc.
	}
}


////////////////////////////////////////////////////////////////////
void
DeleteMover::minimize_after_delete( pose::Pose & pose ) const{

	using namespace core::pose::full_model_info;

	// normally happens with modeler -- this is important for seeing if terminal phosphates that
	// previously instantiated due to contact with deleted residues need to *disappear* now.
	if ( options_->sampler_perform_phosphate_pack() ) {
		using namespace protocols::stepwise::modeler::rna::phosphate;
		MultiPhosphateSampler phosphate_sampler( pose );
		phosphate_sampler.set_scorefxn( get_phosphate_scorefxn( stepwise_modeler_->scorefxn()->energy_method_options() ) );
		phosphate_sampler.sample_phosphates(); // samples on internal clone of pose.
		phosphate_sampler.copy_phosphates( pose );
	}

	stepwise_modeler_->set_moving_res_and_reset( 0 );
	stepwise_modeler_->set_working_prepack_res( full_to_sub( interface_res_, pose ) );
	stepwise_modeler_->set_working_minimize_res( get_moving_res_from_full_model_info( pose ) );
	stepwise_modeler_->apply( pose );
}

//////////////////////////////////////////////////////////////////////
// following is not really in use anymore -- perhaps should deprecate.
void
DeleteMover::wipe_out_moving_residues( pose::Pose & pose ) {

	// don't do any minimizing -- just get rid of everything...
	bool const minimize_after_delete_save( minimize_after_delete_ );
	minimize_after_delete_ = false;

	utility::vector1< StepWiseMove > swa_moves;
	StepWiseMoveSelector swa_move_selector;
	swa_move_selector.get_intramolecular_delete_move_elements( pose, swa_moves);

	if ( swa_moves.size() > 0 ) { // recursively delete all residues.
		apply( pose, swa_moves[1].move_element() );
		wipe_out_moving_residues( pose );
	}

	minimize_after_delete_ = minimize_after_delete_save;
}

///////////////////////////////////////////////////////////////////
void
DeleteMover::set_stepwise_modeler( protocols::stepwise::modeler::StepWiseModelerOP stepwise_modeler ){
	stepwise_modeler_ = stepwise_modeler;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
DeleteMover::set_options( protocols::stepwise::monte_carlo::options::StepWiseMonteCarloOptionsCOP options ){
	options_ = options;
}

//////////////////////////////////////////////////////////////////////
std::string
DeleteMover::get_name() const {
	return mover_name();
}

std::string
DeleteMover::mover_name() {
	return "DeleteMover";
}

std::string DeleteMoverCreator::keyname() const {
	return DeleteMover::mover_name();
}

protocols::moves::MoverOP
DeleteMoverCreator::create_mover() const {
	return utility::pointer::make_shared< DeleteMover >();
}

void DeleteMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	DeleteMover::provide_xml_schema( xsd );
}


} //mover
} //monte_carlo
} //stepwise
} //protocols
