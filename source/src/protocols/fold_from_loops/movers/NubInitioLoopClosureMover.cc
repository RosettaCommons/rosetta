// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/NubInitioLoopClosureMover.cc
/// @brief Applies CCD closure to a NubInitioMover results being aware of its restrictions through labels.
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#include <protocols/fold_from_loops/movers/NubInitioLoopClosureMover.hh>
#include <protocols/fold_from_loops/movers/NubInitioLoopClosureMoverCreator.hh>
#include <protocols/fold_from_loops/selectors/CutpointResidueSelector.hh>
#include <protocols/fold_from_loops/movers/LabelPoseFromResidueSelectorMover.hh>
#include <protocols/fold_from_loops/movers/DisplayPoseLabelsMover.hh>

// Protocol headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/loops/loop_closure/ccd/SlidingWindowLoopClosure.hh>
#include <protocols/loops/Exceptions.hh>
#include <protocols/jumping/util.hh>
#include <protocols/checkpoint/CheckPointer.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResConnID.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/AndResidueSelector.hh>
#include <core/select/residue_selector/OrResidueSelector.hh>
#include <core/select/residue_selector/PrimarySequenceNeighborhoodSelector.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/ResLvlTaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.fold_from_loops.NubInitioLoopClosureMover" );

namespace protocols {
namespace fold_from_loops {
namespace movers {

NubInitioLoopClosureMover::NubInitioLoopClosureMover():
	protocols::moves::Mover( mover_name() ),
	centroid_scorefxn_( default_centroid_scorefxn() ),
	fullatom_scorefxn_( default_fullatom_scorefxn() ),
	break_side_( default_break_side() ),
	break_side_ramp_value_( 0 ),
	break_side_ramp_( default_break_side_ramp() ),
	max_trials_( default_max_trials() ),
	fragments_( new core::fragment::ConstantLengthFragSet ),
	altered_residues_( core::select::residue_selector::ReturnResidueSubsetSelector() ),
	trust_( default_trust() ),
	label_( default_label() ),
	centroid_( default_centroid() )
{}

NubInitioLoopClosureMover::~NubInitioLoopClosureMover()= default;

void
NubInitioLoopClosureMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	if ( tag->hasOption( "centroid_scorefxn" ) ) {
		centroid_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, "centroid_scorefxn", data ) );
	}
	if ( tag->hasOption( "fullatom_scorefxn" ) ) {
		fullatom_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, "fullatom_scorefxn", data ) );
	}
	max_trials( tag->getOption< core::Size >( "max_trials", default_max_trials() ) );
	break_side( tag->getOption< core::Size >( "break_side", default_break_side() ) );
	break_side_ramp( tag->getOption< bool >( "break_side_ramp", default_break_side_ramp() ) );
	trust( tag->getOption< bool >( "trust", default_trust() ) );
	label( tag->getOption< bool >( "label", default_label() ) );
	design( tag->getOption< bool >( "design", default_design() ) );

	fragments_id( tag->getOption< std::string >( "fragments_id", std::string() ) );
	if ( data.has( "FragSet", fragments_id_ + "small" ) ) {
		fragments_ = data.get_ptr<core::fragment::FragSet>( "FragSet", fragments_id_ + "small");

		TR.Trace << TR.Green << "Fragments loaded through the basic::datacache::DataMap" << TR.Reset << std::endl;
		TR.Trace << TR.Green << "Small: " << data.get_ptr<core::fragment::FragSet>( "FragSet", fragments_id_ + "small") << TR.Reset << std::endl;
	}
}

void NubInitioLoopClosureMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "max_trials", xs_integer,
		"Defines how many times we should try to close the breakpoints before giving up.", std::to_string( default_max_trials() ))
		+ XMLSchemaAttribute::required_attribute( "fragments_id", xs_string,
		"Fragments are necessary to close the loop. They need to be included/created with StructFragmentMover. "
		"The value set in the 'prefix' attribute of that Mover needs to be provided here again.")
		+ XMLSchemaAttribute::attribute_w_default( "break_side", xs_integer,
		"Defines number of sequence neighbors around breakpoint residues (multisegment motif) allowed to minimize in order to achieve the chain break closure.",
		std::to_string( default_break_side() ))
		+ XMLSchemaAttribute::attribute_w_default( "break_side_ramp", xsct_rosetta_bool,
		"In each breakpoint closure trial, increase by one the number of side residues allowed to move.",
		std::to_string( default_break_side_ramp() ))
		+ XMLSchemaAttribute::attribute_w_default( "trust", xsct_rosetta_bool,
		"Do nothing on those cutpoints that seem to have been closed before reaching this Mover.",
		std::to_string( default_trust() ))
		+ XMLSchemaAttribute::attribute_w_default( "label", xsct_rosetta_bool,
		"If set to true (default), labels residues changed by loop_closure to LOOPCLOSURE and repacked residues as REPACKED.",
		std::to_string( default_label() ))
		+ XMLSchemaAttribute::attribute_w_default( "design", xsct_rosetta_bool,
		"If set to true (not default), redesigns residues of the closed loops and around as longs as they are labeled as TEMPLATE or COLDSPOT.",
		std::to_string( default_design() ));
	protocols::rosetta_scripts::attributes_for_parse_score_function_w_description( attlist, "centroid_scorefxn", "Centroid Score function to use." );
	protocols::rosetta_scripts::attributes_for_parse_score_function_w_description( attlist, "fullatom_scorefxn", "Centroid Score function to use." );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Applies CCD closure to a NubInitioMover result, being aware of its restrictions through labels.", attlist );
}

void
NubInitioLoopClosureMover::apply( core::pose::Pose & pose )
{
	using namespace core::kinematics;
	using namespace core::select::residue_selector;
	using namespace protocols::loops::loop_closure::ccd;

	// Check if the pose has the expected labels to be closed; otherwise there is nothing that
	// can be done here.
	if ( !has_movable_residues( pose ) ) {
		set_last_move_status( moves::FAIL_BAD_INPUT );
		TR << TR.Bold << TR.Red << "[ERROR]" << TR.Reset << std::endl;
		TR << get_name() << " is prepared to work with NubInitioMover-derived Poses." << std::endl;
		TR << "It only allows to move residues labeled as TEMPLATE or FLEXIBLE. Thus, without residues ";
		TR << "labeled as such, it will never succeed in closing the loops." << std::endl;
		TR << "If you want to use this Mover you must add one of those labels to your residues.";
		TR << " Otherwise, you should use some other loop closure Mover." << std::endl;
		return;
	}

	// If the Pose has been tried to be
	if ( trust_ ) {
		TR << "Trusting previous closure modifications: adapting FoldTree to Residue Variants" << std::endl;
		make_foldtree_coherent_to_cutpoints( pose );
	} else {
		TR << "NOT trusting previous closure modifications: adapting Residue Variants to  FoldTree" << std::endl;
		make_cutpoints_coherent_to_foldtree( pose );
	}

	// Check if there is any cutpoint, otherwise just skip alltogether
	if ( !has_cutpoints( pose ) ) {
		TR << TR.Bold << TR.Green << "[WARNING]" << TR.Reset << std::endl;
		TR << "No cutpoints detected. " << get_name() << " has nothing to do here." << std::endl;
		core::pose::setPoseExtraScore( pose, "niccd_break_side", 0 );
		core::pose::setPoseExtraScore( pose, "niccd_trials", -1 );
		set_last_move_status( moves::MS_SUCCESS );
		return;
	}

	// Initial expected work:
	show_cutpoints( pose );

	// Setup Centroid Score to guide loop closure
	core::scoring::ScoreFunctionOP closure_score = set_centroid_default_chainclosure_weight( centroid_scorefxn_ );
	TR << std::endl << "Closing loops with Centroid Score:" << std::endl;
	closure_score->show_pretty( TR );

	// Setup Closing Conditions
	FoldTree final_tree = make_final_tree( pose );
	FoldTree init_tree  = pose.fold_tree();
	// apply_closure_trust( pose );

	// If we have trusted previous closure procedures, we might not need to do anything
	// Check again if there is anything to do.
	if ( !has_cutpoints( pose ) ) {
		TR << "Cutpoints were already closed" << std::endl;
		core::pose::setPoseExtraScore( pose, "niccd_break_side", 0 );
		core::pose::setPoseExtraScore( pose, "niccd_trials", 0 );
		pose.fold_tree( final_tree );
		set_last_move_status( moves::MS_SUCCESS );
		return;
	}

	// Save constraints: we are going to take them out, as sometimes they complain
	// (i think this happens on the full atom to centroin conversion)
	core::scoring::constraints::ConstraintSetOP saved_cst( new core::scoring::constraints::ConstraintSet() );
	saved_cst->detached_copy( *pose.constraint_set() );
	pose.remove_constraints();
	TR.Trace << "constraints saved" << std::endl;

	// Setup Loop Closure
	bool success = false;
	core::Size trials = 1;
	SlidingWindowLoopClosureOP closure_protocol( new SlidingWindowLoopClosure );
	closure_protocol->scorefxn( closure_score );
	closure_protocol->scored_frag_cycle_ratio( 0.2 );
	closure_protocol->short_frag_cycle_ratio( 0.1 );
	closure_protocol->set_bIdealLoopClosing( false );
	closure_protocol->set_chainbreak_max( 2 );
	closure_protocol->fragments( fragments_ );

	// Save sidechains a go back to centroid
	core::pose::Pose saved_fa_pose = pose;
	simple_moves::SwitchResidueTypeSetMoverOP switcher( new simple_moves::SwitchResidueTypeSetMover );
	switcher->type_set_tag( core::chemical::CENTROID );
	switcher->apply( pose );

	// Try to close
	while ( !success and trials < max_trials_ ) {
		closure_protocol->movemap( loop_closure_movemap( pose ) );
		// Setup CheckPoint
		protocols::checkpoint::CheckPointer sliding_checkpoint("closing");
		// Run
		try {
			if ( label_ ) {
				LabelPoseFromResidueSelectorMover labeler( active_movable_residues( pose ), "LOOPCLOSURE" );
				labeler.apply( pose );
			}
			jumping::close_chainbreaks( closure_protocol, pose, sliding_checkpoint , "sliding", final_tree );
			success = true;
		} catch ( loops::EXCN_Loop_not_closed& excn ) {
			if ( break_side_ramp_ ) {
				++break_side_ramp_value_;
			}
			++trials;
		}
	}
	// If we cannot close, let's just leave
	if ( !success ) {
		set_last_move_status( moves::FAIL_DO_NOT_RETRY );
		break_side_ramp_value_ = 0;
		return;
	}
	if ( centroid_ ) {
		set_last_move_status( moves::MS_SUCCESS );
		break_side_ramp_value_ = 0;
		return;
	}
	switcher->type_set_tag( core::chemical::FA_STANDARD );
	switcher->apply( pose );
	for ( core::Size i=1; i<=saved_fa_pose.size(); ++i ) {
		pose.replace_residue( i, saved_fa_pose.residue( i ), true );
		if ( pose.conformation().residue( i ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) {
			TR.Trace << "Removing LOWER VARIANT on residue " << i << std::endl;
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_LOWER, i );
		}
		if ( pose.conformation().residue( i ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
			TR.Trace << "Removing UPPER VARIANT on residue " << i << std::endl;
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_UPPER, i );
		}
	}
	core::pose::setPoseExtraScore( pose, "niccd_break_side", break_side_ + break_side_ramp_value_ );
	core::pose::setPoseExtraScore( pose, "niccd_trials", trials );
	pose.fold_tree( final_tree );
	pose.update_residue_neighbors();
	repack( pose );
	if ( label_ ) {
		DisplayPoseLabelsMover shower;
		shower.write( true );
		shower.apply( pose );
	}
	break_side_ramp_value_ = 0;
	if ( saved_cst->has_constraints() ) {
		pose.constraint_set( saved_cst );
	}
	set_last_move_status( moves::MS_SUCCESS );
}

void
NubInitioLoopClosureMover::apply_closure_trust( core::pose::Pose & pose, bool by_variant ) const {
	using namespace core::kinematics;
	FoldTree work_tree  = pose.fold_tree();
	if ( trust_ ) {
		for ( core::Size i = 1; i < pose.size(); i++ ) {
			if ( pose.conformation().residue( i ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) {
				if ( pose.residue( i ).chain() == pose.residue( i + 1 ).chain() ) {
					if ( is_cutpoint_close( pose, i, by_variant ) ) {
						core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_LOWER, i );
						core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_UPPER, i + 1 );
						core::Size edge_res1 = work_tree.boundary_left( i );
						if ( edge_res1 == i ) edge_res1 = work_tree.boundary_right( i );
						Edge edge1 = work_tree.get_residue_edge( i );
						core::Size edge_res2 = work_tree.boundary_left( i + 1 );
						if ( edge_res2 == i + 1 ) edge_res2 = work_tree.boundary_right( i + 1 );
						Edge edge2 = work_tree.get_residue_edge( i + 1 );
						work_tree.delete_edge( edge1 );
						work_tree.delete_edge( edge2 );
						work_tree.add_edge( edge_res1, edge_res2, -1 );
						core::Size jump = work_tree.jump_nr( work_tree.root(), edge_res2 );
						Edge jumpe = work_tree.jump_edge( jump );
						work_tree.delete_edge( jumpe );
						work_tree.renumber_jumps();
						work_tree.delete_extra_vertices();
					}
				}
			}
		}
		pose.fold_tree( work_tree );
		TR << "Expected after accept trusted closed cutpoints:" << std::endl;
		show_cutpoints( pose );
	}
}

core::kinematics::MoveMapOP
NubInitioLoopClosureMover::loop_closure_movemap( core::pose::Pose const & pose ) {
	core::select::movemap::MoveMapFactoryOP closure_mapFCTRY( new core::select::movemap::MoveMapFactory );
	closure_mapFCTRY->all_bb( false );
	closure_mapFCTRY->all_chi( false );
	closure_mapFCTRY->add_chi_action( core::select::movemap::mm_enable, active_movable_residues( pose ) );
	closure_mapFCTRY->add_bb_action( core::select::movemap::mm_enable, active_movable_residues( pose ) );
	core::kinematics::MoveMapOP closure_map = closure_mapFCTRY->create_movemap_from_pose( pose );
	return closure_map;
}

bool
NubInitioLoopClosureMover::is_cutpoint_close( core::pose::Pose const & pose, core::Size const & residue, bool by_variant ) const {
	if ( by_variant ) {
		return !pose.residue(residue).has_variant_type( core::chemical::CUTPOINT_LOWER );
	}
	core::Real const CNdist( pose.residue( residue ).xyz( "C" ).distance( pose.residue( residue + 1 ).xyz( "N" ) ) );
	core::Real const NOdist( pose.residue( residue ).xyz( "O" ).distance( pose.residue( residue + 1 ).xyz( "N" ) ) );
	return ( CNdist <= 1.4 ) and ( NOdist > 2.0 );
}

core::kinematics::FoldTree
NubInitioLoopClosureMover::make_final_tree( core::pose::Pose const & pose ) const {
	using namespace core::kinematics;
	FoldTree final_tree = pose.fold_tree();
	TR.Trace << "Initial tree: " << final_tree << std::endl;
	core::Size i = 1;
	while ( i <= final_tree.num_jump() ) {
		Edge jump = final_tree.jump_edge( i );
		TR.Trace << jump << std::endl;
		if ( pose.residue( jump.start() ).chain() == pose.residue( jump.stop()).chain() ) {
			final_tree.delete_jump_and_intervening_cutpoint( i );
		} else {
			++i;
		}
	}
	final_tree.renumber_jumps();
	TR.Trace << "Final tree: " << final_tree << std::endl;
	return final_tree;
}

bool
NubInitioLoopClosureMover::has_cutpoints( core::pose::Pose const & pose ) const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP chainbreaks( new selectors::CutpointResidueSelector( false ) );
	return has_any_true_selection( chainbreaks->apply(pose) );
}

void
NubInitioLoopClosureMover::show_cutpoints( core::pose::Pose const & pose ) const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP chainbreaks( new selectors::CutpointResidueSelector( false ) );
	utility::vector1< core::Size > residues = selection_positions( chainbreaks->apply( pose ) );
	TR << TR.Bold << "Target Residues:" << std::endl;
	TR << TR.Bold << pose.sequence() << std::endl;
	TR << TR.Bold << represent_residue_selector( chainbreaks->apply( pose ) ) << std::endl;
	for ( auto res : residues ) {
		TR << TR.Bold << "CHAINBREAK RESIDUE: " << res << std::endl;
	}
	TR << TR.Bold << pose.fold_tree() << TR.Reset << std::endl;
	TR << TR.Bold << pose.annotated_sequence() << TR.Reset << std::endl;
}

void
NubInitioLoopClosureMover::make_cutpoints_coherent_to_foldtree( core::pose::Pose & pose ) const {
	using namespace core::kinematics;
	FoldTree tree = pose.fold_tree();
	core::Size i = 1;
	while ( i <= tree.num_jump() ) {
		Edge jump = tree.jump_edge( i );
		if ( pose.residue( jump.start() ).chain() == pose.residue( jump.stop()).chain() ) {
			core::Size cutpoint = tree.cutpoint( i );
			if ( !pose.residue(cutpoint).has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT ) &&
					!pose.residue(cutpoint).has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) ) {
				core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, cutpoint );
				if ( cutpoint < pose.size() ) {
					if ( !pose.residue(cutpoint + 1).has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT ) &&
							!pose.residue(cutpoint + 1).has_variant_type( core::chemical::LOWER_TERMINUS_VARIANT ) ) {
						core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, cutpoint + 1 );
					}
				}
			}
		}
		++i;
	}
}

void
NubInitioLoopClosureMover::make_foldtree_coherent_to_cutpoints( core::pose::Pose & pose ) const {
	apply_closure_trust( pose, true );
}

void
NubInitioLoopClosureMover::repack( core::pose::Pose & pose )
{

	using namespace core::pack::task::operation;
	using namespace core::pack::task;
	using namespace core::select::residue_selector;

	//Selectors
	ResidueSelectorOP active_packable_selector = active_packable_residues();
	ResidueSelectorOP not_to_move_selector( new NotResidueSelector( active_packable_selector) );
	ResidueSelectorOP true_selector( new TrueResidueSelector() );

	// Make task Factory
	TaskFactoryOP postTaskFact( new TaskFactory );
	TaskOperationOP fixMotifOperation( new OperateOnResidueSubset( ResLvlTaskOperationCOP( new PreventRepackingRLT() ), not_to_move_selector ) );
	postTaskFact->push_back( fixMotifOperation );
	if ( not design_ ) {
		TR.Trace << "No sequence design" << std::endl;
		TaskOperationOP templateOperation( new OperateOnResidueSubset( ResLvlTaskOperationCOP( new RestrictToRepackingRLT() ), active_packable_selector ) );
		postTaskFact->push_back( templateOperation );
	} else {
		TR.Trace << "Sequence design allowed" << std::endl;
		DisallowIfNonnativeRLT mutableres;
		mutableres.disallow_aas("CYS");
		TaskOperationOP templateOperation( new OperateOnResidueSubset( mutableres.clone(), active_packable_selector ) );
		postTaskFact->push_back( templateOperation );
	}

	// Repack
	minimization_packing::PackRotamersMover repacker( fullatom_scorefxn_, postTaskFact->create_task_and_apply_taskoperations( pose ), 3 );
	repacker.apply( pose );

	if ( label_ ) {
		LabelPoseFromResidueSelectorMover labeler( active_packable_selector, "REPACKED" );
		labeler.apply( pose );
	}
}

// -- SELECTORS -- //

core::select::residue_selector::ResidueSelectorOP
NubInitioLoopClosureMover::movable_residues() const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP template_selector( new ResiduePDBInfoHasLabelSelector( "TEMPLATE" ) );
	ResidueSelectorOP flexible_selector( new ResiduePDBInfoHasLabelSelector( "FLEXIBLE" ) );
	ResidueSelectorOP movable_selector( new OrResidueSelector( template_selector, flexible_selector ) );
	return movable_selector;
}

core::select::residue_selector::ResidueSelectorOP
NubInitioLoopClosureMover::packable_residues() const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP coldspot_selector( new ResiduePDBInfoHasLabelSelector( "COLDSPOT" ) );
	ResidueSelectorOP template_selector( new ResiduePDBInfoHasLabelSelector( "TEMPLATE" ) );
	ResidueSelectorOP packable_residues( new OrResidueSelector( coldspot_selector, template_selector ) );
	return packable_residues;
}

bool
NubInitioLoopClosureMover::has_movable_residues( core::pose::Pose const & pose ) const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP movable_selector = movable_residues();
	return has_any_true_selection(movable_selector->apply( pose ) ) ;
}

core::select::residue_selector::ResidueSelectorOP
NubInitioLoopClosureMover::active_movable_residues(  core::pose::Pose const & pose  ) {
	using namespace core::select::residue_selector;
	ResidueSelectorOP chainbreaks( new selectors::CutpointResidueSelector( false ) );
	ResidueSelectorOP movable_selector = movable_residues();
	ResidueSelectorOP break_friends_selector( new PrimarySequenceNeighborhoodSelector( break_side_ + break_side_ramp_value_,
		break_side_ + 1 + break_side_ramp_value_,
		chainbreaks) );
	ResidueSelectorOP active_movable_selector( new AndResidueSelector( break_friends_selector, movable_selector ) );
	altered_residues_.set_residue_subset( active_movable_selector->apply( pose ) );
	return active_movable_selector;
}

core::select::residue_selector::ResidueSelectorOP
NubInitioLoopClosureMover::active_packable_residues() const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP torepack_selector( new NeighborhoodResidueSelector( altered_residues_.clone(), 10.01, true ) );
	ResidueSelectorOP active_packable_selector( new AndResidueSelector( torepack_selector, packable_residues() ) );
	return active_packable_selector;
}

// -- SET/GET/DEFAULT: CENTROID SCORE -- //

core::scoring::ScoreFunctionOP
NubInitioLoopClosureMover::centroid_scorefxn() const {
	return centroid_scorefxn_;
}

void
NubInitioLoopClosureMover::centroid_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ) {
	centroid_scorefxn_ = scorefxn;
}

core::scoring::ScoreFunctionOP
NubInitioLoopClosureMover::default_centroid_scorefxn() {
	core::scoring::ScoreFunctionOP tmp = core::scoring::ScoreFunctionFactory::create_score_function( "cen_std", "score4L" );
	return set_centroid_default_chainclosure_weight(tmp);
}

core::scoring::ScoreFunctionOP
NubInitioLoopClosureMover::set_centroid_default_chainclosure_weight(core::scoring::ScoreFunctionOP const & scorefxn ) {
	core::scoring::ScoreFunctionOP tmp = scorefxn->clone();
	tmp->set_weight_if_zero( core::scoring::chainbreak, default_centroid_chainclosure_weight() );
	tmp->set_weight_if_zero( core::scoring::linear_chainbreak, default_centroid_chainclosure_weight() );
	tmp->set_weight_if_zero( core::scoring::overlap_chainbreak, default_centroid_chainclosure_weight() );
	return tmp;
}

core::Real
NubInitioLoopClosureMover::default_centroid_chainclosure_weight() {
	return 1.5;
}

// -- SET/GET/DEFAULT: FULLATOM SCORE -- //

core::scoring::ScoreFunctionOP
NubInitioLoopClosureMover::fullatom_scorefxn() const {
	return fullatom_scorefxn_;
}

void
NubInitioLoopClosureMover::fullatom_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ) {
	fullatom_scorefxn_ = scorefxn;
}

core::scoring::ScoreFunctionOP
NubInitioLoopClosureMover::default_fullatom_scorefxn() {
	return core::scoring::get_score_function();;
}

// -- SET/GET/DEFAULT: BREAK PROPERTIES -- //

core::Size
NubInitioLoopClosureMover::break_side() const {
	return break_side_;
}

void
NubInitioLoopClosureMover::break_side( core::Real value ) {
	break_side_ = value;
}

core::Size
NubInitioLoopClosureMover::default_break_side() {
	return 3;
}

bool
NubInitioLoopClosureMover::break_side_ramp() const {
	return break_side_ramp_;
}

void
NubInitioLoopClosureMover::break_side_ramp( bool pick ) {
	break_side_ramp_ = pick;
}

bool
NubInitioLoopClosureMover::default_break_side_ramp() {
	return false;
}

// -- SET/GET/DEFAULT: TRIALS -- //

void
NubInitioLoopClosureMover::max_trials( core::Size choice ) {
	max_trials_ = choice;
}

core::Size
NubInitioLoopClosureMover::max_trials() const {
	return max_trials_;
}

core::Size
NubInitioLoopClosureMover::default_max_trials() {
	return 10;
}

// -- SET/GET/DEFAULT: TRUST -- //

bool
NubInitioLoopClosureMover::trust() const {
	return trust_;
}
void
NubInitioLoopClosureMover::trust( bool pick ) {
	trust_ = pick;
}

bool
NubInitioLoopClosureMover::default_trust() {
	return false;
}

// -- SET/GET/DEFAULT: LABEL -- //

bool
NubInitioLoopClosureMover::label() const {
	return label_;
}
void
NubInitioLoopClosureMover::label( bool pick ) {
	label_ = pick;
}

bool
NubInitioLoopClosureMover::default_label() {
	return true;
}

// -- SET/GET/DEFAULT: CENTROID -- //

bool
NubInitioLoopClosureMover::centroid() const {
	return centroid_;
}
void
NubInitioLoopClosureMover::centroid( bool pick ) {
	centroid_ = pick;
}

bool
NubInitioLoopClosureMover::default_centroid() {
	return false;
}

// -- SET/GET/DEFAULT: CENTROID -- //

bool
NubInitioLoopClosureMover::design() const {
	return design_;
}
void
NubInitioLoopClosureMover::design( bool pick ) {
	design_ = pick;
}

bool
NubInitioLoopClosureMover::default_design() {
	return false;
}

// -- SET/GET: FRAGMENTS -- //

std::string
NubInitioLoopClosureMover::fragments_id() const {
	return fragments_id_;
}

void
NubInitioLoopClosureMover::fragments_id( std::string const & name ) {
	fragments_id_ = name;
}

core::fragment::FragSetOP
NubInitioLoopClosureMover::fragments() const {
	return fragments_;
}
void
NubInitioLoopClosureMover::fragments( core::fragment::FragSetOP frags) {
	fragments_ = frags;
}

// -- ROSETTASCRIPTS -- //
protocols::moves::MoverOP
NubInitioLoopClosureMover::clone() const
{
	return protocols::moves::MoverOP( new NubInitioLoopClosureMover( *this ) );
}

protocols::moves::MoverOP
NubInitioLoopClosureMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new NubInitioLoopClosureMover );
}

std::string NubInitioLoopClosureMover::get_name() const {
	return mover_name();
}

std::string NubInitioLoopClosureMover::mover_name() {
	return "NubInitioLoopClosureMover";
}

// -- MOVERCREATOR -- //
std::string NubInitioLoopClosureMoverCreator::keyname() const {
	return NubInitioLoopClosureMover::mover_name();
}

protocols::moves::MoverOP
NubInitioLoopClosureMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new NubInitioLoopClosureMover );
}

void NubInitioLoopClosureMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	NubInitioLoopClosureMover::provide_xml_schema( xsd );
}

}
} //protocols
} //fold_from_loops
