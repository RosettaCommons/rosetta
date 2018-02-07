// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   NubInitioMover.hh
/// @brief  Ab Initio with a pre-folded segment (nub - pun intended)
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

// Unit headers
#include <protocols/moves/Mover.hh>
#include <protocols/fold_from_loops/NubInitioMover.hh>
#include <protocols/fold_from_loops/NubInitioMoverCreator.hh>
#include <protocols/fold_from_loops/filters/RmsdFromResidueSelectorFilter.hh>
#include <protocols/fold_from_loops/selectors/ConstraintResidueSelector.hh>
#include <protocols/fold_from_loops/movers/AlignByResidueSelectorMover.hh>
#include <protocols/fold_from_loops/movers/SplitAndMixPoseMover.hh>
#include <protocols/fold_from_loops/movers/ReleaseConstraintFromResidueMover.hh>
#include <protocols/fold_from_loops/movers/LabelPoseFromResidueSelectorMover.hh>
#include <protocols/fold_from_loops/movers/DisplayPoseLabelsMover.hh>
#include <protocols/fold_from_loops/movers/NubInitioLoopClosureMover.hh>
#include <protocols/fold_from_loops/selectors/CutpointResidueSelector.hh>
#include <protocols/fold_from_loops/utils/utils.hh>
#include <protocols/fold_from_loops/utils/Nub.hh>

// Protocol headers
#include <protocols/moves/mover_schemas.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/DsspMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/DumpPdb.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/abinitio/FoldConstraints.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/loops/loop_closure/ccd/SlidingWindowLoopClosure.hh>
#include <protocols/loops/Exceptions.hh>
#include <protocols/jumping/util.hh>
#include <protocols/checkpoint/CheckPointer.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/util/disulfide_util.hh>
#include <core/conformation/util.hh>
#include <core/fragment/FragSet.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/AndResidueSelector.hh>
#include <core/select/residue_selector/OrResidueSelector.hh>
#include <core/select/residue_selector/PrimarySequenceNeighborhoodSelector.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/ResLvlTaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <string>

namespace protocols {
namespace fold_from_loops {

static basic::Tracer TR( "protocols.fold_from_loops.NubInitioMover" );

NubInitioMover::NubInitioMover():
	protocols::moves::Mover( mover_name() ),
	prefix_( default_prefix() ),
	template_pose_( core::pose::Pose() ),
	template_selector_( default_template_selector() ),
	fullatom_scorefxn_( default_fullatom_scorefxn() ),
	nub_( new utils::Nub() ),
	use_cst_( default_use_cst() ),
	clear_motif_cst_( default_clear_motif_cst() ),
	rmsd_threshold_( default_rmsd_threshold() ),
	rmsd_include_motif_( default_rmsd_include_motif() ),
	rmsd_include_unconstrained_( default_rmsd_include_unconstrained() ),
	repack_disulfides_( default_repack_disulfides() ),
	disulfides_bb_( default_disulfides_bb() ),
	disulfides_side_( default_disulfides_side() ),
	fragments_id_( std::string() ),
	trials_( 1 ),
	max_trials_( default_max_trials() ),
	mc_binder_weight_( default_mc_binder_weight() ),
	mc_angle_weight_( default_mc_angle_weight() ),
	mc_dihedral_weight_( default_mc_dihedral_weight() ),
	mc_correction_weights_( default_mc_correction_weights() ),
	use_correction_weights_( default_use_correction_weights() ),
	dump_centroid_( default_dump_centroid() ),
	drop_unfolded_pose_( default_drop_unfolded_pose() ),
	design_( default_design() ),
	residue_type_( default_residue_type() )
{}

NubInitioMover::~NubInitioMover(){}

void
NubInitioMover::apply( core::pose::Pose & pose )
{

	make_unfolded_pose( pose );         // 1. CREATE UNFOLDED POSE
	core::Size abinitio = refold_pose( pose ); // 2. REFOLD & EVALUATE

	if ( abinitio == 1 ) {           // 3. POSTPROCESS
		// We stop here, providing only the unfolded pose.
		set_last_move_status( moves::MS_SUCCESS );
	} else if ( abinitio == 2 ) {
		// We have been unable to obtain a folded pose under the RMSD threshold
		set_last_move_status( moves::FAIL_DO_NOT_RETRY );
	} else {
		// abinitio == 0: success!
		post_process( pose );
		set_last_move_status( moves::MS_SUCCESS );
	}
	(*fullatom_scorefxn_)( pose );
}

/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
void
NubInitioMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & reference_pose )
{

	prefix( tag->getOption< std::string >( "name", default_prefix() ) );
	template_selector( core::select::residue_selector::parse_residue_selector( tag, data, "template_motif_selector" ) );
	if ( tag->hasOption( "fullatom_scorefxn" ) ) {
		fullatom_scorefxn( protocols::rosetta_scripts::parse_score_function( tag, "fullatom_scorefxn", data ) );
	}
	max_trials( tag->getOption< core::Size >( "max_trials", default_max_trials() ) );

	TR.Trace << TR.Green << "Loading constraint configuration..." << TR.Reset << std::endl;
	use_cst( tag->getOption< bool >( "use_cst", default_use_cst() ) );
	clear_motif_cst( tag->getOption< bool >( "clear_motif_cst", default_clear_motif_cst() ) );

	TR.Trace << TR.Green << "Loading RMSD configuration..." << TR.Reset << std::endl;
	rmsd_include_motif( tag->getOption< bool >( "rmsd_include_motif", default_rmsd_include_motif() ) );
	rmsd_include_unconstrained( tag->getOption< bool >( "rmsd_include_unconstrained", default_rmsd_include_unconstrained() ) );
	rmsd_threshold( tag->getOption< core::Real >( "rmsd_threshold", default_rmsd_threshold() ) );

	binder_weight( tag->getOption< core::Real >( "binder_weight", default_mc_binder_weight() ) );
	angle_weight( tag->getOption< core::Real >( "angle_weight", default_mc_angle_weight() ) );
	dihedral_weight( tag->getOption< core::Real >( "dihedral_weight", default_mc_dihedral_weight() ) );
	use_correction_weights( tag->getOption< bool > ( "correction_weights", default_use_correction_weights() ) );

	repack_disulfides( tag->getOption< bool > ( "repack_disulfides", default_repack_disulfides() ) );
	disulfides_bb( tag->getOption< bool > ( "disulfides_bb", default_disulfides_bb() ) );
	disulfides_side( tag->getOption< core::Size >( "disulfides_side", default_disulfides_side() ) );

	design( tag->getOption< bool > ( "design", default_design() ) );
	if ( tag->hasOption( "residue_type" ) ) {
		design( true );
		residue_type( tag->getOption< std::string >( "residue_type", default_residue_type() ) );
	}

	TR.Trace << TR.Green << "Loading Fragments configuration..." << TR.Reset << std::endl;
	fragments_id( tag->getOption< std::string >( "fragments_id", std::string() ) );
	if ( data.has( "FragSet", fragments_id_ + "small" ) && data.has( "FragSet", fragments_id_ + "large" ) ) {
		nub_->small_fragments( data.get_ptr<core::fragment::FragSet>( "FragSet", fragments_id_ + "small") );
		nub_->large_fragments( data.get_ptr<core::fragment::FragSet>( "FragSet", fragments_id_ + "large") );

		TR.Trace << TR.Green << "Fragments loaded through the basic::datacache::DataMap" << TR.Reset << std::endl;
		TR.Trace << TR.Green << "Small: " << data.get_ptr<core::fragment::FragSet>( "FragSet", fragments_id_ + "small") << TR.Reset << std::endl;
		TR.Trace << TR.Green << "Large: " << data.get_ptr<core::fragment::FragSet>( "FragSet", fragments_id_ + "large") << TR.Reset << std::endl;
	}
	dump_centroid( tag->getOption< bool >( "dump_centroid", default_dump_centroid() ) );
	drop_unfolded_pose( tag->getOption< bool >( "drop_unfolded_pose", default_drop_unfolded_pose() ) );

	nub_->parse_tag( tag, data, reference_pose );

}

void
NubInitioMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ optional_name_attribute()
		+ XMLSchemaAttribute::attribute_w_default( "use_cst", xsct_rosetta_bool,
		"Use constraints to guide the folding (RECOMENDED). Basically this is the difference between running ClassicAbinitio or FoldConstraints. "
		"Be aware that the constraints MUST have been set into the Pose beforehand with a ConstraintGenerator.", std::to_string( default_use_cst() ) )
		+ XMLSchemaAttribute::attribute_w_default( "max_trials", xs_integer,
		"Defines how many times we should try to fold the protein under a RMSD threshold with the template before giving up.", std::to_string( default_max_trials() ))
		+ XMLSchemaAttribute::attribute_w_default( "clear_motif_cst", xsct_rosetta_bool,
		"This will clear any constraint that includes residues of the template's motif; i.e. the part of the template that will be "
		"substituted by the Nub. This option is relevant when the 3D configuration of the motif and the Nub differ greatly. "
		"Be aware that this change is applied to the pose, and, thus, these constraints will not be available in further Movers. "
		"(WARNING!) When the size of motif and Nub do not match, this should be kept as TRUE.", std::to_string( default_clear_motif_cst() ) )
		+ XMLSchemaAttribute::attribute_w_default( "rmsd_threshold", xsct_real,
		"Maximum allow difference between the obtained decoy and the original template.", std::to_string( default_rmsd_threshold() ) )
		+ XMLSchemaAttribute::attribute_w_default( "rmsd_include_motif", xsct_rosetta_bool,
		"Use all the residues of the template, including those of the motif, to calculate RMSD."
		"(WARNING!) When the size of motif and Nub do not match, this should be kept as FALSE.", std::to_string( default_rmsd_include_motif() ) )
		+ XMLSchemaAttribute::attribute_w_default( "rmsd_include_unconstrained", xsct_rosetta_bool,
		"Use all the residues of the template, including those that had no assigned constraint (this excludes SequenceComposition constraints). "
		"This is interesting if a segment of the template was left unsconstraint on purpose so that it would move for example, to avoid clashed with a binder. "
		"It would allow to set tighter RMSD filtering while leaving this particular region alone,", std::to_string( default_rmsd_include_unconstrained() ) )
		+ XMLSchemaAttribute::attribute_w_default( "binder_weight", xsct_real,
		"Weight for interchain scores to add to the AbInitio MonteCarlo's evaluator default score (0 - no weight).", std::to_string( default_mc_binder_weight() ) )
		+ XMLSchemaAttribute::attribute_w_default( "angle_weight", xsct_real,
		"Weight for angle constraints to add to the AbInitio MonteCarlo's evaluator default score (0 - no weight).", std::to_string( default_mc_angle_weight() ) )
		+ XMLSchemaAttribute::attribute_w_default( "dihedral_weight", xsct_real,
		"Weight for dihedral constraints to add to the AbInitio MonteCarlo's evaluator default score (0 - no weight).", std::to_string( default_mc_dihedral_weight() ) )
		+ XMLSchemaAttribute::attribute_w_default( "correction_weights", xsct_rosetta_bool,
		"When true (default), it will autodetect the presence of helix/sheet in the template and add scoring terms appropiately.", std::to_string( default_use_correction_weights() ) )
		+ XMLSchemaAttribute::attribute_w_default( "repack_disulfides", xsct_rosetta_bool,
		"Force disulfide awarenes to guide ab initio.", std::to_string( default_repack_disulfides() ) )
		+ XMLSchemaAttribute::attribute_w_default( "disulfides_bb", xsct_rosetta_bool,
		"Allow backbone movements to try to fix disulfides.", std::to_string( default_disulfides_bb() ) )
		+ XMLSchemaAttribute::attribute_w_default( "disulfides_side", xs_integer,
		"Defines number of sequence neighbors around CYS residues allowed to pack/minimize in order to achieve the disulfide bridge.",
		std::to_string( default_disulfides_side() ))
		+ XMLSchemaAttribute::required_attribute( "fragments_id", xs_string,
		"Fragments are necessary to perform the ab initio run. They need to be included/created with StructFragmentMover. "
		"The value set in the 'prefix' attribute of that Mover needs to be provided here again.")
		+ XMLSchemaAttribute::attribute_w_default( "dump_centroid", xsct_rosetta_bool,
		"Output centroid level structures in a outputname_CENTROID silent file.", std::to_string( default_dump_centroid() ) )
		+ XMLSchemaAttribute::attribute_w_default( "drop_unfolded_pose", xsct_rosetta_bool,
		"For testing and imaging. Stops the process just after creating the unfolded pose (still full atom)", std::to_string( default_drop_unfolded_pose() ) )
		+ XMLSchemaAttribute::attribute_w_default( "design", xsct_rosetta_bool,
		"If true, run design on the template residues", std::to_string( default_design() ) )
		+ XMLSchemaAttribute::attribute_w_default( "residue_type", xs_string, "Change all template residues to the specified residue name (name1)", default_residue_type() );

	protocols::rosetta_scripts::attributes_for_parse_score_function_w_description( attlist, "fullatom_scorefxn",
		"Full atom score for disulfide optimization and repacking, if needed. If not specified calls Rosetta's default full atom score." );
	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "template_motif_selector",
		"Selector specifying the residues from the template that will be substituted for the Nub; i.e. everything from the template that"
		"we WON'T KEEP. This region is dubed as the template's motif.");

	XMLSchemaSimpleSubelementList subelements;
	utils::Nub().provide_xml_definition( xsd, subelements ); // Add the Nub Definition to the xsd

	// Adding everything to the xsd
	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & protocols::moves::complex_type_name_for_mover )
		.element_name( mover_name() )
		.description( "The Nub Initio mover generates a low definition (centroid level) ab initio folding (with or withou constraints) "
		"around a foreing static piece of structure called Nub. Effectively, it grafts the Nub inside a template protein "
		"and then folds that protein around to better fit the inserted Nub. The region of the template that will be substituted "
		"(i.e. not kept) is the 'motif'. The pose that comes at apply time is considered the template. "
		"If the template pose is discontinued, this might produce unexpected and possibly fatal results. "
		"Finally a full atom pose is provided keeping the original residue types of the template. "
		"This pose should not be consider OK to terminate the run at that point; at least the template residues need to be redesigned.")
		.add_attributes( attlist )
		.set_subelements_single_appearance_required( subelements )
		.write_complex_type_to_schema( xsd );
}

////////////////////////////////////////////////
////////////////// MAIN STEPS //////////////////
////////////////////////////////////////////////

/// @details STEP1 of NubInitioMover
/// This checks that all inputs are as expected and that the starting pose (unfolded) is properly generated.
/// This function is only runned once, regardless of how many structures are requested; thus, the NubInitioMover
/// in a given run can only be runned for the same funtional motif - template pair.
void
NubInitioMover::make_unfolded_pose( core::pose::Pose & pose ) {

	if ( nub_->unfolded_pose()->empty() ) {
		protocols::moves::DsspMover dssp;
		dssp.apply( pose );

		TR << "Evaluation of input data's coherence" << std::endl;
		sanity_check( pose );

		TR << "Detecting disulfides" << std::endl;
		utility::vector1< std::pair< core::Size, core::Size > > disulfides;
		core::conformation::disulfide_bonds( pose.conformation(), disulfides );

		TR << "Processing the Template Segments to keep" << std::endl;
		utility::vector1< core::pose::PoseOP > template_pieces = get_template_pieces( pose );
		if ( TR.Trace.visible() ) {
			core::Size count_piece = 1;
			for ( auto piece: template_pieces ) {
				TR.Trace << "Piece " << count_piece << ": " << piece->sequence() << std::endl;
				++count_piece;
			}
		}

		TR << "Creating the Unfolded Pose." << std::endl;
		nub_->apply( pose, template_selector_->apply( pose ), template_pieces, disulfides );

		TR << "Managing Template Constraints" << std::endl;
		manage_constraints( pose );    // Make sure Constraints are available if requested.
	}
}

/// @details STEP2 of NubInitioMover
/// This is the main job of the mover. It will execute the ab initio
/// and it will evaluate the result. 0: success; 1: undolded pose (debug & showcase); 2: fail
core::Size
NubInitioMover::refold_pose( core::pose::Pose & pose )
{
	trials_ = 1;
	movers::DisplayPoseLabelsMover display;
	simple_moves::SwitchResidueTypeSetMoverOP switcher( new simple_moves::SwitchResidueTypeSetMover );
	while ( trials_ < max_trials_ ) {
		TR << "Transfer Pose and Label Residues" << std::endl;
		nub_->transfer_unfolded_conformation( pose );

		TR << "Unfolded Pose Summary" << std::endl;
		display.movemap_factory( nub_->movemapfactory() );
		display.apply( pose );

		// This cuts the process in half. For debuging and presentation-creating purposes only!!
		if ( drop_unfolded_pose_ ) {
			TR << "Abruptly stopping the process by user's request!" << std::endl;
			return 1;
		}

		TR << "Running Ab Initio protocol" << std::endl;
		switcher->type_set_tag( core::chemical::CENTROID );
		switcher->apply( pose );
		abinitio_score_ = apply_abinitio( pose );

		TR << "Evaluating RMSD-threshold success" << std::endl;
		core::Real rmsd = template_rmsd( pose );

		if ( rmsd > rmsd_threshold_ ) {
			TR << TR.Red << "RMSD of " << rmsd << ": over threshold" << TR.Reset << std::endl;
			TR << max_trials_ - trials_ << " trials left" << std::endl;
			pose.detached_copy( template_pose_ );  // reset pose
			++trials_;
		} else {
			TR << TR.Green << "RMSD of " << rmsd << ": under threshold" << TR.Reset << std::endl;
			core::pose::setPoseExtraScore( pose, nub_->design_chain() + "_ni_trials", trials_ );
			display.apply( pose );
			return 0;
		}
	}
	return 2;
}

/// @details STEP3 of the NubInitioMover
/// Close loops, score, recover full atoms, pack and/or design
void
NubInitioMover::post_process( core::pose::Pose & pose ) {

	using namespace core::chemical;
	protocols::moves::DsspMover dssp;
	movers::DisplayPoseLabelsMover display;
	simple_moves::SwitchResidueTypeSetMoverOP switcher( new simple_moves::SwitchResidueTypeSetMover );

	// Save this for the end
	core::kinematics::FoldTree working_tree = pose.fold_tree();

	// Close Loops, if any
	if ( nub_->is_multisegment() ) {
		TR.Debug << "LOOP CLOSURE: START" << std::endl;
		movers::NubInitioLoopClosureMover lclose;
		lclose.break_side_ramp( true );
		lclose.label( false );
		lclose.centroid( true );
		lclose.fragments( nub_->small_fragments() );
		lclose.apply( pose );
		TR.Debug << "LOOP CLOSURE: END" << std::endl;
		display.apply( pose ); // Verbose
	}
	pose.update_residue_neighbors();

	// Score and count contacts of the inserted motif with the scaffold.
	(*abinitio_score_)( pose );
	core::pose::setPoseExtraScore( pose, nub_->design_chain() + "_ni_mtcontacts", count_contacts( pose ) );

	// Dump centroid if requested.
	pose.fold_tree( working_tree );
	dump_centroid( pose, abinitio_score_ );

	TR << "Refitting side chains" << std::endl;
	switcher->type_set_tag( core::chemical::FA_STANDARD );
	switcher->apply( pose );
	nub_->refit_motif_sidechains( pose );
	repack_minimize_disulfides( pose );
	display.apply( pose ); // Verbose
	pose.update_residue_neighbors();
	repack( pose );
	pose.update_residue_neighbors();

	// Refit the open FoldTree
	TR << "Refit Working FoldTree" << std::endl;
	for ( core::Size i = 1; i < pose.size(); ++i ) {  // First Clean all variants
		core::pose::remove_variant_type_from_pose_residue( pose,  CUTPOINT_UPPER, i );
		core::pose::remove_variant_type_from_pose_residue( pose,  CUTPOINT_LOWER, i );
	}
	for ( core::Size i = 1; i < pose.size(); ++i ) { // Less than because cutpoints are between i and i+1
		if ( working_tree.is_cutpoint( i ) ) {
			if ( !pose.residue(i).has_variant_type( UPPER_TERMINUS_VARIANT ) ) {
				core::pose::correctly_add_cutpoint_variants( pose, i, false );
			}
		}
	}
	std::ostringstream ft;
	ft << pose.fold_tree();
	core::pose::add_comment( pose, "WORKING_FOLDTREE", ft.str() );
	core::pose::add_comment( pose, "POST_NUBINITIO_SEQ", pose.sequence() );

	TR << "After AbInitio Pose Summary" << std::endl;
	dssp.apply( pose );
	display.write( true );
	display.apply( pose );

}

////////////////////////////////////////////////
/////////////// SUPPORT FUNCTIONS //////////////
////////////////////////////////////////////////

void
NubInitioMover::sanity_check( core::pose::Pose const & pose ) {

	using namespace core::select::residue_selector;

	TR << TR.Blue <<  "Is the template a single chain?" << std::endl;
	runtime_assert_msg(pose.num_chains() == 1 , "Template can only have a single chain." );
	TR << TR.Green << " * Yes!" << std::endl;

	TR << TR.Blue << "Are there segments of the template that we will keep?" << std::endl;
	ResidueSubset subset = used_template_selector()->apply( pose );
	runtime_assert_msg( count_selected( subset ) > 0, "No residues were selected to be kept from the template." );
	TR << TR.Green << " * Yes!" << std::endl;

	TR << TR.Blue << "Are there insertion regions selected for the template?" << std::endl;
	subset = template_selector()->apply( pose );
	runtime_assert_msg( count_selected( subset ) > 0, "No residues were selected in the template as insertion segments." );
	TR << TR.Green << " * Yes!" << std::endl;

	TR << TR.Blue << "Seems template a Protein Sequence" << std::endl;
	TR << TR.Green << pose.sequence() << std::endl;
	TR << TR.Green << pose.secstruct() << std::endl;

	TR << TR.Blue << "Does the template have alpha content?" << std::endl;
	has_alphas_ = (pose.secstruct().find("H") != std::string::npos);
	if ( has_alphas_ ) {
		TR << TR.Green << " * Yes!" << std::endl;
	} else {
		TR << TR.Green << " * No!" << std::endl;
	}
	TR << TR.Blue << "Does the template have beta content?" << std::endl;
	has_betas_  = (pose.secstruct().find("E") != std::string::npos);
	if ( has_betas_ ) {
		TR << TR.Green << " * Yes!" << std::endl;
	} else {
		TR << TR.Green << " * No!" << std::endl;
	}
	// Copy Pose. To obtain further info once it has changed.
	template_pose_.detached_copy( pose );

	if ( use_cst_ ) {
		TR << TR.Blue << "Does the template have constraints?" << std::endl;
		runtime_assert_msg( pose.constraint_set()->has_residue_pair_constraints(),
			"If constraints are requested for the Pose, at least atom pair constraints must be provided." );
		TR << TR.Green << " * Yes!" << std::endl;
		template_pose_.transfer_constraint_set( pose );
	}
}

/// @details This provides the pieces of the template that will be kept and attached to the functional
/// motif as the unfolded segments. As a rule, the number of these pieces has to be twice the number of
/// functional motifs; thus, if N-terminal or C-terminal motif insertions exist, empty Poses can be provided
/// here in order to make this work properly.
utility::vector1< core::pose::PoseOP >
NubInitioMover::get_template_pieces( core::pose::Pose const & pose ) const
{
	movers::SplitAndMixPoseMover splitter;
	splitter.set_ranges( make_template_ranges( used_template_ranges( true ) ) );
	return splitter.apply_without_merge( pose );
}

/// @details This calculates if the original ranges in the template provided by the user are enough or not.
/// If needs be, it will split ranges between to insertion motifs to allow for cutpoinst or it will add empty
/// ranges to keep the 2 to 1 proportion between the insertion motifs and the unfolded template segments.
core::select::residue_selector::ResidueRangesOP
NubInitioMover::make_template_ranges( core::select::residue_selector::ResidueRanges original ) const {
	using namespace core::select::residue_selector;
	ResidueRange empty_range( 0, 0 );
	ResidueRangesOP new_ranges( new ResidueRanges );
	new_ranges->push_back( original[1] );

	if ( original.size() > 2 ) { // this is only needed for multisegment motifs
		for ( core::Size i = 2; i <= (original.size()-1); ++i ) {
			ResidueRange range  = original[i];
			std::string wstruct = template_pose_.secstruct().substr( range.start() - 1, range.stop() - range.start() + 1 );
			// If the difference is only 1 residue, we put one on one side and an empty pose in the other.
			if ( wstruct.size() == 1 ) {
				new_ranges->push_back( range );
				new_ranges->push_back( empty_range );
			} else {
				core::Size cutpoint = utils::find_cutpoint_from_secondary_structure( wstruct ) + range.start();
				ResidueRange range1( range.start(), cutpoint );
				new_ranges->push_back( range1 );
				ResidueRange range2( cutpoint + 1, range.stop() );
				new_ranges->push_back( range2 );
			}
		}
	}

	new_ranges->push_back( original[ original.size() ] );
	return new_ranges;
}

/// @details Here the constraint options are set. Also, the score terms that inform in the
/// context in which the constraints are used is added.
/// @remark This function can change the ConstraintSet of the Pose by: (1) removing them from
/// the nub (static) segments and (2) bt remaping them to the new legnth of the unfolded pose if needed.
void
NubInitioMover::manage_constraints( core::pose::Pose & pose )
{
	using namespace core::scoring::constraints;
	if ( use_cst_ ) {
		core::pose::add_score_line_string( pose, nub_->design_chain() + "_ni_cst", "true" );
		core::pose::setPoseExtraScore( pose, nub_->design_chain() + "_ni_cst", "true" );
		if ( clear_motif_cst_ ) {
			movers::ReleaseConstraintFromResidueMover release( template_selector_ );
			release.apply( pose );
			core::pose::add_score_line_string( pose, nub_->design_chain() + "_ni_motif_cst", "false" );
			core::pose::setPoseExtraScore( pose, nub_->design_chain() + "_ni_motif_cst", "false" );
		} else {
			core::pose::add_score_line_string( pose, nub_->design_chain() + "_ni_motif_cst", "true" );
			core::pose::setPoseExtraScore( pose, nub_->design_chain() + "_ni_motif_cst", "true" );
		}
		if ( not nub_->template_to_unfolded_mapping().is_identity() ) {
			ConstraintSetOP constraints( new ConstraintSet( *pose.constraint_set() ) );
			constraints->remap_residue_positions( nub_->template_to_unfolded_mapping() );
			pose.constraint_set( constraints );
		}
	} else {
		core::pose::add_score_line_string( pose, nub_->design_chain() + "_ni_cst", "false" );
		core::pose::setPoseExtraScore( pose, nub_->design_chain() + "_ni_cst", "false" );
	}
}

/// @details Sets up the extra weights that seem appropiate for each case and runs the ab initio simulation.
core::scoring::ScoreFunctionOP
NubInitioMover::apply_abinitio( core::pose::Pose & pose )
{
	using namespace protocols::abinitio; // ClassicAbinitio, FoldConstraints

	ClassicAbinitioOP abinitio( new ClassicAbinitio( nub_->small_fragments(), nub_->large_fragments(), nub_->movemap()->clone() ) );
	if ( use_cst_ ) {
		abinitio = ClassicAbinitioOP( new FoldConstraints( nub_->small_fragments(), nub_->large_fragments() , nub_->movemap()->clone() ) );
	}
	abinitio->init( pose );
	if ( use_correction_weights_ and has_alphas_ ) {
		abinitio->set_score_weight( core::scoring::hbond_sr_bb, mc_correction_weights_ );
	}
	if ( use_correction_weights_ and has_betas_ ) {
		abinitio->set_score_weight( core::scoring::hbond_lr_bb, mc_correction_weights_ );
		abinitio->set_score_weight( core::scoring::rsigma,      mc_correction_weights_ );
		abinitio->set_score_weight( core::scoring::sheet,       mc_correction_weights_ );
		abinitio->set_score_weight( core::scoring::ss_pair,     mc_correction_weights_ );
	}
	if ( use_correction_weights_ and has_alphas_ and has_betas_ ) {
		abinitio->set_score_weight( core::scoring::hs_pair, mc_correction_weights_ );
	}
	if ( use_correction_weights_ and nub_->has_binder() and mc_binder_weight_ > 0 ) {
		abinitio->set_score_weight( core::scoring::interchain_pair,    mc_binder_weight_ );
		abinitio->set_score_weight( core::scoring::interchain_env,     mc_binder_weight_ );
		abinitio->set_score_weight( core::scoring::interchain_contact, mc_binder_weight_ );
		abinitio->set_score_weight( core::scoring::interchain_vdw,     mc_binder_weight_ );
	}
	if ( nub_->disulfides().size() > 0 and repack_disulfides_ ) {
		abinitio->set_score_weight( core::scoring::dslfc_cen_dst, 5.0 );
		abinitio->set_score_weight( core::scoring::dslfc_cb_dst,  4.0 );
		abinitio->set_score_weight( core::scoring::dslfc_ang,     2.5 );
		abinitio->set_score_weight( core::scoring::dslfc_cb_dih,  0.5 );
		abinitio->set_score_weight( core::scoring::dslfc_bb_dih,  0.5 );
	}
	if ( mc_angle_weight_ > 0 ) {
		abinitio->set_score_weight( core::scoring::angle_constraint, mc_angle_weight_ );
	}
	if ( mc_dihedral_weight_ > 0 ) {
		abinitio->set_score_weight( core::scoring::dihedral_constraint, mc_dihedral_weight_ );
	}
	abinitio->apply( pose );
	//utils::report_unfolded( pose, nub_->movemap() );
	return abinitio->mc().score_function().clone();
}

/// @details Compare the new folded pose with the template.
/// The function compares according to whatever options the user has picked and sets
/// the appropiate score terms to log those decisions.
core::Real
NubInitioMover::template_rmsd( core::pose::Pose & pose )
{
	using namespace core::select::residue_selector;

	filters::RmsdFromResidueSelectorFilterOP rmsd( new filters::RmsdFromResidueSelectorFilter );
	rmsd->reference_pose( template_pose_ );

	ResidueSelectorOP cstsele( new selectors::ConstraintResidueSelector() ) ;
	std::string infovalue = "full";

	if ( !rmsd_include_motif() ) {
		if ( rmsd_include_unconstrained() ) {
			infovalue = "no_motif";
			rmsd->selectors( used_template_selector(), work_template_selector() );
		} else {
			infovalue = "no_motif_no_constrained";
			ResidueSelectorOP designcst( new AndResidueSelector( work_template_selector(), cstsele) );
			ResidueSelectorOP originalcst( new AndResidueSelector( used_template_selector(), cstsele) );
			rmsd->selectors( originalcst, designcst );
		}
	} else {
		if ( !rmsd_include_unconstrained() ) {
			infovalue = "no_constrained";
			rmsd->selectors( cstsele, cstsele );
		} else {
			ResidueSelectorOP all_template( new TrueResidueSelector() );
			rmsd->selectors( all_template, design_chain_selector());
		}
	}
	core::Real rmsd_num = rmsd->compute( pose );
	core::pose::add_score_line_string( pose, nub_->design_chain() + "_ni_rmsd_type",  infovalue );
	core::pose::setPoseExtraScore( pose, nub_->design_chain() + "_ni_rmsd_type",      infovalue );
	core::pose::setPoseExtraScore( pose, nub_->design_chain() + "_ni_rmsd",           rmsd_num );
	core::pose::setPoseExtraScore( pose, nub_->design_chain() + "_ni_rmsd_threshold", rmsd_threshold_ );

	TR << "RMSD with template is: " << rmsd_num << std::endl;
	return rmsd_num;
}

core::Size
NubInitioMover::count_contacts( core::pose::Pose & pose ) const
{
	using namespace core::select::residue_selector;

	TR << "Calculating number of contacts of the MOTIF with the SCAFFOLD" << std::endl;
	// We'll exclude binder contacts and the 2 residues on each side of the insertion.
	ResidueSelectorOP seqsel( new PrimarySequenceNeighborhoodSelector( 2, 2, work_motif_selector() ) );
	ResidueSelectorOP exclude_inv( new OrResidueSelector( seqsel, context_chain_selector() ) );
	ResidueSelectorOP exclude( new NotResidueSelector( exclude_inv ) );
	TR.Trace << "EXCLUDE " << represent_residue_selector( exclude->apply(pose) ) << std::endl;

	// Get neighbors, remove excluded and count. Distance is 10.1 to force the calculation instead of
	// using the 10 angstrom matrix, which does not seem to work with CENTROID.
	ResidueSelectorOP contacts( new NeighborhoodResidueSelector( work_motif_selector(), 10, true ) );
	TR.Trace << "CONTACT " << represent_residue_selector( contacts->apply(pose) ) << std::endl;
	ResidueSelectorOP query( new AndResidueSelector( exclude, contacts ) );

	// Label them (here, in the original pose)
	movers::LabelPoseFromResidueSelectorMover labeler( query, nub_->contact_label() );
	labeler.apply( pose );

	// Return count
	return count_selected( query->apply( pose ) );
}

void
NubInitioMover::dump_centroid( core::pose::Pose const & pose, core::scoring::ScoreFunctionOP scorefxn )
{
	using namespace core::io::silent;
	movers::DisplayPoseLabelsMover display;
	std::ostringstream ft;
	if ( dump_centroid_ and jd2::jd2_used() ) {
		TR << "Dump into silent file" << std::endl;
		SilentFileOptionsOP silent_options( new SilentFileOptions );
		silent_score_file_ = SilentFileDataOP( new SilentFileData( *silent_options ) );
		silent_score_file_->set_filename( jd2::current_output_filename() + "_CENTROID" );
		core::pose::Pose cpose = pose;
		(*scorefxn)( cpose );
		display.write( true );
		display.apply( cpose );
		ft << cpose.fold_tree();
		core::pose::add_comment( cpose, "WORKING_FOLDTREE", ft.str() );
		SilentStructOP ss = SilentStructFactory::get_instance()->get_silent_struct_out( *silent_options );
		std::string pose_tag( jd2::current_output_name() );
		ss->fill_struct( cpose, pose_tag + "_CENTROID" );
		silent_score_file_->write_silent_struct( *ss,  silent_score_file_->filename() );
	}
}

/// @details If there are disulfides in the scaffold or the motif, it tries to fix them as well as
/// possible.
void
NubInitioMover::repack_minimize_disulfides( core::pose::Pose & pose )
{
	if ( nub_->disulfides().size() > 0 and repack_disulfides_ ) {

		using namespace core::pack::task;
		using namespace core::select::residue_selector;
		TR.Debug << "REPACK DISULFIDES: START" << std::endl;
		TR << "Fixing disulfides... " << std::endl;
		for ( auto dis : nub_->disulfides() ) {
			TR << "CYD pair is: " << dis.first << " " << dis.second << std::endl;
		}

		// Setup Scores
		TR << "Setting up scores" << std::endl;
		core::scoring::ScoreFunctionOP packer_score = fullatom_scorefxn_->clone();
		packer_score->set_weight( core::scoring::dslf_fa13, 1.25 );
		packer_score->set_weight( core::scoring::fa_dun,    0.7  );
		core::scoring::ScoreFunctionOP minimizer_score = packer_score->clone();
		minimizer_score->set_weight( core::scoring::dslf_fa13, 4 );
		minimizer_score->set_weight( core::scoring::fa_dun,    0 );

		// Setup Selectors
		TR << "Setting up selectors" << std::endl;
		ResidueSelectorOP CYS_side_selector( new PrimarySequenceNeighborhoodSelector( disulfides_side_, disulfides_side_, cys_design_selector() ) );
		ResidueSelectorOP CYS_bb_movable( new AndResidueSelector(CYS_side_selector, bb_movable_selector() ) );
		ResidueSelectorOP all_chi_move( new OrResidueSelector(cys_design_selector(), chi_movable_selector()) );
		ResidueSelectorOP CYS_chi_movable( new AndResidueSelector(CYS_side_selector, all_chi_move ) );

		// Setup MoveMap
		TR << "Setting up movemap" << std::endl;
		core::select::movemap::MoveMapFactoryOP minmovemapfact( new core::select::movemap::MoveMapFactory );
		minmovemapfact->all_bb( false );
		minmovemapfact->all_chi( false );
		minmovemapfact->add_chi_action( core::select::movemap::mm_enable, CYS_chi_movable );
		if ( disulfides_bb_ ) {
			minmovemapfact->add_bb_action( core::select::movemap::mm_enable, CYS_bb_movable );
		}
		core::kinematics::MoveMapOP minmovemap = minmovemapfact->create_movemap_from_pose( pose );
		minmovemap->show( TR.Trace );

		// Run
		core::util::rebuild_disulfide( pose, nub_->disulfides(), NULL, packer_score, minmovemap, minimizer_score );
		pose.conformation().fix_disulfides( nub_->disulfides() );
		core::pose::setPoseExtraScore( pose, nub_->design_chain() + "_ni_disulfides_side", disulfides_side_ );
		TR.Debug << "REPACK DISULFIDES: END" << std::endl;
	}
}

/// @brief Repacking is necessary to obtain a accetable scored protein.
/// Here it also takes care of designing of changing all residues to x if requested.
void
NubInitioMover::repack( core::pose::Pose & pose )
{

	using namespace core::pack::task::operation;
	using namespace core::pack::task;
	using namespace core::select::residue_selector;

	TR.Debug << "REPACK: START" << std::endl;

	// Setup aligner
	movers::AlignByResidueSelectorMover aligner;
	aligner.reference_pose( pose );
	aligner.reference_selector( work_motif_selector() );

	// Setup Selectors
	TR.Trace << "Setting up selectors" << std::endl;
	ResidueSelectorOP static_selector( new NotResidueSelector( chi_movable_selector() ) );

	// Make task Factory
	TaskFactoryOP postTaskFact( new TaskFactory );
	TaskOperationOP fixMotifOperation( new OperateOnResidueSubset( ResLvlTaskOperationCOP( new PreventRepackingRLT() ), static_selector ) );
	postTaskFact->push_back( fixMotifOperation );
	if ( not design_ ) {
		TR.Trace << "No sequence design" << std::endl;
		TaskOperationOP templateOperation( new OperateOnResidueSubset( ResLvlTaskOperationCOP( new RestrictToRepackingRLT() ), chi_movable_selector() ) );
		postTaskFact->push_back( templateOperation );
	} else {
		if ( residue_type_ == "" ) {
			TR.Trace << "Sequence design allowed" << std::endl;
			DisallowIfNonnativeRLT mutableres;
			mutableres.disallow_aas("CYS");
			TaskOperationOP templateOperation( new OperateOnResidueSubset( mutableres.clone(), chi_movable_selector() ) );
			postTaskFact->push_back( templateOperation );
		} else {
			TR.Trace << "Template residues set to " << residue_type_ << std::endl;
			RestrictAbsentCanonicalAASRLT singleres;
			singleres.aas_to_keep( residue_type_ );
			TaskOperationOP templateOperation( new OperateOnResidueSubset( singleres.clone(), chi_movable_selector() ) );
			postTaskFact->push_back( templateOperation );
		}
	}

	// Repack
	minimization_packing::PackRotamersMover repacker( fullatom_scorefxn_, postTaskFact->create_task_and_apply_taskoperations( pose ), 3 );
	repacker.apply( pose );

	// Realign
	aligner.query_selector( work_motif_selector() );
	aligner.apply( pose );
	TR.Debug << "REPACK: COMPLETE" << std::endl;

}

// -- SETTERS/GETTERS -- //

std::string
NubInitioMover::prefix() const {
	return prefix_;
}

void
NubInitioMover::prefix( std::string const & prefix ) {
	prefix_ = prefix;
}

core::pose::Pose
NubInitioMover::template_pose() const {
	return template_pose_;
}

core::select::residue_selector::ResidueSelectorCOP
NubInitioMover::template_selector() const {
	return template_selector_;
}

void
NubInitioMover::template_selector( core::select::residue_selector::ResidueSelectorCOP const & selector ) {
	if ( selector ) template_selector_ = selector; // Avoid setting up empty selector pointers.
}

core::scoring::ScoreFunctionOP
NubInitioMover::fullatom_scorefxn() const {
	return fullatom_scorefxn_;
}

void
NubInitioMover::fullatom_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ) {
	fullatom_scorefxn_ = scorefxn;
}

utils::NubOP
NubInitioMover::nub() const {
	return nub_;
}

void
NubInitioMover::nub( utils::NubOP nub ) {
	nub_ = nub;
}

bool
NubInitioMover::use_cst() const {
	return use_cst_;
}

void
NubInitioMover::use_cst( bool pick ) {
	use_cst_ = pick;
}

bool
NubInitioMover::clear_motif_cst() const {
	return clear_motif_cst_;
}

void
NubInitioMover::clear_motif_cst( bool pick ) {
	clear_motif_cst_ = pick;
}

core::Real
NubInitioMover::rmsd_threshold() const {
	return rmsd_threshold_;
}

void
NubInitioMover::rmsd_threshold( core::Real threshold ) {
	rmsd_threshold_ = threshold;
}

bool
NubInitioMover::rmsd_include_motif() const {
	return rmsd_include_motif_;
}

void
NubInitioMover::rmsd_include_motif( bool pick ) {
	rmsd_include_motif_ = pick;
}

bool
NubInitioMover::rmsd_include_unconstrained() const {
	return rmsd_include_unconstrained_;
}

void
NubInitioMover::rmsd_include_unconstrained( bool pick ) {
	rmsd_include_unconstrained_ = pick;
}

bool
NubInitioMover::repack_disulfides() const {
	return repack_disulfides_;
}

void
NubInitioMover::repack_disulfides( bool pick ) {
	repack_disulfides_ = pick;
}

bool
NubInitioMover::disulfides_bb() const {
	return disulfides_bb_;
}

void
NubInitioMover::disulfides_bb( bool pick ) {
	disulfides_bb_ = pick;
}

core::Size
NubInitioMover::disulfides_side() const {
	return disulfides_side_;
}
void
NubInitioMover::disulfides_side( core::Real value ) {
	disulfides_side_ = value;
}

std::string
NubInitioMover::fragments_id() const {
	return fragments_id_;
}

void
NubInitioMover::fragments_id( std::string const & name ) {
	fragments_id_ = name;
}

void
NubInitioMover::max_trials( core::Size choice ) {
	max_trials_ = choice;
}

core::Size
NubInitioMover::max_trials() const {
	return max_trials_;
}

void
NubInitioMover::binder_weight( core::Real value ) {
	mc_binder_weight_ = value;
}
core::Real
NubInitioMover::binder_weight() const {
	return mc_binder_weight_;
}

void
NubInitioMover::angle_weight( core::Real value ) {
	mc_angle_weight_ = value;
}

core::Real
NubInitioMover::angle_weight() const {
	return mc_binder_weight_;
}

void
NubInitioMover::dihedral_weight( core::Real value ) {
	mc_dihedral_weight_ = value;
}

core::Real
NubInitioMover::dihedral_weight() const {
	return mc_dihedral_weight_;
}

void
NubInitioMover::correction_weights( core::Real value ) {
	mc_correction_weights_ = value;
}

core::Real
NubInitioMover::correction_weights() const {
	return mc_correction_weights_;
}

bool
NubInitioMover::use_correction_weights() const {
	return use_correction_weights_;
}
void
NubInitioMover::use_correction_weights( bool pick ) {
	use_correction_weights_ = pick;
}

bool
NubInitioMover::dump_centroid() const {
	return dump_centroid_;
}
void
NubInitioMover::dump_centroid( bool pick ) {
	dump_centroid_ = pick;
}

bool
NubInitioMover::drop_unfolded_pose() const {
	return drop_unfolded_pose_;
}
void
NubInitioMover::drop_unfolded_pose( bool pick ) {
	drop_unfolded_pose_ = pick;
}

bool
NubInitioMover::design() const {
	return design_;
}
void
NubInitioMover::design( bool pick ) {
	design_ = pick;
}
std::string
NubInitioMover::residue_type() const {
	return residue_type_;
}
void
NubInitioMover::residue_type( std::string pick ) {
	residue_type_ = pick;
}

// -- COMMON WORKING SELECTORS -- //

core::select::residue_selector::ResidueSelectorOP
NubInitioMover::design_chain_selector() const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP design( new ChainSelector( nub_->design_chain() ) );
	return design;
}

core::select::residue_selector::ResidueSelectorOP
NubInitioMover::context_chain_selector() const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP context( new NotResidueSelector( design_chain_selector() ) );
	return context;
}

core::select::residue_selector::ResidueSelectorOP
NubInitioMover::used_template_selector() const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP noinsert( new NotResidueSelector( template_selector_ ) );
	return noinsert;
}

core::select::residue_selector::ResidueSelectorOP
NubInitioMover::work_motif_selector() const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP label( new ResiduePDBInfoHasLabelSelector( nub_->motif_label() ) );
	ResidueSelectorOP design( new ChainSelector( nub_->design_chain() ) );
	ResidueSelectorOP last_selector( new AndResidueSelector( label, design ) );
	return last_selector;
}

core::select::residue_selector::ResidueSelectorOP
NubInitioMover::work_template_selector() const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP label( new ResiduePDBInfoHasLabelSelector( nub_->template_label() ) );
	ResidueSelectorOP design( new ChainSelector( nub_->design_chain() ) );
	ResidueSelectorOP last_selector( new AndResidueSelector( label, design ) );
	return last_selector;
}

core::select::residue_selector::ResidueSelectorOP
NubInitioMover::bb_movable_selector() const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP label1( new ResiduePDBInfoHasLabelSelector( nub_->template_label() ) );
	ResidueSelectorOP label2( new ResiduePDBInfoHasLabelSelector( nub_->flexible_label() ) );
	ResidueSelectorOP labels( new OrResidueSelector( label1, label2 ) );
	ResidueSelectorOP last_selector( new AndResidueSelector( labels, design_chain_selector() ) );
	return last_selector;
}

core::select::residue_selector::ResidueSelectorOP
NubInitioMover::chi_movable_selector() const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP label1( new ResiduePDBInfoHasLabelSelector( nub_->template_label() ) );
	ResidueSelectorOP label2( new ResiduePDBInfoHasLabelSelector( nub_->coldspot_label() ) );
	ResidueSelectorOP labels( new OrResidueSelector( label1, label2 ) );
	ResidueSelectorOP last_selector( new AndResidueSelector( labels, design_chain_selector() ) );
	return last_selector;
}

core::select::residue_selector::ResidueSelectorOP
NubInitioMover::cys_design_selector() const {
	using namespace core::select::residue_selector;
	ResidueSelectorOP label( new ResiduePDBInfoHasLabelSelector( nub_->disulfide_label() ) );
	ResidueSelectorOP last_selector( new AndResidueSelector( label, design_chain_selector() ) );
	return last_selector;
}

core::select::residue_selector::ResidueRanges
NubInitioMover::template_insertion_ranges() const {
	using namespace core::select::residue_selector;
	ResidueSubset subset = template_selector_->apply( template_pose_ );
	runtime_assert_msg( count_selected( subset ) > 0, "No insertion regions were selected in the template." );
	ResidueRanges ranges( subset );
	return ranges;
}

core::select::residue_selector::ResidueRanges
NubInitioMover::used_template_ranges( bool terminals ) const {
	using namespace core::select::residue_selector;
	ResidueSubset subset = used_template_selector()->apply( template_pose_ );
	TR << represent_residue_selector( subset ) << std::endl;
	TR << represent_residue_selector( template_selector_->apply( template_pose_)) << std::endl;
	runtime_assert_msg( count_selected( subset ) > 0, "No residues were selected to be kept from the template." );
	ResidueRanges ranges( subset );
	if ( terminals ) {
		ResidueRange empty_range( 0, 0 );
		ResidueRanges insert_ranges = template_insertion_ranges();
		TR.Trace << "r1: " << insert_ranges[1].start();
		TR.Trace << " r2: " << insert_ranges[ insert_ranges.size() ].stop();
		TR.Trace << " of " << insert_ranges.size() << " pieces ";
		TR.Trace << " len: " << template_pose_.size() << std::endl;
		// Do we have a N-terminal insertion? Add an empty range at the beginning of the range list!
		runtime_assert_msg( insert_ranges[1].start() <= template_pose_.size(),
			"Insertion ranges start out of the template's length!"  );
		if ( insert_ranges[1].start() == 1 ) {
			TR << "  N-terminal insertion found!" << std::endl;
			ranges.insert( ranges.begin(), empty_range );
		}
		// Do we have a C-terminal insertion? Add an empty range at the end of the range list!
		runtime_assert_msg( insert_ranges[ insert_ranges.size() ].stop() <= template_pose_.size(),
			"Insertion ranges finish out of the template's length!"  );
		if ( insert_ranges[ insert_ranges.size() ].stop() == template_pose_.size() ) {
			TR << "  C-terminal insertion found!" << std::endl;
			ranges.push_back( empty_range );
		}
	}
	return ranges;
}

// -- DEFAULTS -- //

std::string
NubInitioMover::default_prefix() {
	return "FFL";
}

core::select::residue_selector::ResidueSelectorCOP
NubInitioMover::default_template_selector() {
	using namespace core::select::residue_selector;
	ResidueSelectorCOP true_selector( new TrueResidueSelector );
	return true_selector;
}

core::scoring::ScoreFunctionOP
NubInitioMover::default_fullatom_scorefxn() {
	return core::scoring::get_score_function();
}

bool
NubInitioMover::default_use_cst() {
	return true;
}

bool
NubInitioMover::default_clear_motif_cst() {
	return true;
}

core::Real
NubInitioMover::default_rmsd_threshold() {
	return 5.0;
}

bool
NubInitioMover::default_rmsd_include_motif() {
	return false;
}

bool
NubInitioMover::default_rmsd_include_unconstrained() {
	return true;
}

bool
NubInitioMover::default_repack_disulfides() {
	return true;
}

bool
NubInitioMover::default_disulfides_bb() {
	return false;
}

core::Size
NubInitioMover::default_disulfides_side() {
	return 0;
}

core::Size
NubInitioMover::default_max_trials() {
	return 10;
}

core::Real
NubInitioMover::default_mc_binder_weight() {
	return 0.0;
}

core::Real
NubInitioMover::default_mc_angle_weight() {
	return 0.0;
}

core::Real
NubInitioMover::default_mc_dihedral_weight() {
	return 0.0;
}

core::Real
NubInitioMover::default_mc_correction_weights() {
	return 1.17;
}

bool
NubInitioMover::default_use_correction_weights() {
	return true;
}

bool
NubInitioMover::default_dump_centroid() {
	return false;
}

bool
NubInitioMover::default_drop_unfolded_pose() {
	return false;
}

bool
NubInitioMover::default_design() {
	return false;
}

std::string
NubInitioMover::default_residue_type() {
	return std::string();
}

// -- ROSETTASCRIPTS -- //

protocols::moves::MoverOP
NubInitioMover::fresh_instance() const {
	return protocols::moves::MoverOP( new NubInitioMover );
}

protocols::moves::MoverOP
NubInitioMover::clone() const {
	return protocols::moves::MoverOP( new NubInitioMover( *this ) );
}

std::string
NubInitioMover::get_name() const {
	return mover_name();
}

std::string
NubInitioMover::mover_name() {
	return "NubInitioMover";
}

// -- MOVERCREATOR -- //

moves::MoverOP
NubInitioMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new NubInitioMover );
}

std::string
NubInitioMoverCreator::keyname() const {
	return NubInitioMover::mover_name();
}

void
NubInitioMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	NubInitioMover::provide_xml_schema( xsd );
}

}
}
