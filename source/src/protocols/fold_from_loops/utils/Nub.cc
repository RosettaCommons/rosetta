// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   Nub.cc
/// @brief
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

// Unit headers
#include <protocols/fold_from_loops/utils/Nub.hh>
#include <protocols/fold_from_loops/utils/utils.hh>
#include <protocols/fold_from_loops/movers/SplitAndMixPoseMover.hh>
#include <protocols/fold_from_loops/movers/LabelPoseFromResidueSelectorMover.hh>

// Protocols headers
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/moves/DsspMover.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/chains_util.hh>
#include <core/conformation/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIterator.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/AndResidueSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ResiduePDBInfoHasLabelSelector.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

#include <string>
#include <algorithm>

namespace protocols {
namespace fold_from_loops {
namespace utils {

static basic::Tracer TR( "protocols.fold_from_loops.Nub", basic::t_debug );

/// @brief Empty Constructor
Nub::Nub():
	reference_pose_( new core::pose::Pose ),
	unfolded_pose_( new core::pose::Pose ),
	small_( new core::fragment::ConstantLengthFragSet ),
	large_( new core::fragment::ConstantLengthFragSet ),
	selector_( default_selector() ),
	binder_selector_( default_selector() ),
	insertion_selector_( default_selector() ),
	unfolded_jump_( 0 ),
	has_binder_( false ),
	movemap_( new core::select::movemap::MoveMapFactory ),
	design_size_( 0 ),
	working_cst_( new core::scoring::constraints::ConstraintSet ),
	default_flexibility_( 0 ),
	nub_disulfides_( 0 )
{
	movemap_->all_bb( false );
	movemap_->all_chi( false );
	movemap_->all_nu( false );
	movemap_->all_jumps( false );
}

/// @brief Destructor
Nub::~Nub(){}

// -- APPLY -- //

/// Produces the unfolded pose
void
Nub::apply(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & template_selection,
	utility::vector1< core::pose::PoseOP > const & template_pieces,
	utility::vector1< std::pair< core::Size, core::Size > > const & disulfides )
{

	core::pose::Pose work_pose;
	work_pose.detached_copy( *reference_pose_ );
	runtime_assert_msg( !work_pose.empty(), "The provided pose is empty. Check the given PDB file or Reference Pose." );

	TR << "Assigning Secondary Structure to Pose with DSSP." << std::endl;
	protocols::moves::DsspMover dssp; // Secondary Structure Assignment.
	dssp.apply( work_pose );

	TR << "Obtaining Nub regions of interest." << std::endl;
	get_nub_pieces( work_pose );

	TR << "Creating the unfolded pose." << std::endl;
	join_pieces( template_pieces );

	TR << "Mapping unfolded pose to template." << std::endl;
	seqmap_ = utils::map_by_residue_subsets( pose, template_selection, *unfolded_pose_, insertion_selector_->apply( *unfolded_pose_ ) );

	TR << "Fixing Disulfides" << std::endl;
	assign_disulfides( disulfides );

	TR << "Fixing fragments if template's size/disposition has changed." << std::endl;
	fix_fragments();

	TR << "Adding the binders (if any)." << std::endl;
	add_binders();

	TR << "Label the residues if the unfolded pose." << std::endl;
	add_labels();

	TR << "Creating the Movemap" << std::endl;
	make_movemap();
}

// -- ROSETTASCRIPTS -- //

void
Nub::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &data, core::pose::Pose const & reference_pose )
{

	TR.Trace << TR.Green << "Parsing Nub data" << TR.Reset << std::endl;
	utility::tag::TagCOP const& nubtag = tag->getTag( object_name() );

	if ( nubtag->hasOption( "reference_name" ) ) {
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose(nubtag, data, "reference_name" );
		TR.Trace << TR.Green << "Loaded reference pose: " << nubtag->getOption< std::string >( "reference_name" ) << " with " << reference_pose_->size() << " residues" << TR.Reset << std::endl;
	} else {
		reference_pose_ = core::pose::PoseOP( new core::pose::Pose( reference_pose ) );
		if ( nubtag->hasOption( "pose_file" ) ) {
			core::import_pose::pose_from_file( *reference_pose_, nubtag->getOption< std::string >( "pose_file" ) , core::import_pose::PDB_file);
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "No pose has been provided from which to extract the functional motif.");
		}
	}

	TR.Trace << TR.Green << "Getting Nub selection" << TR.Reset << std::endl;
	selector( core::select::residue_selector::parse_residue_selector( nubtag, data, "residue_selector" ) );
	if ( nubtag->hasOption( "binder_selector" ) ) {
		TR.Trace << TR.Green << "Getting Binder(s) selection" << TR.Reset << std::endl;
		binder_selector( core::select::residue_selector::parse_residue_selector( nubtag, data, "binder_selector" ) );
	} else {
		TR.Trace << TR.Green << "No binder selected" << TR.Reset << std::endl;
	}

	// Process Motifs Properties (if any)
	utility::vector1< utility::tag::TagCOP > const & segments = nubtag->getTags( segment_object_name() );
	if ( segments.size() != 0 ) {
		for ( auto segtag : segments ) {
			NubSegment seg;
			seg.order  = segtag->getOption< core::Size >( "order" );
			seg.c_flex = segtag->getOption< core::Size >( "c_term_flex", default_flexibility_ );
			seg.n_flex = segtag->getOption< core::Size >( "n_term_flex", default_flexibility_ );
			if ( segtag->hasOption( "editable" ) ) {
				seg.coldspots = utility::string_split( segtag->getOption< std::string >( "editable", "" ), ',', core::Size() );
			}
			piece_properties_.push_back( seg );
		}
	}
}

void
Nub::provide_xml_definition( utility::tag::XMLSchemaDefinition & xsd, utility::tag::XMLSchemaSimpleSubelementList & elements )
{
	using namespace utility::tag;
	XMLSchemaSimpleSubelementList subelements;
	AttributeList segattr;
	segattr
		+ XMLSchemaAttribute::required_attribute(  "order", xs_integer, "Specifies to which motif segment this definition is targetting" )
		+ XMLSchemaAttribute::attribute_w_default( "n_term_flex", xs_integer, "Number of residues in the N-terminus of the motif segment allowed to move", "0" )
		+ XMLSchemaAttribute::attribute_w_default( "c_term_flex", xs_integer, "Number of residues in the C-terminus of the motif segment allowed to move", "0" )
		+ XMLSchemaAttribute( "editable", xsct_int_cslist, "Residues in the motif allowed to mutate (COLDSPOTS)" );
	subelements.add_simple_subelement( segment_object_name() , segattr , "Describes the properties of a given segment of the nub motif" );
	subelements.complex_type_naming_func( & complex_object_name );


	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "pose_file", xs_string,
		"File containing the structure with the static segment to keep (Nub). This option will be ignored if a reference_pose "
		"loaded with SavePoseMover is provided." );

	protocols::rosetta_scripts::attributes_for_saved_reference_pose_w_description( attlist,
		"Reference pose containing the structure with the static segment to keep (Nub). This option will "
		"override the pose_file option.", "reference_name" );
	core::select::residue_selector::attributes_for_parse_residue_selector_when_required( attlist, "residue_selector",
		"Selector specifying residues that will be kept as the Nub structure. "
		"The number of segments this selector specifies has to be the same as the number of motifs that "
		"the template selector will specify in the NubInitio Mover, but the sizes do not have to be the same.");
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "binder_selector",
		"If provided, it selects a binder to be placed together with the region of interest in the Nub. "
		"The binder(s) will be static during the Ab initio process.");


	XMLSchemaComplexTypeGenerator nub_ctgen;
	nub_ctgen.element_name( object_name() )
		.add_attributes( attlist )
		.complex_type_naming_func( & complex_object_name )
		.set_subelements_repeatable( subelements )
		.description( "Describes the static segment in the NubInitio protocol" )
		.write_complex_type_to_schema( xsd );

	elements.add_already_defined_subelement( object_name(), & complex_object_name );

}

// -- SETTERS/GETTERS, STATICS AND UTILITY FUNCTIONS -- //

void
Nub::refit_motif_sidechains( core::pose::Pose & pose )
{
	using namespace core::select::residue_selector; //ResiduePDBInfoHasLabelSelector
	using namespace core::chemical;
	ResiduePDBInfoHasLabelSelectorOP selector( new ResiduePDBInfoHasLabelSelector( motif_label() ) );
	ResidueSubset motifs = selector->apply( pose );

	for ( core::Size i = 1; i <= motifs.size(); ++i ) {
		if ( motifs[i] ) {
			pose.replace_residue( i, unfolded_pose_->residue( i ), true );
			// if ( unfolded_pose_->conformation().residue( i ).has_variant_type( CUTPOINT_UPPER ) ) {
			//  core::pose::remove_variant_type_from_pose_residue( pose, CUTPOINT_UPPER, i );
			// }
		}
	}
	if ( has_binder_ ) {
		ResidueSelectorCOP design_chain_selector( new ChainSelector( design_chain_ ) );
		ResidueSelectorCOP bindselector( new NotResidueSelector( design_chain_selector ) );
		ResidueSubset binders = bindselector->apply( pose );
		for ( core::Size i = 1; i <= binders.size(); ++i ) {
			if ( binders[i] ) {
				pose.replace_residue( i, unfolded_pose_->residue( i ), true );
			}
		}
	}

}

void
Nub::assign_disulfides( utility::vector1< std::pair< core::Size, core::Size > > const & disulfides )
{
	using namespace core::select::residue_selector;
	for ( auto disf : disulfides ) {
		if ( seqmap_[ disf.first ] != 0 and seqmap_[ disf.second ] != 0 ) {
			if ( unfolded_pose_->conformation().residue( seqmap_[ disf.first ] ).name3() == "CYS" and
					unfolded_pose_->conformation().residue( seqmap_[ disf.second ] ).name3() == "CYS" ) {
				core::pose::add_variant_type_to_pose_residue( *unfolded_pose_, core::chemical::DISULFIDE, seqmap_[ disf.first ] );
				core::pose::add_variant_type_to_pose_residue( *unfolded_pose_, core::chemical::DISULFIDE, seqmap_[ disf.second ] );
				std::pair< core::Size, core::Size > thispair( seqmap_[ disf.first ], seqmap_[ disf.second ] );
				disulfides_.push_back( thispair );
			}
		}
	}

	if ( nub_disulfides_ > 0 ) {
		movers::LabelPoseFromResidueSelectorMover labeler;
		ResiduePDBInfoHasLabelSelector labsel;
		for ( core::Size i = 0; i < nub_disulfides_; ++i ) {
			labsel.set_label( "DIS" + std::to_string(i) + "NB" );
			utility::vector1< core::Size > picks = selection_positions( labsel.apply(*unfolded_pose_) );
			if ( picks.size() > 1 ) { // only for nub to nub disulfides.
				core::pose::add_variant_type_to_pose_residue( *unfolded_pose_, core::chemical::DISULFIDE, picks[1] );
				core::pose::add_variant_type_to_pose_residue( *unfolded_pose_, core::chemical::DISULFIDE, picks[2] );
				std::pair< core::Size, core::Size > thispair( picks[1], picks[2] );
				disulfides_.push_back( thispair );
			}
			std::string label2remove = "DIS" + std::to_string(i) + "NB";
			labeler.unlabel( label2remove );
			labeler.residue_selector( labsel );
			labeler.apply( *unfolded_pose_ );
		}
	}

	for ( auto disf: disulfides_ ) {
		TR << disf.first << "-" << disf.second << std::endl;
		unfolded_pose_->pdb_info()->add_reslabel( disf.first, disulfide_label() );
		unfolded_pose_->pdb_info()->add_reslabel( disf.second, disulfide_label() );
	}

	nub_disulfides_ = 0;
}

void
Nub::transfer_unfolded_conformation( core::pose::Pose & pose )
{
	using namespace core::scoring::constraints;
	TR << pose.size() << std::endl;
	TR << unfolded_pose_->annotated_sequence() << std::endl;
	TR << unfolded_pose_->secstruct() << std::endl;

	if ( working_cst_->is_empty() ) {
		ConstraintSetOP thiscst(new core::scoring::constraints::ConstraintSet( *pose.constraint_set() ));
		*working_cst_ = *thiscst;
	}
	// unfolded_pose_->transfer_constraint_set( pose );
	pose.clear();
	TR << "Transfer Confromation" << std::endl;
	pose.set_new_conformation( unfolded_pose_->conformation_ptr() );
	pose.conformation().detect_disulfides();
	pose.conformation().fix_disulfides( disulfides_ );
	TR << "Transfer FoldTree" << std::endl;
	pose.fold_tree( unfolded_pose_->fold_tree() );
	TR << "Transfer PDBInfo" << std::endl;
	pose.pdb_info( unfolded_pose_->pdb_info() );
	TR << "Transfer Constraints" << std::endl;
	// pose.transfer_constraint_set( *unfolded_pose_ );
	pose.constraint_set( working_cst_ );

}

void
Nub::get_nub_pieces( core::pose::Pose & pose )
{
	using namespace core::select::residue_selector; // ResidueSubset, ResidueRanges
	ResidueSubset segments = selector_->apply( pose );
	ResidueRanges ranges( segments );
	runtime_assert_msg( has_any_true_selection( segments ), "There is no segment in the Reference Pose selected to be kept!!" );
	if ( TR.visible() ) {
		TR << "Regions of interest are: " << std::endl;
		for ( auto const & range : ranges ) {
			TR << "  N >  " << range.to_string() << std::endl;
		}
	}

	utility::vector1< std::pair< core::Size, core::Size > > nub_disulfides;
	core::conformation::disulfide_bonds( pose.conformation(), nub_disulfides );
	for ( auto disf: nub_disulfides ) {
		TR << disf.first << "-" << disf.second << std::endl;
		pose.pdb_info()->add_reslabel( disf.first,  "DIS" + std::to_string(nub_disulfides_) + "NB" );
		pose.pdb_info()->add_reslabel( disf.second, "DIS" + std::to_string(nub_disulfides_) + "NB" );
		nub_disulfides_ += 1;
	}

	movers::SplitAndMixPoseMover splitter;
	splitter.set_residue_selector( selector_ );
	if ( piece_properties_.size() > 0 ) {
		utility::vector1< core::Size > order = utility::vector1< core::Size >();
		for ( core::Size i = 1; i <= piece_properties_.size(); ++i ) {
			order.push_back( piece_properties_[i].order );
		}
		splitter.set_order( order );
	}
	pieces_ = splitter.apply_without_merge( pose );
	design_chain_ = core::pose::get_chain_from_chain_id( 1, *pieces_[1] );
	TR << "Design chain is: " << design_chain_ << std::endl;

	if ( default_flexibility_ != 0 and piece_properties_.size() == 0 ) {
		for ( core::Size p = 1; p <= pieces_.size(); ++p ) {
			NubSegment seg;
			seg.order  = p;
			seg.c_flex = default_flexibility_;
			seg.n_flex = default_flexibility_;
			piece_properties_.push_back( seg );
		}
	}

	if ( piece_properties_.size() > 0 ) {
		runtime_assert_msg( piece_properties_.size() == pieces_.size(),
			"Specified number of Segments should be the same as motif fragments. ");
		movers::LabelPoseFromResidueSelectorMover labeler;
		for ( core::Size p = 1; p <= piece_properties_.size(); ++p ) {
			// FLEXIBILITY
			core::Size n = std::min(piece_properties_[p].n_flex, pieces_[p]->size() );
			for ( core::Size i = 1; i <= n; ++i ) {
				pieces_[p]->pdb_info()->add_reslabel(i, flexible_label() );
			}
			core::Size c = std::min(piece_properties_[p].c_flex, pieces_[p]->size() );
			for ( core::Size i = pieces_[p]->size() - c + 1; i <= pieces_[p]->size(); ++i ) {
				pieces_[p]->pdb_info()->add_reslabel(i, flexible_label() );
			}
			// EDITABILITY
			ResidueIndexSelector colds;
			for ( core::Size i = 1; i <=  piece_properties_[p].coldspots.size(); ++i ) {
				colds.append_index( piece_properties_[p].coldspots[i] );
			}
			ResidueSelectorCOP coldscop( new ResidueIndexSelector( colds ) );
			ResidueSelectorCOP hots( new NotResidueSelector( coldscop ) );
			labeler.residue_selector( coldscop );
			labeler.label( coldspot_label() );
			labeler.apply( *pieces_[p] );
			labeler.residue_selector( hots );
			labeler.label( hotspot_label() );
			labeler.apply( *pieces_[p] );
		}
	} else {
		for ( core::Size p = 1; p <= pieces_.size(); ++p ) {
			for ( core::Size i = 1; i <= pieces_[p]->size(); ++i ) {
				pieces_[p]->pdb_info()->add_reslabel(i, hotspot_label() );
			}
		}
	}
}

// @brief Attaches the template unfolded regions to the Nub to create the unfolded pose.
void
Nub::join_pieces( utility::vector1< core::pose::PoseOP > const & template_pieces )
{
	using namespace core::select::residue_selector;
	using namespace core::chemical; // UPPER/LOWER_TERMINUS_VARIANT, CUTPOINT_LOWER/UPPER
	runtime_assert_msg( ( template_pieces.size() == ( pieces_.size() * 2 ) ),
		"The number of template pieces ( " + std::to_string( template_pieces.size() ) + " ) "
		"has to be twice the number of Nub regions ( " + std::to_string( pieces_.size() * 2 ) + " )" );

	core::Size length = 0;
	std::string str_ranges;
	for ( core::Size i=1; i<=pieces_.size(); ++i ) {
		core::Size j = i * 2;
		length = template_pieces[ j - 1 ]->size() + length;
		ResidueRange current_range( length + 1, length + pieces_[i]->size() );
		if ( TR.Trace.visible() ) {
			TR.Trace << TR.Green << "Attaching to static motif " << i << std::endl;
			TR.Trace << TR.Green << "   N-terminal piece of length " << template_pieces[ j - 1 ]->size() << std::endl;
			if ( template_pieces[ j - 1 ]->size() > 0 ) {
				TR.Trace << TR.Green << "   Sequence: " << template_pieces[ j - 1 ]->sequence() << std::endl;
			}
			TR.Trace << TR.Green << "   C-terminal piece of length " << template_pieces[ j ]->size() << std::endl;
			if ( template_pieces[ j ]->size() ) {
				TR.Trace << TR.Green << "   Sequence " << template_pieces[ j ]->sequence() << std::endl;
			}
		}
		utils::attach_n_and_c_unfolded_poses_to_pose( *template_pieces[ j - 1 ], *pieces_[ i ], *template_pieces[ j ] );
		if ( i != 1 ) {
			core::pose::add_variant_type_to_pose_residue( *pieces_[ i ], CUTPOINT_UPPER, 1 );
		}
		if ( i != pieces_.size() ) {
			core::pose::add_variant_type_to_pose_residue( *pieces_[ i ], CUTPOINT_LOWER, pieces_[ i ]->size() );
		}
		utils::append_pose_to_pose_keep_fold_tree( *unfolded_pose_, *pieces_[i], false );
		length = unfolded_pose_->size();
		str_ranges.append( current_range.to_string() );
		if ( i < pieces_.size() ) {
			str_ranges.append( "," );
		}
	}
	TR << "Final ranges of the nub are: " << str_ranges << std::endl;

	core::pose::remove_variant_type_from_pose_residue( *unfolded_pose_, CUTPOINT_UPPER, 1 );
	core::pose::add_variant_type_to_pose_residue( *unfolded_pose_, LOWER_TERMINUS_VARIANT, 1);
	core::pose::remove_variant_type_from_pose_residue( *unfolded_pose_, CUTPOINT_LOWER, unfolded_pose_->size() );
	core::pose::add_variant_type_to_pose_residue( *unfolded_pose_, UPPER_TERMINUS_VARIANT, unfolded_pose_->size() );
	core::pose::correctly_add_cutpoint_variants( *unfolded_pose_ );

	ResidueSelectorCOP index_selector( new ResidueIndexSelector( str_ranges ) );
	insertion_selector_ = index_selector;
	unfolded_jump_   = unfolded_pose_->fold_tree().num_jump();
	design_size_    = unfolded_pose_->size();
	unfolded_pose_->pdb_info()->obsolete(false);
	unfolded_pose_->pdb_info()->rebuild_pdb2pose();
	for ( core::Size i = 1; i <= unfolded_pose_->size(); i++ ) {
		utility::vector1< std::string > labels = unfolded_pose_->pdb_info()->get_reslabels( i );
		unfolded_pose_->pdb_info()->set_resinfo( i, *design_chain_.c_str(), i );
		for ( auto label: labels ) {
			unfolded_pose_->pdb_info()->add_reslabel( i, label );
		}
	}
}

/// @brief Add selecter binders if provided.
void
Nub::add_binders()
{
	using namespace core::select::residue_selector;
	using namespace core::kinematics;
	protocols::moves::DsspMover dssp;
	ResidueSubset binder_subset = binder_selector_->apply( *reference_pose_ );
	if ( core::select::residue_selector::has_any_true_selection( binder_subset ) ) {
		core::pose::Pose binder;
		binder.detached_copy( *reference_pose_ );
		binder.remove_constraints();
		protocols::grafting::simple_movers::DeleteRegionMover deleter;
		ResidueSelectorCOP notindex( new NotResidueSelector( binder_selector_ ) );
		deleter.set_residue_selector( notindex );
		deleter.apply( binder );
		utility::vector1< core::pose::PoseOP > binders = binder.split_by_chain();
		core::Size current_length = unfolded_pose_->size();
		for ( auto const & bind : binders ) {
			FoldTree bind_tree = FoldTree( bind->size() );
			TR.Trace << "ADD BINDER TREE " << bind_tree  << std::endl;
			core::Size closest = closest_binder( *bind );
			TR.Trace << "SPLIT TREE AT   " << closest  << std::endl;
			if ( closest != 1 and closest != bind->size() ) {
				bind_tree.split_existing_edge_at_residue( closest );
				bind_tree.reorder( closest );
			}
			bind->fold_tree( bind_tree );
			TR.Trace <<"ADD BINDER TREE " << bind->fold_tree() << std::endl;
			core::pose::add_lower_terminus_type_to_pose_residue ( *bind, 1 );
			core::pose::add_upper_terminus_type_to_pose_residue ( *bind , bind->size() );
			dssp.apply( *bind );
			utils::append_pose_to_pose_keep_fold_tree( *unfolded_pose_, *bind, true );
			unfolded_pose_->pdb_info()->copy(*(bind->pdb_info()), 1, bind->size(), current_length + 1);
			char bind_chain = bind->pdb_info()->chain(1);
			for ( core::Size i = 1; i <= bind->size(); ++i ) {
				// Renumber chain from 1, as silent_files don't play well with insertion codes
				unfolded_pose_->pdb_info()->set_resinfo( current_length + i, bind_chain, i);
			}
			current_length = unfolded_pose_->size();
		}
		has_binder_ = true;
		unfolded_pose_->pdb_info()->obsolete(false);
	} else {
		TR << "No binders were requested" << std::endl;
	}
}

/// @brief Get closest binder residue to the FoldTree root
core::Size
Nub::closest_binder( core::pose::Pose const & pose )
{
	core::Size root = unfolded_pose_->fold_tree().root();
	core::Real distance = 1000;
	core::Size residue = 0;
	for ( core::Size i = 1; i <= pose.size(); i++ ) {
		core::Real const current = unfolded_pose_->residue( root ).nbr_atom_xyz().distance(pose.residue(i).nbr_atom_xyz()) ;
		if ( current < distance ) {
			distance = current;
			residue = i;
		}
	}
	return residue;
}

/// @brief Add labels to the unfolded pose
void
Nub::add_labels()
{
	using namespace core::select::residue_selector;
	movers::LabelPoseFromResidueSelectorMover labeler;
	// 1. Label NUB region
	labeler.residue_selector( insertion_selector_ );
	labeler.label( motif_label() );
	labeler.apply( *unfolded_pose_ );

	// 2. Label Binder (if any)
	if ( has_binder_ ) {
		ResidueSelectorCOP design_chain_selector( new ChainSelector( design_chain_ ) );
		ResidueSelectorCOP binder( new NotResidueSelector( design_chain_selector ) );
		labeler.residue_selector( binder );
		labeler.label( context_label() );
		labeler.apply( *unfolded_pose_ );
	}

	// 3. Label Template
	ResidueSelectorCOP scaffold( new ChainSelector( design_chain_ ) );
	ResidueSelectorCOP notinsert( new NotResidueSelector( insertion_selector_ ) );
	ResidueSelectorCOP scaffnoins( new AndResidueSelector( scaffold, notinsert ) );
	labeler.residue_selector( scaffnoins );
	labeler.label( template_label() );
	labeler.apply( *unfolded_pose_ );
}

/// @brief Creates the movemap to submit to ab initio protocol.
void
Nub::make_movemap()
{
	using namespace core::select::residue_selector;
	ResidueSelectorCOP yes_template( new ResiduePDBInfoHasLabelSelector( template_label() ) );
	movemap_->add_bb_action( core::select::movemap::mm_enable, yes_template );
	movemap_->add_chi_action( core::select::movemap::mm_enable, yes_template );
	ResidueSelectorCOP flex_template( new ResiduePDBInfoHasLabelSelector( flexible_label() ) );
	movemap_->add_bb_action( core::select::movemap::mm_enable, flex_template );
}

/// @brief Fixes fragments in case the insertion has changed the count of the residues of the template.
void
Nub::fix_fragments()
{
	runtime_assert_msg( small_->size() > 0, "FragSet for small fragments provided is empty");
	runtime_assert_msg( large_->size() > 0, "FragSet for large fragments provided is empty");
	core::select::residue_selector::ResidueSubset unfolded_selection = insertion_selector_->apply( *unfolded_pose_ );
	if ( not seqmap_.is_identity() ) {
		small_ = fix_fragment( small_ );
		large_ = fix_fragment( large_ );
	} else {
		TR << "Global Size and disposition of template residues has not changed. Fragments are OK." << std::endl;
	}
}

/// @brief Fixes a FragSet
core::fragment::FragSetOP
Nub::fix_fragment( core::fragment::FragSetOP fragset )
{
	using namespace core::fragment;
	FragSetOP myFragSet = FragSetOP( new ConstantLengthFragSet( fragset->max_frag_length() ) );
	core::Size i = 1;
	for ( FrameIterator it=fragset->nonconst_begin(), eit=fragset->nonconst_end(); it!=eit; ++it, ++i ) {
		TR << i << " " << seqmap_[ i ];  // Size of the map == Size of fset
		if ( seqmap_[ i ] != 0 ) {
			FrameOP frame( new Frame( seqmap_[ i ] ) );
			TR << " new frame created with " << (*it)->nr_frags() << std::endl;
			for ( core::Size j = 1; j<= (*it)->nr_frags(); ++j ) {
				FragData old_fragment = (*it)->fragment( j );
				FragDataOP current_fragment( NULL );
				current_fragment = FragDataOP( new AnnotatedFragData( old_fragment.pdbid(), seqmap_[ i ] ) );
				current_fragment->copy( (*it)->fragment( j ) );
				frame->add_fragment( current_fragment );
			}
			myFragSet->add( frame );
		} else if ( seqmap_[ i ] == 0 ) {
			TR << " skip!" << std::endl;
		}
	}
	TR << "original fragset length: " << fragset->max_pos() << std::endl;
	TR << "new fragset length: " << myFragSet->max_pos() << std::endl;
	return myFragSet;
}

core::pose::PoseOP
Nub::reference_pose() const
{
	return reference_pose_;
}

void
Nub::reference_pose( core::pose::PoseOP const & pose )
{
	reference_pose_->detached_copy( *pose );
}

core::pose::PoseOP
Nub::unfolded_pose() const
{
	return unfolded_pose_;
}

core::fragment::FragSetOP
Nub::small_fragments() const
{
	return small_;
}

void
Nub::small_fragments( core::fragment::FragSetOP fragset )
{
	small_ = fragset;
}

core::fragment::FragSetOP
Nub::large_fragments() const
{
	return large_;
}

void
Nub::large_fragments( core::fragment::FragSetOP fragset )
{
	large_ = fragset;
}

core::select::residue_selector::ResidueSelectorCOP
Nub::selector() const
{
	return selector_;
}

core::id::SequenceMapping
Nub::template_to_unfolded_mapping() const
{
	return seqmap_;
}

void
Nub::selector( core::select::residue_selector::ResidueSelectorCOP const & selector )
{
	selector_ = selector;
}

core::select::residue_selector::ResidueSelectorCOP
Nub::binder_selector() const
{
	return binder_selector_;
}

void
Nub::binder_selector( core::select::residue_selector::ResidueSelectorCOP const & selector )
{
	binder_selector_ = selector;
}

bool
Nub::has_binder() const
{
	return has_binder_;
}

bool
Nub::is_multisegment() const
{
	return pieces_.size() > 1;
}

core::Size
Nub::unfolded_jump() const
{
	return unfolded_jump_;
}

core::Size
Nub::design_size() const
{
	return design_size_;
}

utility::vector1< core::pose::PoseOP >
Nub::pieces() const
{
	return pieces_;
}

core::Size
Nub::default_flexibility() const {
	return default_flexibility_;
}
void
Nub::default_flexibility( core::Size pick ) {
	default_flexibility_ = pick;
}

core::select::residue_selector::ResidueSubset
Nub::apply_selector( core::pose::Pose const & pose ) const
{
	return selector_->apply( pose );
}

core::select::residue_selector::ResidueRanges
Nub::apply_selector_as_ranges( core::pose::Pose const & pose ) const
{
	core::select::residue_selector::ResidueRanges result( selector_->apply( pose ) );
	return result;
}

core::select::residue_selector::ResidueSubset
Nub::apply_binder_selector( core::pose::Pose const & pose ) const
{
	return binder_selector_->apply( pose );
}

core::select::residue_selector::ResidueRanges
Nub::apply_binder_selector_as_ranges( core::pose::Pose const & pose ) const
{
	core::select::residue_selector::ResidueRanges result( binder_selector_->apply( pose ) );
	return result;
}

core::kinematics::MoveMapOP
Nub::movemap() const
{
	return movemap_->create_movemap_from_pose( *unfolded_pose_ );
}
core::select::movemap::MoveMapFactoryOP
Nub::movemapfactory() const
{
	using namespace core::select::movemap;
	return MoveMapFactoryOP( new MoveMapFactory( *movemap_ ) );
}

std::string
Nub::design_chain() const
{
	return design_chain_;
}

utility::vector1< std::pair< core::Size, core::Size > >
Nub::disulfides() const
{
	return disulfides_;
}

std::string
Nub::object_name()
{
	return "Nub";
}

std::string
Nub::segment_object_name()
{
	return "Segment";
}

std::string
Nub::complex_object_name( std::string tag_name )
{
	return "nub_initio_" + tag_name + "_complex_type";
}

core::select::residue_selector::ResidueSelectorCOP
Nub::default_selector()
{
	using namespace core::select::residue_selector;
	ResidueSelectorCOP true_selector( new TrueResidueSelector );
	ResidueSelectorCOP false_selector( new NotResidueSelector( true_selector) );
	return false_selector;
}

std::string
Nub::motif_label()
{
	return "MOTIF";
}

std::string
Nub::flexible_label()
{
	return "FLEXIBLE";
}

std::string
Nub::template_label()
{
	return "TEMPLATE";
}

std::string
Nub::hotspot_label()
{
	return "HOTSPOT";
}
std::string
Nub::coldspot_label()
{
	return "COLDSPOT";
}
std::string
Nub::context_label()
{
	return "CONTEXT";
}

std::string
Nub::disulfide_label()
{
	return "DISULFIDIZE";
}

std::string
Nub::contact_label()
{
	return "CONTACT";
}

}
}
}
