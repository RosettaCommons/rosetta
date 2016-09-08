// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/RemodelLoopMoverPoseFolder.cc
/// @brief Folds residues in a pose using RemodelLoopMover
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/components/RemodelLoopMoverPoseFolder.hh>

// Protocol headers
#include <protocols/denovo_design/components/Picker.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/forge/remodel/RemodelLoopMover.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/conformation/Conformation.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// Basic/Utililty headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.components.RemodelLoopMoverPoseFolder" );

namespace protocols {
namespace denovo_design {
namespace components {

RemodelLoopMoverPoseFolder::RemodelLoopMoverPoseFolder():
	PoseFolder( RemodelLoopMoverPoseFolder::class_name() ),
	scorefxn_()
{
}

RemodelLoopMoverPoseFolder::~RemodelLoopMoverPoseFolder()
{}

RemodelLoopMoverPoseFolder::PoseFolderOP
RemodelLoopMoverPoseFolder::clone() const
{
	return RemodelLoopMoverPoseFolder::PoseFolderOP( new RemodelLoopMoverPoseFolder( *this ) );
}

std::string
RemodelLoopMoverPoseFolder::class_name()
{
	return "RemodelLoopMoverPoseFolder";
}

void
RemodelLoopMoverPoseFolder::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
}

/// @brief performs folding
/// @param pose    - The pose to be folded, with all residues added.  The pose should be prepared with
///                  any necessary cutpoints added before giving to the PoseFolder. Torsions in the pose
///                  should be adjusted, and no residues should be added or removed.
/// @param movable - Subset of residues for which new backbone conformations will be sampled. Residues
///                  specified as 'True' in movable must also be present in one or more Loops in order
///                  to be folded. Movable's size must match pose.size()
/// @param loops   - Loops to be folded.  Cutpoints specified here must be match the cutpoints found in
///                  the pose. Residues not within any loop should not be folded. Residues contained
///                  in a loop but not in the movable set should not be folded.
/// @throws EXCN_Fold if anything goes wrong in folding. Derived classes should throw this.
void
RemodelLoopMoverPoseFolder::apply(
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueSubset const & movable,
	protocols::loops::Loops const & loops ) const
{
	components::StructureData sd = StructureDataFactory::get_instance()->get_from_pose( pose );
	protocols::moves::MoverOP loop_mover = create_remodel_loop_mover( pose, sd, movable, loops );
	if ( !loop_mover ) {
		throw EXCN_Fold( type() + "::apply(): mover to remodel loops could not be constructed. Fail." );
	}

	// run remodel loop mover
	loop_mover->apply( pose );
	TR << "After remodel, FT=" << pose.fold_tree() << std::endl;

	// check status
	if ( loop_mover->get_last_move_status() != protocols::moves::MS_SUCCESS ) {
		throw EXCN_Fold( type() + "::apply(): RemodelLoopMover failed" );
	}

	// save modified SD
	remove_cutpoints( sd, loops );
	StructureDataFactory::get_instance()->save_into_pose( pose, sd );
}

void
RemodelLoopMoverPoseFolder::remove_cutpoints( StructureData & sd, protocols::loops::Loops const & loops ) const
{
	using protocols::loops::Loops;
	for ( Loops::const_iterator l=loops.begin(); l!=loops.end(); ++l ) {
		if ( l->cut() == 0 ) continue;
		SegmentName const seg_name = sd.segment_name( l->cut() );
		sd.set_cutpoint( seg_name, 0 );
		TR.Debug << "Removed cutpoint in StructureData from residue " << l->cut()
			<< " segment " << seg_name << std::endl;
	}
}

protocols::moves::MoverOP
RemodelLoopMoverPoseFolder::create_remodel_loop_mover(
	core::pose::Pose const & pose,
	StructureData const & sd,
	core::select::residue_selector::ResidueSubset const & movable,
	protocols::loops::Loops const & loops ) const
{
	// setup remodel for folding loops
	protocols::loops::LoopsOP loops_ptr( new protocols::loops::Loops( loops ) );
	protocols::forge::remodel::RemodelLoopMoverOP remodel =
		protocols::forge::remodel::RemodelLoopMoverOP( new protocols::forge::remodel::RemodelLoopMover( loops_ptr ) );
	remodel->set_keep_input_foldtree( true );
	remodel->set_seal_foldtree( false );
	remodel->false_movemap( *create_false_movemap( pose, movable ) );

	//if ( ( max_chainbreak_ > 0 ) && loops_to_close ) {
	//  remodel->max_linear_chainbreak( max_chainbreak_*loops_to_close );
	//}

	// select fragments
	components::Picker & picker = *components::Picker::get_instance();
	core::fragment::ConstantLengthFragSetCOP frag9 =
		picker.pick_and_cache_fragments( sd.ss(), sd.abego(), loops, pose.conformation().chain_endings(), 9 );
	core::fragment::ConstantLengthFragSetCOP frag3 =
		picker.pick_and_cache_fragments( sd.ss(), sd.abego(), loops, pose.conformation().chain_endings(), 3 );
	remodel->add_fragments( frag9 );
	remodel->add_fragments( frag3 );

	core::scoring::ScoreFunctionOP sfx = create_scorefxn( pose );
	remodel->scorefunction( *sfx );

	return remodel;
}

core::kinematics::MoveMapOP
RemodelLoopMoverPoseFolder::create_false_movemap(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueSubset const & movable ) const
{
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	for ( core::Size resid=1; resid<=movable.size(); ++resid ) {
		mm->set_bb( resid, movable[resid] );
		mm->set_chi( resid, movable[resid] );
	}

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		core::pose::symmetry::make_symmetric_movemap( pose, *mm );
	}

	TR.Debug << "Created false movemap: " << std::endl;
	mm->show( TR.Debug );
	TR.Debug << std::endl;
	return mm;
}

core::scoring::ScoreFunctionOP
RemodelLoopMoverPoseFolder::default_score_function()
{
	using core::scoring::ScoreFunctionOP;

	ScoreFunctionOP sfx = core::scoring::ScoreFunctionFactory::create_score_function( "fldsgn_cen.wts" );
	sfx->set_weight( core::scoring::vdw, 3.0 );
	sfx->set_weight( core::scoring::atom_pair_constraint, 1.0 );
	sfx->set_weight( core::scoring::coordinate_constraint, 1.0 );
	sfx->set_weight( core::scoring::angle_constraint, 1.0 );
	sfx->set_weight( core::scoring::dihedral_constraint, 1.0 );
	sfx->set_weight( core::scoring::hbond_sr_bb, 1.0 );
	sfx->set_weight( core::scoring::hbond_lr_bb, 2.0 );
	sfx->set_weight( core::scoring::rg, 0.1 );
	return sfx;
}

core::scoring::ScoreFunctionOP
RemodelLoopMoverPoseFolder::create_scorefxn( core::pose::Pose const & pose ) const
{
	core::scoring::ScoreFunctionOP sfx;

	if ( scorefxn_ ) {
		sfx = scorefxn_->clone();
	} else {
		sfx = default_score_function();
	}

	core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose( pose, sfx );

	return sfx;
}

void
RemodelLoopMoverPoseFolder::set_scorefxn( core::scoring::ScoreFunction const & sfxn )
{
	scorefxn_ = sfxn.clone();
}

} //protocols
} //denovo_design
} //components
