// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/domain_assembly/RedesignDomainAssemblyMover.cc
/// @brief

// Unit Headers
#include <devel/domain_assembly/RedesignDomainAssemblyMover.hh>

// Package Headers
#include <devel/domain_assembly/domain_assembly_setup.hh>
#include <devel/domain_assembly/domain_assembly.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/selection.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/util/kinematics_util.hh>
#include <core/scoring/constraints/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/pack/task/residue_selector/OrResidueSelector.hh>
#include <core/pack/task/residue_selector/ResidueIndexSelector.hh>
#include <core/pack/task/residue_selector/InterGroupInterfaceByVectorSelector.hh>
#include <core/pack/task/residue_selector/NotResidueSelector.hh>
#include <core/pack/task/residue_selector/ChainSelector.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// Utility Headers
#include <utility/io/izstream.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

// option key headers
#include <basic/options/keys/DomainAssembly.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ headers
#include <iostream>
#include <fstream>

namespace devel{
namespace domain_assembly{

RedesignDomainAssemblyMover::RedesignDomainAssemblyMover()
{
// basic initialization here
}

RedesignDomainAssemblyMover::~RedesignDomainAssemblyMover() {}
void
RedesignDomainAssemblyMover::apply( core::pose::Pose & pose )
{
	initialize();

	if( basic::options::option[ basic::options::OptionKeys::DomainAssembly::run_centroid ]() )
		run_centroid_stage( pose );

	if( get_last_move_status() != protocols::moves::MS_SUCCESS ) return;

	if( basic::options::option[ basic::options::OptionKeys::DomainAssembly::run_centroid_abinitio ]() )
		run_abinitio_centroid_stage( pose );

	if( basic::options::option[ basic::options::OptionKeys::DomainAssembly::run_fullatom ]() )
		run_fullatom_relax( pose );

	if( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() )
		evaluate_pose( pose );
	}

protocols::moves::MoverOP
RedesignDomainAssemblyMover::clone() const
{
	return new RedesignDomainAssemblyMover( *this );
}

protocols::moves::MoverOP
RedesignDomainAssemblyMover::fresh_instance() const
{
	return new RedesignDomainAssemblyMover;
}

std::string
RedesignDomainAssemblyMover::get_name() const
{
	return "RedesignDomainAssemblyMover";
}


void
RedesignDomainAssemblyMover::run_fullatom_stage( core::pose::Pose & pose )
{
// This is a part of Oana's original protocol refactored to use the ResidueSelector Framework and to
// redesign residues in the interface. However, the newly implemented version using FastRelax
// or side chain minimization while ramping up fa_rep

	// recover sidechains if pose has been loaded from centriod PDB...
	// rather check if pose is in centroid mode and try recovering sidechains.
	if( basic::options::option[ basic::options::OptionKeys::DomainAssembly::da_start_pdb ].user() )
		recover_sidechains( pose );
 	// otherwise warn about inserting random sidechains

	core::scoring::ScoreFunctionOP scorefxn( core::scoring::getScoreFunction() );

	std::vector< std::string > domain_definitions;
	get_domain_definition( pose, domain_definitions );

	using namespace core::pack::task;
	residue_selector::NotResidueSelectorOP not_rs = new core::pack::task::residue_selector::NotResidueSelector;
	residue_selector::OrResidueSelectorOP or_rs = new core::pack::task::residue_selector::OrResidueSelector;

	 // add all possible interdomain interfaces to the ORResidueSelector for repacking
	for( core::Size ii = 0; ii < domain_definitions.size(); ++ii ) {
		for( core::Size jj = ii+1; jj < domain_definitions.size(); ++jj ) {
			residue_selector::InterGroupInterfaceByVectorSelectorOP vector_rs = new residue_selector::InterGroupInterfaceByVectorSelector;
			vector_rs->group1_resstring( domain_definitions[ii] );
			vector_rs->group2_resstring( domain_definitions[jj] );
			or_rs->add_residue_selector( vector_rs );
		}
	}

	or_rs->add_residue_selector( new residue_selector::ResidueIndexSelector( get_linker_definition( pose ) ) );
	not_rs->set_residue_selector( or_rs );
	operation::OperateOnResidueSubsetOP repack_operation = new operation::OperateOnResidueSubset(
		new operation::PreventRepackingRLT(),
		not_rs	);

	// global repack of the side chains
	core::pack::task::PackerTaskOP base_packer_task( core::pack::task::TaskFactory::create_packer_task( pose ));
	base_packer_task -> initialize_from_command_line();
	base_packer_task -> restrict_to_repacking();
	base_packer_task -> set_bump_check( true );
	base_packer_task -> or_include_current( true );
	core::pack::pack_rotamers( pose, *scorefxn, base_packer_task );

	// Add the constraints, if any, to the score function
	core::scoring::constraints::add_fa_constraints_from_cmdline(pose, *scorefxn);

	protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *scorefxn, 0.8 /*temperature*/ ) );

	core::kinematics::MoveMapOP movemap_local = new core::kinematics::MoveMap( move_map() );
	// MOVER: small moves
	protocols::simple_moves::SmallMoverOP small_mover( new protocols::simple_moves::SmallMover( movemap_local , 0.8/*temp*/, 1 ) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 3.0 );
	small_mover->angle_max( 'L', 4.0 );

	// MOVER: rotamer trials
	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial_mover( new protocols::simple_moves::EnergyCutRotamerTrialsMover( scorefxn, *base_packer_task, mc, 0.01 /*energycut*/ ) );

	// MOVER minimization
	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( movemap_local, scorefxn, "dfpmin", 0.001, true /*use_nblist*/ );

	//// STAGE1 ////
	protocols::moves::SequenceMoverOP stage1_seq = new protocols::moves::SequenceMover;
	stage1_seq -> add_mover( small_mover );
	stage1_seq -> add_mover( pack_rottrial_mover );
	protocols::moves::TrialMoverOP stage1_trial = new protocols::moves::TrialMover( stage1_seq, mc );
	protocols::moves::RepeatMoverOP stage1_inner_loop = new protocols::moves::RepeatMover( stage1_trial, 25 /*100 cycles*/ );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= 20; ++i ) {
		stage1_inner_loop -> apply( pose );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
		core::pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
		// only repack at this point and redesign later?
		repack_operation->apply( pose, *this_packer_task );
		core::pack::pack_rotamers( pose, *scorefxn, this_packer_task );
	}
	std::cout << "end of stage one fullatom refinement" << std::endl;
	//// END STAGE1 ////

	//// STAGE2 ////
	protocols::moves::SequenceMoverOP stage2_seq = new protocols::moves::SequenceMover;
	stage2_seq->add_mover( small_mover );
	stage2_seq->add_mover( pack_rottrial_mover );
	stage2_seq->add_mover( min_mover );
	protocols::moves::TrialMoverOP stage2_trial = new protocols::moves::TrialMover( stage2_seq, mc );
	protocols::moves::RepeatMoverOP stage2_inner_loop = new protocols::moves::RepeatMover( stage2_trial, 50 /*cycles*/ );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= 20; ++i ) {
		stage2_inner_loop -> apply( pose );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
		core::pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
		// redesign here using a new taskop and selector.
		repack_operation->apply( pose, *this_packer_task );
		core::pack::pack_rotamers( pose, *scorefxn, this_packer_task );
	}

	std::cout << "end of stage two fullatom refinement" << std::endl;
	//// END STAGE2 ////
	mc -> recover_low( pose );

}

void RedesignDomainAssemblyMover::run_fullatom_relax( core::pose::Pose & pose ) {
using namespace core::pack::task;

	// recover sidechains if pose has been loaded from centriod PDB...
	// rather check if pose is in centroid mode and try recovering sidechains.
	if( basic::options::option[ basic::options::OptionKeys::DomainAssembly::da_start_pdb ].user() )
		recover_sidechains( pose );
 	// otherwise warn about inserting random sidechains

	core::scoring::ScoreFunctionOP scorefxn( core::scoring::getScoreFunction() );

	std::vector< std::string > domain_definitions;
	get_domain_definition( pose, domain_definitions );

	residue_selector::NotResidueSelectorOP not_interface_or_linker_rs = new core::pack::task::residue_selector::NotResidueSelector;
	residue_selector::OrResidueSelectorOP interface_or_linker_rs = new core::pack::task::residue_selector::OrResidueSelector;

	// add all possible interdomain interfaces to the ORResidueSelector for repacking
	for( core::Size ii = 0; ii < domain_definitions.size(); ++ii ) {
		for( core::Size jj = ii+1; jj < domain_definitions.size(); ++jj ) {
			residue_selector::InterGroupInterfaceByVectorSelectorOP vector_rs = new residue_selector::InterGroupInterfaceByVectorSelector;
			vector_rs->group1_resstring( domain_definitions[ii] );
			vector_rs->group2_resstring( domain_definitions[jj] );
			interface_or_linker_rs->add_residue_selector( vector_rs );
		}
	}

	interface_or_linker_rs->add_residue_selector( new residue_selector::ResidueIndexSelector( get_linker_definition( pose ) ) );
	not_interface_or_linker_rs->set_residue_selector( interface_or_linker_rs );
	operation::OperateOnResidueSubsetOP block_outside_interface_operation = new operation::OperateOnResidueSubset(
			new operation::PreventRepackingRLT(),
			not_interface_or_linker_rs	);

	residue_selector::ResidueIndexSelectorOP repack_only_rs = new residue_selector::ResidueIndexSelector( residues_to_repack_only_ );
	operation::OperateOnResidueSubsetOP block_design_operation = new operation::OperateOnResidueSubset(
		new operation::RestrictToRepackingRLT(),
		repack_only_rs
	);

	// global repack of the side chains
	core::pack::task::PackerTaskOP base_packer_task( core::pack::task::TaskFactory::create_packer_task( pose ));
	base_packer_task -> initialize_from_command_line();
	base_packer_task -> restrict_to_repacking();
	base_packer_task -> set_bump_check( true );
	base_packer_task -> or_include_current( true );
	core::pack::pack_rotamers( pose, *scorefxn, base_packer_task );

	// Add the constraints, if any, to the score function
	core::scoring::constraints::add_fa_constraints_from_cmdline(pose, *scorefxn);

	utility::vector1< bool > to_repack(false, pose.total_residue() );
	interface_or_linker_rs->apply( pose, to_repack );

	// allow minimazion of residues in the linker and the interdomain interface
	core::kinematics::MoveMapOP movemap_local = new core::kinematics::MoveMap( move_map() );
	movemap_local->set_chi( to_repack );
//	protocols::relax::FastRelaxOP frlx = new protocols::relax::FastRelax( scorefxn, 5);
//	frlx->set_movemap( movemap_local );
//	frlx->apply( pose );

 // Poor man's relax: Minimization, pack_rotamers, while ramping up fa_rep energy term
	core::Size outer_iterations = 4;
	core::Size inner_iterations = 10;
	core::Real final_fa_rep = 0.44;

	for( Size i_outer = 1; i_outer <= outer_iterations; ++i_outer ) {
		for( Size i_inner = 1; i_inner <= inner_iterations; ++i_inner ) {
			core::Real fa_rep_weight = ( 0.1 + 0.9/(inner_iterations-1) * (i_inner-1) ) * final_fa_rep;
			scorefxn->set_weight( core::scoring::fa_rep , fa_rep_weight );
			interface_or_linker_rs->apply( pose, to_repack );
			movemap_local->set_chi( to_repack );
			movemap_local->set_bb( to_repack );

			protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( movemap_local, scorefxn, "dfpmin_armijo_nonmonotone", 0.001, true /*use_nblist*/ );
			min_mover->apply( pose );

			// redesign of residues on the LOV-side of the interface
			core::pack::task::PackerTaskOP design_packer_task( core::pack::task::TaskFactory::create_packer_task( pose ));
			design_packer_task -> initialize_from_command_line();
			design_packer_task -> set_bump_check( true );
			design_packer_task -> or_include_current( true );
			block_outside_interface_operation -> apply( pose, *design_packer_task );
			block_design_operation -> apply( pose, *design_packer_task );
			core::pack::pack_rotamers( pose, *scorefxn, design_packer_task );
		}
	}

	// rewrite this or put this in a function at some point. Just want to know if it works now.
	base_packer_task = core::pack::task::TaskFactory::create_packer_task( pose );
	base_packer_task -> initialize_from_command_line();
	base_packer_task -> restrict_to_repacking();
	base_packer_task -> set_bump_check( true );
	base_packer_task -> or_include_current( true );
	core::pack::pack_rotamers( pose, *scorefxn, base_packer_task );
}

void
RedesignDomainAssemblyMover::initialize()
{
	if ( ! movemap_set() )
		initialize_movemap_from_commandline();

	if ( ! fragsets_set() )
		initialize_fragments_from_commandline();

	if ( basic::options::option [ basic::options::OptionKeys::DomainAssembly::da_eval_pose_map ].user() )
		initialize_pose_map_from_commandline();

	if ( target_pose().empty() && basic::options::option[ basic::options::OptionKeys::in::file::native ].user() )
		initialize_target_pose();

	if ( starting_pose().empty() && basic::options::option[ basic::options::OptionKeys::DomainAssembly::da_start_pdb ].user() )
		initialize_start_pose_from_commandline();

	if( buried_residues().empty() && basic::options::option[ basic::options::OptionKeys::DomainAssembly::da_require_buried ].user() )
		initialize_buried_from_commandline();

	if( residues_to_repack_only_.empty() && basic::options::option[ basic::options::OptionKeys::DomainAssembly::residues_repack_only ].user() )
		residues_to_repack_only_ = basic::options::option[ basic::options::OptionKeys::DomainAssembly::residues_repack_only ]();
}

} // namespace domain_assembly
} // namespace devel

