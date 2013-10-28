// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/domain_assembly/DomainAssemblyMover.cc
/// @brief

// Unit Headers
#include <devel/domain_assembly/DomainAssemblyMover.hh>

// Package Headers
#include <devel/domain_assembly/domain_assembly_setup.hh>
#include <devel/domain_assembly/domain_assembly.hh>

// Project Headers
#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ResidueConnection.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/fragment/ConstantLengthFragSet.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/util/SwitchResidueTypeSet.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/FragmentMover.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

// Utility Headers

// Basic headers
#include <basic/options/option.hh>

// option key headers
#include <basic/options/keys/DomainAssembly.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>


namespace devel{
namespace domain_assembly{

DomainAssemblyMover::DomainAssemblyMover() :
	movemap_set_( false )
{
	// basic initialization here
}

void
DomainAssemblyMover::apply( core::pose::Pose & pose )
{
	// barebones protocol
	// 0. Intialize all data that hasn't been initialized programmatically
	// 1. run a low-resolution protocol w/ fragment insertion (if available)
	// 2. run a high-resolution refinement w/ small moves & repacking & minimization
	initialize();
	run_centroid_stage( pose );
	run_fullatom_stage( pose );
}

protocols::moves::MoverOP
DomainAssemblyMover::clone() const
{
	return new DomainAssemblyMover( *this );
}

protocols::moves::MoverOP
DomainAssemblyMover::fresh_instance() const
{
	return new DomainAssemblyMover;
}

std::string
DomainAssemblyMover::get_name() const
{
	return "DomainAssemblyMover";
}

void
DomainAssemblyMover::run_centroid_stage( core::pose::Pose & pose )
{
	Size inside_steps_stage1 ( 50 );
	Size outside_steps_stage1 ( 15 );
	Size inside_steps_stage2 ( 50 );
	Size outside_steps_stage2 ( 5 );

	core::pose::Pose const saved_input_pose( pose ); //used to return sidechains later
	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );
	core::scoring::ScoreFunctionOP scorefxn_centroid( core::scoring::ScoreFunctionFactory::create_score_function( core::scoring::CENTROID_WTS ) );
	protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *scorefxn_centroid, 0.8 /*temperature*/ ) );

	// more traditional small moves
	protocols::simple_moves::SmallMoverOP small_mover( new protocols::simple_moves::SmallMover( movemap_, 0.8/*temp*/, 1/*nmoves*/ ) );
	small_mover->angle_max( 'H', 2.0 );  // max angle displacement 180 degrees
	small_mover->angle_max( 'E', 4.0 );
	small_mover->angle_max( 'L', 4.0 );

	//// STAGE 1 /////
	protocols::moves::TrialMoverOP centroid_trial_mover = NULL;
	if ( fragset3mer_ ) {
		protocols::simple_moves::FragmentMoverOP frag_mover = new protocols::simple_moves::ClassicFragmentMover(fragset3mer_, movemap_ );
		centroid_trial_mover = new protocols::moves::TrialMover( frag_mover, mc );
	} else {
		// if no fragments, use coarse small moves
		Size nmoves ( 1 );
		protocols::simple_moves::SmallMoverOP coarse_small_mover( new protocols::simple_moves::SmallMover( movemap_, 0.8/*temp*/, nmoves ) );
		coarse_small_mover->angle_max( 'H', 180.0 );  // max angle displacement 180 degrees
		coarse_small_mover->angle_max( 'E', 180.0 );
		coarse_small_mover->angle_max( 'L', 180.0 );
		centroid_trial_mover = new protocols::moves::TrialMover( coarse_small_mover, mc );
	}

	protocols::moves::RepeatMoverOP inner_centroid_loop = new protocols::moves::RepeatMover( centroid_trial_mover, inside_steps_stage1 );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= outside_steps_stage1; ++i ) {
		inner_centroid_loop -> apply( pose );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
	}
	std::cout << "end of stage one centroid refinement" << std::endl;
	//// END STAGE 1 ////

	////// STAGE 2 ///////
	protocols::moves::TrialMoverOP stage2_trial = new protocols::moves::TrialMover( small_mover, mc );
	protocols::moves::RepeatMoverOP stage2 = new protocols::moves::RepeatMover( stage2_trial, inside_steps_stage2 );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= outside_steps_stage2; ++i ) {
		stage2 -> apply( pose );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
	}
	std::cout << "end of stage two centroid refinement" << std::endl;
	///////  END_STAGE 2 ////////

	mc -> recover_low( pose );
	protocols::simple_moves::ReturnSidechainMover return_sidechains( saved_input_pose );
	return_sidechains.apply( pose );

}

//protocols::toolbox::task_operations::RestrictInterGroupVectorOperationOP
//DomainAssemblyMover::prepare_domain_interface_restiction_operation(
//	core::pose::Pose const & pose
//)
//{
//	// both flexible and rigid residue ranges
//	std::list< std::pair< core::Size, core::Size > > region_ranges;
//	core::Size region_start( 0 ), region_end( 0 );
//	bool last_flexible( false );
//	for ( core::Size ii = 1; ii <= pose.total_residue(); ++ii ) {
//		if ( movemap_.get_bb( ii ) ) {
//			if ( ii == 1 ) {
//				region_start = ii;
//				last_flexible = true;
//			} else if ( ! last_flexible ) {
//				// aha! end of range of an inflexible region
//				region_ranges.push_back( std::make_pair( region_start, region_end ) );
//				last_flexible = true;
//				region_start = ii;
//			} else {
//				region_end = ii;
//			}
//		} else {
//			if ( ii == 1 ) {
//				region_start == ii;
//				last_flexible = false;
//			} else if ( last_flexible ) {
//				// aha! end of a range of a flexible region
//				region_ranges.push_back( std::make_pair( region_start, region_end ) );
//				last_flexible = false;
//				region_start = ii;
//			} else {
//				region_end = ii;
//			}
//		}
//	}
//
//	utility::vector1< std::set< core::Size > > region_sets;
//	region_sets.reserve( region_ranges.size() );
//	for ( std::list< std::pair< core::Size, core::Size > >::const_iterator
//			iter = region_ranges.begin(),
//			iter_end = region_ranges.end();
//			iter != iter_end; ++iter ) {
//		std::set< core::Size > region_set;
//		for ( core::Size ii = iter->first; ii <= iter->second; ++ii ) region_set.insert( ii );
//		region_sets.push_back( region_set );
//	}
//
//	protocols::toolbox::task_operations::RestrictInterGroupVectorOperationOP operation = new
//		protocols::toolbox::task_operations::RestrictInterGroupVectorOperation;
//
//	// add all pairs of region ranges
//	for ( core::Size ii = 1; ii <= region_sets.size(); ++ii ) {
//		for ( core::Size jj = ii+1; jj <= region_sets.size(); ++jj ) {
//			op->insert_pair( std::makepair( region_sets[ ii ], region_sets[ jj ] ) );
//		}
//	}
//
//	return op;
//
//}


void
DomainAssemblyMover::run_fullatom_stage( core::pose::Pose & pose )
{
	//core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD );
	core::scoring::ScoreFunctionOP scorefxn( core::scoring::getScoreFunction() );




	// for each residue - identify the nearest movable residue forward and backwards in sequence
	// this will be used to help determine which side chains to repack
	utility::vector1 < std::pair < Size, Size > > nearest_movable_residues;
	find_nearest_movable_residues( movemap_, pose, nearest_movable_residues );

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

	// MOVER: small moves
  protocols::simple_moves::SmallMoverOP small_mover( new protocols::simple_moves::SmallMover( movemap_, 0.8/*temp*/, 1 ) );
  small_mover->angle_max( 'H', 2.0 );
  small_mover->angle_max( 'E', 3.0 );
  small_mover->angle_max( 'L', 4.0 );

	// MOVER: rotamer trials
	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial_mover( new protocols::simple_moves::EnergyCutRotamerTrialsMover( scorefxn, *base_packer_task, mc, 0.01 /*energycut*/ ) );

	// MOVER minimization
	protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover( movemap_, scorefxn, "dfpmin", 0.001, true /*use_nblist*/ );
	//TrialMoverOP min_trial_mover( new TrialMover( min_mover, mc ) );

	// Initial Minimization //
	//std::cout << " Score before initial fullatom minimization " << mc->last_accepted_score() << std::endl;
	//min_mover -> apply (pose);
	//std::cout << " Score after initial fullatom minimization " << mc->last_accepted_score() << std::endl;

	//// STAGE1 ////
	protocols::moves::SequenceMoverOP stage1_seq = new protocols::moves::SequenceMover;
	stage1_seq -> add_mover( small_mover );
	stage1_seq -> add_mover( pack_rottrial_mover );
	protocols::moves::TrialMoverOP stage1_trial = new protocols::moves::TrialMover( stage1_seq, mc );
	protocols::moves::RepeatMoverOP stage1_inner_loop = new protocols::moves::RepeatMover( stage1_trial, 15 /*100 cycles*/ );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= 10; ++i ) {
		stage1_inner_loop -> apply( pose );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
		utility::vector1<bool> repack_residues( pose.total_residue(), false );
		da_residues_to_repack( movemap_, nearest_movable_residues, pose, repack_residues );
		core::pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
		this_packer_task -> restrict_to_residues( repack_residues );
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
	protocols::moves::RepeatMoverOP stage2_inner_loop = new protocols::moves::RepeatMover( stage2_trial, 10 /*cycles*/ );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= 10; ++i ) {
		stage2_inner_loop -> apply( pose );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
		utility::vector1<bool> repack_residues( pose.total_residue(), false );
		da_residues_to_repack( movemap_, nearest_movable_residues, pose, repack_residues );
		core::pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
		this_packer_task -> restrict_to_residues( repack_residues );
		core::pack::pack_rotamers( pose, *scorefxn, this_packer_task );
	}
	std::cout << "end of stage two fullatom refinement" << std::endl;
	//// END STAGE2 ////
	mc -> recover_low( pose );

}



void
DomainAssemblyMover::initialize()
{
	if ( ! movemap_set_ ) initialize_movemap_from_commandline();
	if ( ! fragset3mer_ ) initialize_fragments_from_commandline();
}

void
DomainAssemblyMover::initialize_movemap_from_commandline()
{
	movemap_ = read_movemap_from_da_linker_file();
	movemap_set_ = true;
}

void
DomainAssemblyMover::initialize_fragments_from_commandline()
{
	// read fragments file
	core::fragment::ConstantLengthFragSetOP fragset3mer = NULL;
	if (basic::options::option[ basic::options::OptionKeys::in::file::frag3].user()){
		fragset3mer = new core::fragment::ConstantLengthFragSet( 3 );
		fragset3mer->read_fragment_file( basic::options::option[ basic::options::OptionKeys::in::file::frag3 ]() );
	}
}


}
}
