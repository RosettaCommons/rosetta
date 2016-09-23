// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/selection.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/util/kinematics_util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/fragment/ConstantLengthFragSet.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/rms_util.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/OperateOnResidueSubset.hh>
#include <core/pack/task/operation/ResLvlTaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <core/select/residue_selector/OrResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/InterGroupInterfaceByVectorSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/ChainSelector.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/ChainSplitMover.hh>
#include <protocols/simple_filters/DomainInterfaceFilter.hh>

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

namespace devel {
namespace domain_assembly {

DomainAssemblyMover::DomainAssemblyMover() :
	movemap_set_( false )
{
	// basic initialization here
}

DomainAssemblyMover::~DomainAssemblyMover() = default;

void
DomainAssemblyMover::apply( core::pose::Pose & pose )
{
	initialize();

	if ( basic::options::option[ basic::options::OptionKeys::DomainAssembly::run_centroid ]() ) {
		run_centroid_stage( pose );
	}

	if ( get_last_move_status() != protocols::moves::MS_SUCCESS ) return;

	if ( basic::options::option[ basic::options::OptionKeys::DomainAssembly::run_centroid_abinitio ]() ) {
		run_abinitio_centroid_stage( pose );
	}

	if ( basic::options::option[ basic::options::OptionKeys::DomainAssembly::run_fullatom ]() ) {
		run_fullatom_relax( pose );
	}

	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		evaluate_pose( pose );
	}
}

void
DomainAssemblyMover::evaluate_pose( core::pose::Pose const & pose) const {
	basic::Tracer TR( "devel.domain_assembly.DomainAssemblyMover" );

	core::Real sum_sq = target_rmsd( pose );

	TR << "input, job, env, pair, cbeta, vdw, rg, cenpack, hs_pair, sheet, total, rms" << std::endl;
	TR << pose.pdb_info()->name() << ", "  << protocols::jd2::get_current_job()->nstruct_index()
		<< ", " << pose.energies().total_energies()[core::scoring::env]
		<< ", " << pose.energies().total_energies()[core::scoring::pair]
		<< ", " << pose.energies().total_energies()[core::scoring::cbeta]
		<< ", " << pose.energies().total_energies()[core::scoring::vdw]
		<< ", " << pose.energies().total_energies()[core::scoring::rg]
		<< ", " << pose.energies().total_energies()[core::scoring::cenpack]
		<< ", " << pose.energies().total_energies()[core::scoring::hs_pair]
		// << ", " << pose.energies().total_energies()[core::scoring::ss_pair]
		// << ", " << pose.energies().total_energies()[core::scoring::rsigma]
		<< ", " << pose.energies().total_energies()[core::scoring::sheet]
		<< ", " <<  pose.energies().total_energies()[core::scoring::total_score]
		<< ", " << sum_sq << std::endl;

	// interface analyzer
	// hard-coded split point!
	// we should either use the middle of the flexible region
	// OR even better, the actual domain split point. This we only get if we
	// run domain assembly as one step protocol without separate calls to
	// setup and optimization phase.
	// wanted behavior:
	// look for pose split points in datacache, else take linker midpoints
	protocols::simple_moves::ChainSplitMover split_mover( 137 );

	core::pose::Pose pose_to_eval( pose );
	split_mover.apply( pose_to_eval );

	protocols::analysis::InterfaceAnalyzerMover interface_mover;
	interface_mover.apply( pose_to_eval );

	interface_mover.set_use_tracer(true);
	interface_mover.report_data();
}

protocols::moves::MoverOP
DomainAssemblyMover::clone() const
{
	return protocols::moves::MoverOP( new DomainAssemblyMover( *this ) );
}

protocols::moves::MoverOP
DomainAssemblyMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new DomainAssemblyMover );
}

std::string
DomainAssemblyMover::get_name() const
{
	return "DomainAssemblyMover";
}

void
DomainAssemblyMover::run_centroid_stage( core::pose::Pose & pose )
{
	Size inside_steps_stage1 ( 100 );
	Size outside_steps_stage1 ( 30 );
	Size inside_steps_stage2 ( 100 );
	Size outside_steps_stage2 ( 20 );

	basic::Tracer TR( "DomainAssemblyMover::run_centroid_stage" );

	core::pose::Pose const saved_input_pose( pose ); //used to return sidechains later
	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID_t );
	core::scoring::ScoreFunctionOP scorefxn_centroid( core::scoring::ScoreFunctionFactory::create_score_function( "score3" ) );
	scorefxn_centroid->set_weight( core::scoring::rama, 0.2 );
	protocols::moves::MonteCarloOP mc( new protocols::moves::MonteCarlo( pose, *scorefxn_centroid, 0.8 /*temperature*/ ) );

	// more traditional small moves
	protocols::simple_moves::SmallMoverOP small_mover( new protocols::simple_moves::SmallMover( movemap_, 0.8/*temp*/, 1/*nmoves*/ ) );
	small_mover->angle_max( 'H', 2.0 );  // max angle displacement 180 degrees
	small_mover->angle_max( 'E', 4.0 );
	small_mover->angle_max( 'L', 4.0 );

	//// STAGE 1 /////
	protocols::moves::TrialMoverOP centroid_trial_mover = nullptr;
	if ( fragset3mer_ ) {
		protocols::simple_moves::FragmentMoverOP frag_mover( new protocols::simple_moves::ClassicFragmentMover(fragset3mer_, movemap_ ) );
		centroid_trial_mover = protocols::moves::TrialMoverOP( new protocols::moves::TrialMover( frag_mover, mc ) );
	} else {
		// if no fragments, use coarse small moves
		Size nmoves ( 1 );
		protocols::simple_moves::SmallMoverOP coarse_small_mover( new protocols::simple_moves::SmallMover( movemap_, 0.8/*temp*/, nmoves ) );
		coarse_small_mover->angle_max( 'H', 180.0 );  // max angle displacement 180 degrees
		coarse_small_mover->angle_max( 'E', 180.0 );
		coarse_small_mover->angle_max( 'L', 180.0 );
		centroid_trial_mover = protocols::moves::TrialMoverOP( new protocols::moves::TrialMover( coarse_small_mover, mc ) );
	}

	protocols::moves::RepeatMoverOP inner_centroid_loop( new protocols::moves::RepeatMover( centroid_trial_mover, inside_steps_stage1 ) );

	for ( Size i = 1; i <= outside_steps_stage1; ++i ) {
		inner_centroid_loop -> apply( pose );
	}
	//// END STAGE 1 ////

	////// STAGE 2 ///////
	protocols::moves::TrialMoverOP stage2_trial( new protocols::moves::TrialMover( small_mover, mc ) );
	protocols::moves::RepeatMoverOP stage2( new protocols::moves::RepeatMover( stage2_trial, inside_steps_stage2 ) );

	for ( Size i = 1; i <= outside_steps_stage2; ++i ) {
		stage2 -> apply( pose );
	}
	///////  END_STAGE 2 ////////

	mc -> recover_low( pose );

	//////// STAGE 3: FILTERING ////////
	if ( !buried_residues_.empty() ) {
		std::vector< std::string > domain_definitions;
		get_domain_definition( pose, domain_definitions );
		std::set< core::Size > want_buried( core::pose::get_resnum_list( buried_residues_, pose ) );

		//std::vector< std::pair< protocols::filters::FilterOP, protocols::filters::boolean_operations > > compound_filter_input;

		core::Size num_true = 0;
		// step 1: go through all domains
		// step 2: find all residues requested to be buried
		// step 3: add a filter with that potential interface
		// step 4: or-combine the filters, at least one of the domains should bury any
		for ( core::Size buried_it : want_buried ) {
			std::set< core::Size > residues_in_other_domains;
			for ( std::vector< std::string >::const_iterator domain_it = domain_definitions.begin();
					domain_it != domain_definitions.end();
					++domain_it ) {
				std::set< core::Size > domain_set( core::pose::get_resnum_list( *domain_it, pose ) );
				if ( domain_set.find( buried_it ) == domain_set.end() ) { // residue not contained in domain! Add!
					residues_in_other_domains.insert( domain_set.begin(), domain_set.end() );
				}
			}
			if ( !residues_in_other_domains.empty() ) { // this can only happen if we're looking at exactly one domain...
				protocols::simple_filters::DomainInterfaceFilterOP di_filter( new protocols::simple_filters::DomainInterfaceFilter );
				std::set< core::Size > one_residue_dummy_set;
				one_residue_dummy_set.insert( buried_it );
				di_filter->query_region( one_residue_dummy_set );
				di_filter->target_region( residues_in_other_domains );
				if ( di_filter->apply( pose ) ) num_true++;
				//    compound_filter_input.push_back( std::make_pair( di_filter, protocols::filters::AND ) );
			}
		}

		//  protocols::filters::CompoundFilter all_and_filter( compound_filter_input );

		//  if( !all_and_filter.apply( pose ) ) {
		//   set_last_move_status( protocols::moves::FAIL_RETRY );
		//  }
		if ( (core::Real)num_true / want_buried.size() < 0.5 ) {
			set_last_move_status( protocols::moves::FAIL_RETRY );
		}

	}
	//////// END STAGE 3 ////////

}
void DomainAssemblyMover::run_abinitio_centroid_stage( core::pose::Pose & pose ) {

	if ( !fragsets_set() ) {
		throw utility::excn::EXCN_Msg_Exception( "No FragmentSet specified - exiting!\n" );
	}

	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID_t );
	core::scoring::ScoreFunctionOP scorefxn_wts3 = core::scoring::ScoreFunctionFactory::create_score_function( "score3" );
	protocols::abinitio::AbrelaxApplication abrelax_app;
	do {
		protocols::abinitio::ClassicAbinitio abinit( fragset9mer_, fragset3mer_, movemap_ );
		abinit.init( pose );
		abinit.apply( pose );
	} while( !abrelax_app.check_filters( pose ) );

	// re-score although we should be up to date.
	scorefxn_wts3->score(pose);
	//recover sidechains from starting structures
	// to_fullatom.apply( pose );
	// recover_sidechains.apply( pose );

}

void
DomainAssemblyMover::run_fullatom_stage( core::pose::Pose & pose )
{
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using core::select::residue_selector::ResidueSelectorOP;
	using core::select::residue_selector::ResidueSelectorCOP;

	// recover sidechains if pose has been loaded from centriod PDB...
	// FIX THIS to get wanted behavior
	if ( basic::options::option[ basic::options::OptionKeys::DomainAssembly::da_start_pdb ].user() ) {
		recover_sidechains( pose );
	}

	core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );

	std::vector< std::string > domain_definitions;
	get_domain_definition( pose, domain_definitions );

	core::select::residue_selector::NotResidueSelectorOP not_rs( new core::select::residue_selector::NotResidueSelector );
	core::select::residue_selector::OrResidueSelectorOP or_rs( new core::select::residue_selector::OrResidueSelector );

	// add all possible interdomain interfaces to the ORResidueSelector for repacking
	for ( core::Size ii = 0; ii < domain_definitions.size(); ++ii ) {
		for ( core::Size jj = ii+1; jj < domain_definitions.size(); ++jj ) {
			core::select::residue_selector::InterGroupInterfaceByVectorSelectorOP vector_rs( new core::select::residue_selector::InterGroupInterfaceByVectorSelector );
			vector_rs->group1_resstring( domain_definitions[ii] );
			vector_rs->group2_resstring( domain_definitions[jj] );
			or_rs->add_residue_selector( vector_rs );
		}
	}

	or_rs->add_residue_selector( ResidueSelectorCOP( ResidueSelectorOP( new core::select::residue_selector::ResidueIndexSelector( get_linker_definition( pose ) ) ) ) );
	not_rs->set_residue_selector( or_rs );
	operation::OperateOnResidueSubsetOP rsOperation( new operation::OperateOnResidueSubset( ResLvlTaskOperationCOP( new operation::PreventRepackingRLT() ), not_rs ) );

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
	protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( movemap_, scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true /*use_nblist*/ ) );

	//// STAGE1 ////
	protocols::moves::SequenceMoverOP stage1_seq( new protocols::moves::SequenceMover );
	stage1_seq -> add_mover( small_mover );
	stage1_seq -> add_mover( pack_rottrial_mover );
	protocols::moves::TrialMoverOP stage1_trial( new protocols::moves::TrialMover( stage1_seq, mc ) );
	protocols::moves::RepeatMoverOP stage1_inner_loop( new protocols::moves::RepeatMover( stage1_trial, 25 /*100 cycles*/ ) );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= 20; ++i ) {
		stage1_inner_loop -> apply( pose );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
		//utility::vector1<bool> repack_residues( pose.size(), false );
		//da_residues_to_repack( movemap_, nearest_movable_residues, pose, repack_residues );
		core::pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
		rsOperation->apply( pose, *this_packer_task );
		core::pack::pack_rotamers( pose, *scorefxn, this_packer_task );
	}
	std::cout << "end of stage one fullatom refinement" << std::endl;
	//// END STAGE1 ////

	//// STAGE2 ////
	protocols::moves::SequenceMoverOP stage2_seq( new protocols::moves::SequenceMover );
	stage2_seq->add_mover( small_mover );
	stage2_seq->add_mover( pack_rottrial_mover );
	stage2_seq->add_mover( min_mover );
	protocols::moves::TrialMoverOP stage2_trial( new protocols::moves::TrialMover( stage2_seq, mc ) );
	protocols::moves::RepeatMoverOP stage2_inner_loop( new protocols::moves::RepeatMover( stage2_trial, 50 /*cycles*/ ) );

	std::cout << "   Current  Low " << std::endl;
	for ( Size i = 1; i <= 20; ++i ) {
		stage2_inner_loop -> apply( pose );
		std::cout << i << "  " << mc->last_accepted_score() << "  " << mc->lowest_score() << std::endl;
		// utility::vector1<bool> repack_residues( pose.size(), false );
		// da_residues_to_repack( movemap_, nearest_movable_residues, pose, repack_residues );
		core::pack::task::PackerTaskOP this_packer_task( base_packer_task->clone() );
		rsOperation->apply( pose, *this_packer_task );
		core::pack::pack_rotamers( pose, *scorefxn, this_packer_task );
	}
	std::cout << "end of stage two fullatom refinement" << std::endl;
	//// END STAGE2 ////
	mc -> recover_low( pose );

}

void DomainAssemblyMover::run_fullatom_relax( core::pose::Pose & pose ) {
	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	using core::select::residue_selector::ResidueSelectorOP;
	using core::select::residue_selector::ResidueSelectorCOP;

	// recover sidechains if pose has been loaded from centriod PDB...
	// FIX THIS to get wanted behavior
	if ( basic::options::option[ basic::options::OptionKeys::DomainAssembly::da_start_pdb ].user() ) {
		recover_sidechains( pose );
	}

	core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );

	std::vector< std::string > domain_definitions;
	get_domain_definition( pose, domain_definitions );

	core::select::residue_selector::NotResidueSelectorOP not_rs( new core::select::residue_selector::NotResidueSelector );
	core::select::residue_selector::OrResidueSelectorOP or_rs( new core::select::residue_selector::OrResidueSelector );

	// add all possible interdomain interfaces to the ORResidueSelector for repacking
	for ( core::Size ii = 0; ii < domain_definitions.size(); ++ii ) {
		for ( core::Size jj = ii+1; jj < domain_definitions.size(); ++jj ) {
			core::select::residue_selector::InterGroupInterfaceByVectorSelectorOP vector_rs( new core::select::residue_selector::InterGroupInterfaceByVectorSelector );
			vector_rs->group1_resstring( domain_definitions[ii] );
			vector_rs->group2_resstring( domain_definitions[jj] );
			or_rs->add_residue_selector( vector_rs );
		}
	}

	or_rs->add_residue_selector( ResidueSelectorCOP( ResidueSelectorOP( new core::select::residue_selector::ResidueIndexSelector( get_linker_definition( pose ) ) ) ) );
	not_rs->set_residue_selector( or_rs );

	operation::OperateOnResidueSubsetOP repack_interface_operation( new operation::OperateOnResidueSubset( operation::ResLvlTaskOperationOP( new operation::PreventRepackingRLT() ), not_rs ) );

	// global repack of the side chains
	core::pack::task::PackerTaskOP base_packer_task( core::pack::task::TaskFactory::create_packer_task( pose ));
	base_packer_task -> initialize_from_command_line();
	base_packer_task -> restrict_to_repacking();
	base_packer_task -> set_bump_check( true );
	base_packer_task -> or_include_current( true );
	core::pack::pack_rotamers( pose, *scorefxn, base_packer_task );

	// Add the constraints, if any, to the score function
	core::scoring::constraints::add_fa_constraints_from_cmdline(pose, *scorefxn);

	utility::vector1< bool > to_repack = or_rs->apply( pose );


	// allow minimazion of residues in the linker and the interdomain interface
	movemap_->set_chi( to_repack );
	// movemap_->set_bb ( to_repack ); // seems like a fun idea, didnt give good results though.

	// almost in all cases very similar to fa_rep ramp x {min,repack}
	// protocols::relax::FastRelaxOP frlx = new protocols::relax::FastRelax( scorefxn, 5);
	// frlx->set_movemap( movemap_ );
	// frlx->apply( pose );

	// Poor man's relax: Minimization, pack_rotamers, while ramping up fa_rep energy term
	core::Size outer_iterations = 5;
	core::Size inner_iterations = 10;
	core::Real final_fa_rep = 0.44;

	for ( Size i_outer = 1; i_outer <= outer_iterations; ++i_outer ) {
		for ( Size i_inner = 1; i_inner <= inner_iterations; ++i_inner ) {
			core::Real fa_rep_weight = (0.1 + 0.9/(inner_iterations-1) * (i_inner-1)) * final_fa_rep;
			scorefxn->set_weight( core::scoring::fa_rep , fa_rep_weight );
			to_repack = or_rs->apply( pose );
			movemap_->set_chi( to_repack );
			//movemap_->set_bb( to_repack );

			protocols::simple_moves::MinMoverOP min_mover( new protocols::simple_moves::MinMover( movemap_, scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true /*use_nblist*/ ) );
			min_mover->apply( pose );

			// repack residues at the domain interface
			core::pack::task::PackerTaskOP interface_repack_task( core::pack::task::TaskFactory::create_packer_task( pose ));
			interface_repack_task -> initialize_from_command_line();
			interface_repack_task -> set_bump_check( true );
			interface_repack_task -> or_include_current( true );
			interface_repack_task -> restrict_to_repacking();

			repack_interface_operation -> apply( pose, *interface_repack_task );
			core::pack::pack_rotamers( pose, *scorefxn, interface_repack_task );
		}
	}

	// rewrite this or put this in a function at some point. Just want to now if it works now.
	// end with global repack
	base_packer_task =  core::pack::task::TaskFactory::create_packer_task( pose );
	base_packer_task -> initialize_from_command_line();
	base_packer_task -> restrict_to_repacking();
	base_packer_task -> set_bump_check( true );
	base_packer_task -> or_include_current( true );
	core::pack::pack_rotamers( pose, *scorefxn, base_packer_task );
}

// both the domain and the linker definition are fixed at initialization.
// recomputing them for each run of the fullatom protocol seems a little wasteful
// but this is nothing compared to the actual computational cost of the protocol.
void DomainAssemblyMover::get_domain_definition( core::pose::Pose const & pose, std::vector< std::string > & domains ) const {
	core::Size last_domain_start = 1;
	domains.clear();
	bool in_linker = movemap_->get_bb( 1 );
	for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( movemap_->get_bb( ii ) && !in_linker ) { // first residue of a linker
			in_linker = true;
			std::stringstream current_domain;
			current_domain << last_domain_start << "-" << ii-1;
			domains.push_back( current_domain.str() );
		} else if ( !movemap_->get_bb( ii ) && in_linker ) { // first residue after a linker
			in_linker = false;
			last_domain_start = ii-1;
		}
	}
	if ( !in_linker ) {
		std::stringstream current_domain;
		current_domain << last_domain_start << "-" << pose.size();
		domains.push_back( current_domain.str() );
	}
}

std::string DomainAssemblyMover::get_linker_definition( core::pose::Pose const & pose ) const {
	std::stringstream linker_str;
	core::Size last_linker_start = 1;
	bool in_linker = movemap_->get_bb( 1 );
	for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( movemap_->get_bb( ii ) && !in_linker ) {
			in_linker = true;
			last_linker_start = ii;
		} else if ( !movemap_->get_bb( ii ) && in_linker ) {
			in_linker = false;
			linker_str << last_linker_start << "-" << ii-1 << ",";
		}
	}
	if ( in_linker ) {
		linker_str << last_linker_start << "-" << pose.size();
	}

	return linker_str.str();
}

core::Real DomainAssemblyMover::target_rmsd( core::pose::Pose const & pose ) const {
	if ( target_pose_map_.empty() ) {
		// define corresponding residues to be compared
		// this is only useful for LOV-RAC and I should get rid of this asap.
		std::map< core::Size, core::Size > correspondence;
		int offset = pose.size() - target_pose_.size() - 1; // !!!!! -1 is a 2wkp vs 1mh1 specific hack because 2wkp is shorter by 1
		for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
			if ( (int)ii > offset && (ii - offset) <= target_pose_.size() ) {
				if ( pose.residue( ii ).name3() == target_pose_.residue( ii - offset ).name3() ) {
					correspondence.insert( std::make_pair( ii, ii-offset ) );
				}
			}
		}
		return core::scoring::CA_rmsd( pose, target_pose_, correspondence );
	} else {
		return core::scoring::CA_rmsd( pose, target_pose_, target_pose_map_);
	}
}

// calculate RMSD without previous alignment.
core::Real DomainAssemblyMover::target_rmsd_no_align( core::pose::Pose const & pose ) const {
	core::Size num_match = 0;
	core::Real sum_sq = 0;
	if ( target_pose_map_.empty() ) {
		// define corresponding residues to be compared
		// get rid of this asap
		int offset = pose.size() - target_pose_.size()  - 1; // !!!!! -1 is a 2wkp vs 1mh1 specific hack because 2wkp is shorter by 1
		for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
			if ( (int)ii > offset && (ii - offset) <= target_pose_.size() ) {
				if ( pose.residue( ii ).name3() == target_pose_.residue( ii - offset ).name3() ) {
					sum_sq += pose.residue( ii ).xyz( "CA" ).distance_squared( target_pose_.residue( ii-offset  ).xyz( "CA" ) );
					num_match++;
				}
			}
		}
	} else {
		core::Size num_match = 0;
		core::Real sum_sq = 0;
		for ( core::Size ii = 1; ii <= pose.size(); ++ii ) {
			auto map_pair = target_pose_map_.find( ii );
			if ( map_pair != target_pose_map_.end() && map_pair->second <= target_pose_.size() ) {
				sum_sq += pose.residue( ii ).xyz( "CA" ).distance_squared( target_pose_.residue( map_pair->second ).xyz( "CA" ) );
				num_match++;
			}
		}
	}

	if ( num_match > 0 ) {
		sum_sq = std::sqrt( sum_sq/num_match );
	}

	return sum_sq;

}


void
DomainAssemblyMover::initialize()
{
	if ( ! movemap_set() ) initialize_movemap_from_commandline();
	if ( ! fragsets_set() ) initialize_fragments_from_commandline();
	if ( basic::options::option [ basic::options::OptionKeys::DomainAssembly::da_eval_pose_map ].user() ) {
		initialize_pose_map_from_commandline();
	}
	if ( target_pose_.empty() && basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) initialize_target_pose();
	if ( starting_pose_.empty() && basic::options::option[ basic::options::OptionKeys::DomainAssembly::da_start_pdb ].user() ) {
		initialize_start_pose_from_commandline();
	}
	if ( buried_residues_.empty() && basic::options::option[ basic::options::OptionKeys::DomainAssembly::da_require_buried ].user() ) {
		initialize_buried_from_commandline();
	}
}

void DomainAssemblyMover::initialize_target_pose() {
	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		core::import_pose::pose_from_file( target_pose_,
			*(core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )),
			basic::options::option[ basic::options::OptionKeys::in::file::native ](), false, core::import_pose::PDB_file );
	} else {
		throw utility::excn::EXCN_Msg_Exception( "No native pose specified! Use -in::file::native\n" );
	}
}

void DomainAssemblyMover::initialize_pose_map_from_commandline() {
	target_pose_map_.clear();
	std::string filename( basic::options::option[ basic::options::OptionKeys::DomainAssembly::da_eval_pose_map ]() );
	utility::io::izstream data( filename );
	if ( !data ) {
		std::stringstream err_msg;
		err_msg << " cannot open map file: " << filename << std::endl;
		throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
	}

	std::string line;
	while ( getline( data, line ) ) {
		std::istringstream line_stream( line );
		core::Size pose_number, native_pose_number;
		line_stream >> pose_number >> native_pose_number;
		if ( line_stream.fail() ) {
			std::stringstream err_msg;
			err_msg << " cannot parse line in map file: " << line << std::endl;
			throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
		}
		target_pose_map_.insert( std::make_pair( pose_number, native_pose_number ) );
	}

	data.close();
	data.clear();
}

void DomainAssemblyMover::initialize_buried_from_commandline() {
	buried_residues_.clear();
	std::string filename( basic::options::option[ basic::options::OptionKeys::DomainAssembly::da_require_buried ]() );
	utility::io::izstream data( filename );
	if ( !data ) {
		std::stringstream err_msg;
		err_msg << " cannot open buried residue definition file: " << filename << std::endl;
		throw utility::excn::EXCN_Msg_Exception( err_msg.str() );
	}

	std::stringstream buried_string;
	std::string line;
	while ( getline( data, line ) ) {
		buried_string << line << ",";
	}

	buried_residues_ = buried_string.str();
	data.close();
	data.clear();
}

void DomainAssemblyMover::initialize_start_pose_from_commandline() {
	core::import_pose::pose_from_file( starting_pose_, basic::options::option[ basic::options::OptionKeys::DomainAssembly::da_start_pdb ]() , core::import_pose::PDB_file);
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
	core::fragment::ConstantLengthFragSetOP fragset3mer = nullptr;
	core::fragment::ConstantLengthFragSetOP fragset9mer = nullptr;
	if ( basic::options::option[ basic::options::OptionKeys::in::file::frag3].user() ) {
		fragset3mer = core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( 3 ) );
		fragset3mer->read_fragment_file( basic::options::option[ basic::options::OptionKeys::in::file::frag3 ]() );
		fragset3mer_ = fragset3mer;
	}
	if ( basic::options::option[ basic::options::OptionKeys::in::file::frag9].user() ) {
		fragset9mer = core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( 9 ) );
		fragset9mer->read_fragment_file( basic::options::option[ basic::options::OptionKeys::in::file::frag9 ]() );
		fragset9mer_ = fragset9mer;
	}
}

void DomainAssemblyMover::recover_sidechains( core::pose::Pose & pose ) const
{
	if ( starting_pose_.empty() ) {
		throw utility::excn::EXCN_Msg_Exception( "No starting pose specified! Use -da_start_pdb\n" );
	}
	protocols::simple_moves::ReturnSidechainMover return_sidechains( starting_pose_ );
	return_sidechains.apply( pose );
}

} // namespace domain_assembly
} // namespace devel

