// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   rosetta/benchmark/InteractionGraph.bench.cc
/// @brief  Performance benchmark for testing the speed of the interaction graphs in simulated annealing
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_apps_benchmark_InteractionGraph_bench_hh
#define INCLUDED_apps_benchmark_InteractionGraph_bench_hh


#include <apps/benchmark/performance/performance_benchmark.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/graph/Graph.fwd.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#include <core/pack/packer_neighbors.hh>

#include <core/pack/annealer/DebuggingAnnealer.hh>

#include <core/pack/interaction_graph/DensePDInteractionGraph.hh>
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.hh>
#include <core/pack/interaction_graph/PDInteractionGraph.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/mistakes.OptionKeys.gen.hh>
#include <basic/options/keys/multithreading.OptionKeys.gen.hh>


#include <basic/Tracer.hh>

#include <core/pack/rotamer_set/RotamerSets.hh> // AUTO IWYU For RotamerSets
#include <core/pack/task/ResidueLevelTask.hh> // AUTO IWYU For ResidueLevelTask


static basic::Tracer TR2( "apps.benchmark.performance.InteractionGraphPerformanceBenchmark" );

enum interaction_graph_perf_benchmark {
	interaction_graph_perfbench_linmemig_score12,
	interaction_graph_perfbench_linmemig_sc12sp2,
	interaction_graph_perfbench_linmemig_sc12he,
	interaction_graph_perfbench_linmemig_mmstd,
	interaction_graph_perfbench_pdig_score12,
	interaction_graph_perfbench_denseig_score12,
	interaction_graph_perfbench_linmemig_current_default_sfxn,
	interaction_graph_perfbench_pdig_current_default_sfxn,
	interaction_graph_perfbench_denseig_current_default_sfxn,
};

class InteractionGraphPerformanceBenchmark : public PerformanceBenchmark
{

public:
	InteractionGraphPerformanceBenchmark(
		std::string name,
		interaction_graph_perf_benchmark benchtype,
		core::Size base_scale = 1
	) :
		PerformanceBenchmark( name ),
		benchtype_( benchtype ),
		base_scale_( base_scale )
	{}

	void setUp() override {
		configured_for_score12_ = basic::options::option[ basic::options::OptionKeys::mistakes::restore_pre_talaris_2013_behavior ]();
		pose_ = utility::pointer::make_shared< core::pose::Pose >();
		core::import_pose::pose_from_file(*pose_, "test_in2.pdb", core::import_pose::PDB_file);
		scorefxn_ = utility::pointer::make_shared< core::scoring::ScoreFunction >();
		switch ( benchtype_ ) {
		case interaction_graph_perfbench_linmemig_score12 :
			setup_for_score12();
			trajectory_fname_ = "interaction_graph_perfbench_linmemig_score12.traj.gz";
			setup_for_linmemig();
			break;
		case interaction_graph_perfbench_linmemig_sc12sp2 :
			setup_for_sc12sp2();
			trajectory_fname_ = "interaction_graph_perfbench_linmemig_sc12sp2.traj";
			setup_for_linmemig();
			break;
		case interaction_graph_perfbench_linmemig_sc12he :
			setup_for_sc12he();
			trajectory_fname_ = "interaction_graph_perfbench_linmemig_sc12he.traj";
			setup_for_linmemig();
			break;
		case interaction_graph_perfbench_linmemig_mmstd :
			setup_for_mmstd();
			trajectory_fname_ = "interaction_graph_perfbench_linmemig_mmstd.traj";
			setup_for_linmemig();
			break;
		case interaction_graph_perfbench_pdig_score12 :
			setup_for_score12();
			trajectory_fname_ = "interaction_graph_perfbench_pdig_score12.traj";
			setup_for_pdig();
			break;
		case interaction_graph_perfbench_denseig_score12 :
			setup_for_score12();
			trajectory_fname_ = "interaction_graph_perfbench_denseig_score12.traj";
			setup_for_denseig();
			break;
		case interaction_graph_perfbench_linmemig_current_default_sfxn :
			setup_for_current_default_sfxn();
			trajectory_fname_ = "interaction_graph_perfbench_linmemig_current_default_sfxn.traj";
			setup_for_linmemig();
			break;
		case interaction_graph_perfbench_pdig_current_default_sfxn :
			setup_for_current_default_sfxn();
			trajectory_fname_ = "interaction_graph_perfbench_pdig_current_default_sfxn.traj";
			setup_for_pdig();
			break;
		case interaction_graph_perfbench_denseig_current_default_sfxn :
			setup_for_current_default_sfxn();
			trajectory_fname_ = "interaction_graph_perfbench_denseig_current_default_sfxn.traj";
			setup_for_denseig();
			break;
		default :
			utility_exit_with_message( "Error in InteractionGraphPerformanceBenchmark::setUp(): An unrecognized benchtype_ value was passed to the constructor." );
		}
	}

	void run(core::Real scaleFactor) override {
		using namespace core;
		using namespace core::pack;
		using namespace core::pack::annealer;
		using namespace core::pack::interaction_graph;
		using namespace core::pack::rotamer_set;

		PackerEnergy bestenergy;
		FixbbRotamerSetsCOP rsets = utility::pointer::static_pointer_cast< core::pack::rotamer_set::FixbbRotamerSets const > ( rotsets_ );
		ObjexxFCL::FArray1D_int bestrotamer_at_seqpos( pose_->size(), 0 );
		ObjexxFCL::FArray1D_int current_rot_index( pose_->size(), 0 );
		ObjexxFCL::FArray1D_float rot_freq( rotsets_->nrotamers(), 0.0f );

		/// SETUP RUN CODE
		//FixbbSimAnnealer fixbbsa(
		// bestrotamer_at_seqpos,
		// bestenergy,
		// false,
		// ig_,
		// rsets,
		// current_rot_index,
		// false,
		// rot_freq
		//);
		//fixbbsa.record_annealer_trajectory( true );
		//fixbbsa.trajectory_file_name( trajectory_fname_ );
		//fixbbsa.run();
		/// END SETUP RUN

		DebuggingAnnealer dba(
			bestrotamer_at_seqpos,
			bestenergy,
			false,
			ig_,
			rsets,
			current_rot_index,
			false,
			rot_freq
		);

		dba.annealer_file( trajectory_fname_ );

		core::Size reps( (core::Size)(base_scale_ * scaleFactor) );
		if ( reps == 0 ) { reps = 1; } // Do at least one rep, regardless of scaling factor
		for ( core::Size ii = 1; ii <= reps; ++ii ) {
			dba.run(); // the main piece of code to benchmark
		}
	}

	void tearDown() override {
		pose_.reset();
		scorefxn_.reset();
		ig_.reset();
		task_.reset();
		rotsets_.reset();
		packer_neighbor_graph_.reset();
	}

	void setup_for_score12() {
		using namespace core::scoring;
		if ( !configured_for_score12_ ) {
			TR2.Warning << "Warning!  Benchmarks were run without the \"-restore_pre_talaris_behavior\" flag.  Score12 weights will be used, but this will NOT represent true scoring with the score12 scoring function." << std::endl;
		}

		scorefxn_->set_weight( fa_atr, 0.8 );
		scorefxn_->set_weight( fa_rep, 0.44 );
		scorefxn_->set_weight( fa_sol, 0.65 );
		scorefxn_->set_weight( fa_intra_rep, 0.004 );
		scorefxn_->set_weight( fa_pair, 0.49 );
		//scorefxn_->set_weight( fa_plane, 0 );
		scorefxn_->set_weight( fa_dun, 0.56 );
		scorefxn_->set_weight( ref, 1 );
		scorefxn_->set_weight( hbond_lr_bb, 1.17 );
		scorefxn_->set_weight( hbond_sr_bb, 0.585 );
		scorefxn_->set_weight( hbond_bb_sc, 1.17 );
		scorefxn_->set_weight( hbond_sc, 1.1 );
		scorefxn_->set_weight( p_aa_pp, 0.32 );
		scorefxn_->set_weight( dslf_ss_dst, 0.5 );
		scorefxn_->set_weight( dslf_cs_ang, 2.0 );
		scorefxn_->set_weight( dslf_ss_dih, 5.0 );
		scorefxn_->set_weight( dslf_ca_dih, 5.0 );
		scorefxn_->set_weight( pro_close, 1.0 );
		scorefxn_->set_weight( omega, 0.5 );
		scorefxn_->set_weight( rama, 0.2 );
	}

	void setup_for_current_default_sfxn() {
		runtime_assert_string_msg( !configured_for_score12_, "Error in apps::benchmark::performance::InteractionGraphPerformanceBenchmark::setup_for_current_default_sfxn(): The current default scoring function is incompatible with the \"-restore_pre_talaris_behaviour\" flag." );
		using namespace core::scoring;
		scorefxn_ = core::scoring::get_score_function();
	}

	void setup_for_sc12sp2() {
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::scoring::hbonds;

		if ( !configured_for_score12_ ) {
			TR2.Warning << "Warning!  Benchmarks were run without the \"-restore_pre_talaris_behavior\" flag.  The sc12sp2 weights will be used, but this will NOT represent true scoring with the sc12sp2 scoring function." << std::endl;
		}

		EnergyMethodOptions emo;
		HBondOptions hbo;
		hbo.use_sp2_chi_penalty( true );
		hbo.measure_sp3acc_BAH_from_hvy( true );
		hbo.fade_energy( true );
		emo.hbond_options( hbo );
		scorefxn_->set_energy_method_options( emo );

		scorefxn_->set_weight( fa_atr, 0.8 );
		scorefxn_->set_weight( fa_rep, 0.44 );
		scorefxn_->set_weight( fa_sol, 0.65 );
		scorefxn_->set_weight( fa_intra_rep, 0.004 );
		scorefxn_->set_weight( fa_pair, 0.49 );
		//scorefxn_->set_weight( fa_plane, 0 );
		scorefxn_->set_weight( fa_dun, 0.56 );
		scorefxn_->set_weight( ref, 1 );
		scorefxn_->set_weight( hbond_lr_bb, 1.17 );
		scorefxn_->set_weight( hbond_sr_bb, 1.17 );
		scorefxn_->set_weight( hbond_bb_sc, 1.17 );
		scorefxn_->set_weight( hbond_sc, 1.1 );
		scorefxn_->set_weight( p_aa_pp, 0.32 );
		scorefxn_->set_weight( dslf_ss_dst, 0.5 );
		scorefxn_->set_weight( dslf_cs_ang, 2.0 );
		scorefxn_->set_weight( dslf_ss_dih, 5.0 );
		scorefxn_->set_weight( dslf_ca_dih, 5.0 );
		scorefxn_->set_weight( pro_close, 1.0 );
		scorefxn_->set_weight( omega, 0.5 );
		scorefxn_->set_weight( rama, 0.2 );

	}

	void setup_for_sc12he() {
		using namespace core::scoring;


		if ( !configured_for_score12_ ) {
			TR2.Warning << "Warning!  Benchmarks were run without the \"-restore_pre_talaris_behavior\" flag.  The sc12he weights will be used, but this will NOT represent true scoring with the sc12he scoring function." << std::endl;
		}

		scorefxn_->set_weight( fa_atr, 0.8 );
		scorefxn_->set_weight( fa_rep, 0.44 );
		scorefxn_->set_weight( fa_sol, 0.65 );
		scorefxn_->set_weight( fa_intra_rep, 0.004 );
		scorefxn_->set_weight( fa_pair, 0.49 );
		//scorefxn_->set_weight( fa_plane, 0 );
		scorefxn_->set_weight( fa_dun, 0.56 );
		scorefxn_->set_weight( ref, 1 );
		scorefxn_->set_weight( hbond_lr_bb, 1.17 );
		scorefxn_->set_weight( hbond_sr_bb, 0.585 );
		scorefxn_->set_weight( hbond_bb_sc, 1.17 );
		scorefxn_->set_weight( hbond_sc, 1.1 );
		scorefxn_->set_weight( p_aa_pp, 0.32 );
		scorefxn_->set_weight( dslf_ss_dst, 0.5 );
		scorefxn_->set_weight( dslf_cs_ang, 2.0 );
		scorefxn_->set_weight( dslf_ss_dih, 5.0 );
		scorefxn_->set_weight( dslf_ca_dih, 5.0 );
		scorefxn_->set_weight( pro_close, 1.0 );
		scorefxn_->set_weight( omega, 0.5 );
		scorefxn_->set_weight( rama, 0.2 );
		scorefxn_->set_weight( fa_elec, 0.7 );
	}

	void setup_for_mmstd() {
		using namespace core::scoring;
		using namespace core::scoring::methods;

		EnergyMethodOptions emo;
		emo.unfolded_energies_type( "UNFOLDED_MM_STD" );
		scorefxn_->set_energy_method_options( emo );

		scorefxn_->set_weight( fa_atr, 0.8 );
		scorefxn_->set_weight( fa_rep, 0.634454 );
		scorefxn_->set_weight( fa_sol, 1.16497 );
		scorefxn_->set_weight( mm_lj_intra_rep, 0.324341 );
		scorefxn_->set_weight( mm_lj_intra_atr, 0.537815 );
		scorefxn_->set_weight( mm_twist, 0.2662 );
		scorefxn_->set_weight( pro_close, 1.44777 );
		scorefxn_->set_weight( hbond_sr_bb, 0.656728 );
		scorefxn_->set_weight( hbond_lr_bb, 1.50186 );
		scorefxn_->set_weight( hbond_bb_sc, 1.45367 );
		scorefxn_->set_weight( hbond_sc, 1.18477 );
		scorefxn_->set_weight( dslf_ss_dst, 0.5 );
		scorefxn_->set_weight( dslf_cs_ang, 2 );
		scorefxn_->set_weight( dslf_ss_dih, 5 );
		scorefxn_->set_weight( dslf_ca_dih, 5 );
		scorefxn_->set_weight( unfolded, -0.904283 );

	}

	void setup_for_sp2hecart() {
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::scoring::hbonds;

		if ( !configured_for_score12_ ) {
			TR2.Warning << "Warning!  Benchmarks were run without the \"-restore_pre_talaris_behavior\" flag.  The sp2hecart weights will be used, but this will NOT represent true scoring with the sp2hecart scoring function." << std::endl;
		}


		EnergyMethodOptions emo;
		HBondOptions hbo;
		hbo.use_sp2_chi_penalty( true );
		hbo.measure_sp3acc_BAH_from_hvy( true );
		hbo.fade_energy( true );
		emo.hbond_options( hbo );
		scorefxn_->set_energy_method_options( emo );

		scorefxn_->set_weight( fa_atr, 0.8 );
		scorefxn_->set_weight( fa_rep, 0.44 );
		scorefxn_->set_weight( fa_sol, 0.65 );
		scorefxn_->set_weight( fa_intra_rep, 0.004 );
		scorefxn_->set_weight( fa_pair, 0.49 );
		//scorefxn_->set_weight( fa_plane, 0 );
		scorefxn_->set_weight( fa_dun, 0.56 );
		scorefxn_->set_weight( ref, 1 );
		scorefxn_->set_weight( hbond_lr_bb, 1.17 );
		scorefxn_->set_weight( hbond_sr_bb, 1.17 );
		scorefxn_->set_weight( hbond_bb_sc, 1.17 );
		scorefxn_->set_weight( hbond_sc, 1.1 );
		scorefxn_->set_weight( p_aa_pp, 0.32 );
		scorefxn_->set_weight( dslf_ss_dst, 0.5 );
		scorefxn_->set_weight( dslf_cs_ang, 2.0 );
		scorefxn_->set_weight( dslf_ss_dih, 5.0 );
		scorefxn_->set_weight( dslf_ca_dih, 5.0 );
		scorefxn_->set_weight( pro_close, 1.0 );
		scorefxn_->set_weight( omega, 0.5 );
		scorefxn_->set_weight( rama, 0.2 );
		scorefxn_->set_weight( fa_elec, 0.7 );
		scorefxn_->set_weight( cart_bonded, 0.5 );
	}

	core::pack::task::PackerTaskOP
	redesign_20() {
		// create a packer task to redesign 20 residues on test_in.pdb
		using namespace core;
		using namespace core::pack::task;
		PackerTaskOP task = TaskFactory::create_packer_task( *pose_ );
		for ( Size ii = 21; ii <= pose_->size(); ++ii ) {
			task->nonconst_residue_task( ii ).prevent_repacking();
		}
		for ( Size ii = 1; ii <= 20; ++ii ) {
			task->nonconst_residue_task( ii ).or_ex1( true );
			task->nonconst_residue_task( ii ).or_ex2( true );
		}
		return task;
	}

	void prepare_rotamer_sets() {
		using namespace core;
		using namespace core::pack;
		using namespace core::pack::rotamer_set;
		using namespace core::pack::task;
		(*scorefxn_)( *pose_ );

		scorefxn_->setup_for_packing(
			*pose_,
			task_->repacking_residues(),
			task_->designing_residues()
		);
		packer_neighbor_graph_ = create_packer_graph( *pose_, *scorefxn_, task_ );
		rotsets_ = utility::pointer::make_shared< RotamerSets >();
		rotsets_->set_task( task_ );
		rotsets_->build_rotamers( *pose_, *scorefxn_, packer_neighbor_graph_ );
		rotsets_->prepare_sets_for_packing( *pose_, *scorefxn_ );
	}

	void setup_for_linmemig() {
		using namespace core::pack::interaction_graph;
		task_ = redesign_20();
		prepare_rotamer_sets();

		LinearMemoryInteractionGraphOP lmig( new LinearMemoryInteractionGraph( task_->num_to_be_packed() ) );
		lmig->set_pose( *pose_ );
		lmig->set_score_function( *scorefxn_ );
		lmig->set_recent_history_size( task_->linmem_ig_history_size() );
		ig_ = lmig;
		rotsets_->compute_energies( *pose_, *scorefxn_, packer_neighbor_graph_, ig_, basic::options::option[ basic::options::OptionKeys::multithreading::interaction_graph_threads ]() );
	}

	void setup_for_pdig() {
		using namespace core::pack::interaction_graph;
		task_ = redesign_20();
		prepare_rotamer_sets();
		ig_ = utility::pointer::make_shared< PDInteractionGraph >( task_->num_to_be_packed() );
		rotsets_->compute_energies( *pose_, *scorefxn_, packer_neighbor_graph_, ig_, basic::options::option[ basic::options::OptionKeys::multithreading::interaction_graph_threads ]() );
	}

	void setup_for_denseig() {
		using namespace core;
		using namespace core::pack::interaction_graph;
		task_ = redesign_20();
		for ( Size ii = 1; ii <= 20; ++ii ) task_->nonconst_residue_task( ii ).restrict_to_repacking();
		prepare_rotamer_sets();
		ig_ = utility::pointer::make_shared< DensePDInteractionGraph >( task_->num_to_be_packed() );
		rotsets_->compute_energies( *pose_, *scorefxn_, packer_neighbor_graph_, ig_, basic::options::option[ basic::options::OptionKeys::multithreading::interaction_graph_threads ]() );
	}

private:
	interaction_graph_perf_benchmark benchtype_;
	core::pose::PoseOP pose_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::pack::interaction_graph::InteractionGraphBaseOP ig_;
	core::pack::task::PackerTaskOP task_;
	core::pack::rotamer_set::RotamerSetsOP rotsets_;
	utility::graph::GraphOP packer_neighbor_graph_;
	std::string trajectory_fname_;
	core::Size base_scale_;
	bool configured_for_score12_ = false;
};

InteractionGraphPerformanceBenchmark igpb_lmig_sc12( "core.pack.linmem_ig_score12", interaction_graph_perfbench_linmemig_score12, 1 );
InteractionGraphPerformanceBenchmark igpb_lmig_sc12sp2( "core.pack.linmem_ig_sc12sp2", interaction_graph_perfbench_linmemig_sc12sp2, 1 );
InteractionGraphPerformanceBenchmark igpb_lmig_sc12he( "core.pack.linmem_ig_sc12he", interaction_graph_perfbench_linmemig_sc12he, 1 );
InteractionGraphPerformanceBenchmark igpb_lmig_mmstd( "core.pack.linmem_ig_mmstd", interaction_graph_perfbench_linmemig_mmstd, 1 );
InteractionGraphPerformanceBenchmark igpb_lmig_current_default_sfxn( "core.pack.linmem_ig_current_default_sfxn", interaction_graph_perfbench_linmemig_current_default_sfxn, 1 );
InteractionGraphPerformanceBenchmark igpb_pdig_sc12( "core.pack.pdig_score12", interaction_graph_perfbench_pdig_score12, 4 );
InteractionGraphPerformanceBenchmark igpb_pdig_current_default_sfxn( "core.pack.pdig_current_default_sfxn", interaction_graph_perfbench_pdig_current_default_sfxn, 4 );
InteractionGraphPerformanceBenchmark igpb_denseig_sc12( "core.pack.denseig_score12", interaction_graph_perfbench_denseig_score12, 7000 );
InteractionGraphPerformanceBenchmark igpb_denseig_current_default_sfxn( "core.pack.denseig_current_default_sfxn", interaction_graph_perfbench_denseig_current_default_sfxn, 7000 );

#endif
