// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/FastRelax.bench.cc
/// @brief  Performance benchmark running the fullatom fast-relax protocol.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_apps_benchmark_FastRelax_bench_hh
#define INCLUDED_apps_benchmark_FastRelax_bench_hh


#include <apps/benchmark/performance/performance_benchmark.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>

#include <protocols/relax/FastRelax.hh>

//#include <core/scoring/ScoringManager.hh>
//#include <core/scoring/ScoreTypeManager.hh>
//#include <core/scoring/Energies.hh>
//#include <core/scoring/EnergyGraph.hh>

//#include <core/scoring/methods/OneBodyEnergy.hh>
//#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>
//#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
//#include <core/scoring/LREnergyContainer.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>


enum fast_relax_perf_benchmark {
	fast_relax_perfbench_score12,
	fast_relax_perfbench_sc12sp2,
	fast_relax_perfbench_sc12he,
	fast_relax_perfbench_mmstd,
	fast_relax_perfbench_sp2hecart
};

class FastRelaxPerformanceBenchmark : public PerformanceBenchmark
{

public:
	FastRelaxPerformanceBenchmark(
		std::string name,
		fast_relax_perf_benchmark benchtype
	) :
		PerformanceBenchmark( name ),
		benchtype_( benchtype )
	{}

	virtual void setUp() {
		pose_ = core::pose::PoseOP( new core::pose::Pose );
		// Use smaller pdb to test relax
		core::import_pose::pose_from_pdb(*pose_, "test_in2.pdb");
		scorefxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
		switch ( benchtype_ ) {
		case fast_relax_perfbench_score12 :
			setup_for_score12();
			break;
		case fast_relax_perfbench_sc12sp2 :
			setup_for_sc12sp2();
			break;
		case fast_relax_perfbench_sc12he :
			setup_for_sc12he();
			break;
		case fast_relax_perfbench_mmstd :
			setup_for_mmstd();
			break;
		case fast_relax_perfbench_sp2hecart :
			setup_for_sp2hecart();
			break;
		}
	}

	virtual void run(core::Real scaleFactor) {
		core::Size reps( (core::Size)(1 * scaleFactor) );
		if ( reps == 0 ) { reps = 1; } // do at least one repetition, regardless of scale factor.
		for ( core::Size ii = 1; ii <= reps; ++ii ) {
			core::pose::Pose runpose = *pose_; // don't start from the last iteration's pose
			fr_->apply( runpose );
		}
	}

	virtual void tearDown() {
		pose_.reset();
		scorefxn_.reset();
		fr_.reset();
	}

	void setup_for_score12() {
		using namespace core::scoring;
		using namespace protocols::relax;

		scorefxn_->set_weight( fa_atr, 0.8 );
		scorefxn_->set_weight( fa_rep, 0.44 );
		scorefxn_->set_weight( fa_sol, 0.65 );
		scorefxn_->set_weight( fa_intra_rep, 0.004 );
		scorefxn_->set_weight( fa_pair, 0.49 );
		scorefxn_->set_weight( fa_plane, 0 );
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

		fr_ = protocols::relax::FastRelaxOP( new FastRelax( scorefxn_ ) );
	}

	void setup_for_sc12sp2() {
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::scoring::hbonds;
		using namespace protocols::relax;

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
		scorefxn_->set_weight( fa_plane, 0 );
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

		fr_ = protocols::relax::FastRelaxOP( new FastRelax( scorefxn_ ) );
	}

	void setup_for_sc12he() {
		using namespace core::scoring;
		using namespace protocols::relax;

		scorefxn_->set_weight( fa_atr, 0.8 );
		scorefxn_->set_weight( fa_rep, 0.44 );
		scorefxn_->set_weight( fa_sol, 0.65 );
		scorefxn_->set_weight( fa_intra_rep, 0.004 );
		scorefxn_->set_weight( fa_pair, 0.49 );
		scorefxn_->set_weight( fa_plane, 0 );
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

		fr_ = protocols::relax::FastRelaxOP( new FastRelax( scorefxn_ ) );

	}

	void setup_for_mmstd() {
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace protocols::relax;

		EnergyMethodOptions emo;
		emo.unfolded_energies_type( "UNFOLDED_MM_STD" );
		scorefxn_->set_energy_method_options( emo );

		scorefxn_->set_weight( fa_atr, 0.8 );
		scorefxn_->set_weight( fa_rep, 0.634454 );
		scorefxn_->set_weight( fa_sol, 1.16497 );
		//scorefxn_->set_weight( mm_lj_intra_rep, 0.324341 ); // UNCOMMENT THESE WHEN THEY STOP NOT WORKING
		//scorefxn_->set_weight( mm_lj_intra_atr, 0.537815 );
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

		fr_ = protocols::relax::FastRelaxOP( new FastRelax( scorefxn_ ) );
	}

	void setup_for_sp2hecart() {
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::scoring::hbonds;
		using namespace protocols::relax;

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
		scorefxn_->set_weight( fa_plane, 0 );
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

		fr_ = protocols::relax::FastRelaxOP( new FastRelax( scorefxn_ ) );
		fr_->cartesian( true );

	}

private:
	fast_relax_perf_benchmark benchtype_;
	core::pose::PoseOP pose_;
	core::scoring::ScoreFunctionOP scorefxn_;
	protocols::relax::FastRelaxOP fr_;
};

FastRelaxPerformanceBenchmark frpb_sc12( "protocols.relax.fast_relax_performance_benchmark_score12", fast_relax_perfbench_score12 );
FastRelaxPerformanceBenchmark frpb_sc12sp2( "protocols.relax.fast_relax_performance_benchmark_sc12sp2", fast_relax_perfbench_sc12sp2 );
FastRelaxPerformanceBenchmark frpb_sc12he( "protocols.relax.fast_relax_performance_benchmark_sc12he", fast_relax_perfbench_sc12he );
FastRelaxPerformanceBenchmark frpb_mmstd( "protocols.relax.fast_relax_performance_benchmark_mmstd", fast_relax_perfbench_mmstd );
//FastRelaxPerformanceBenchmark frpb_sp2hecart( "protocols.relax.fast_relax_performance_benchmark_sp2hecart", fast_relax_perfbench_sp2hecart );

#endif
