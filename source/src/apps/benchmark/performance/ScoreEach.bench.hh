// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   rosetta/benchmark/score.bench.cc
///
/// @brief  Scoring Each benchmark
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_apps_benchmark_ScoreEach_bench_hh
#define INCLUDED_apps_benchmark_ScoreEach_bench_hh


#include <apps/benchmark/performance/performance_benchmark.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>

#include <core/scoring/methods/OneBodyEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>
#include <core/scoring/methods/LongRangeTwoBodyEnergy.hh>
#include <core/scoring/LREnergyContainer.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>


class ScoreEachBenchmark : public PerformanceBenchmark
{

public:
	ScoreEachBenchmark(
		std::string name,
		core::scoring::ScoreType score_type,
		core::Size base_scale_factor,
		std::string infile = "test_in.pdb",
		bool centroid = false
	) :
		PerformanceBenchmark(name),
		score_type_(score_type),
		base_scale_factor_(base_scale_factor),
		setup_successful_(false),
		in_file_(infile),
		centroid_(centroid)
	{}


	virtual
	void
	set_scorefxn(
		core::scoring::ScoreFunctionOP scorefxn
	) {
		scorefxn_ = scorefxn;
	}

	virtual void setUp() {
		pose_ = core::pose::PoseOP( new core::pose::Pose );
		if ( centroid_ ) {
			core::chemical::ResidueTypeSetCOP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::CENTROID );
			core::import_pose::pose_from_file(*pose_, *rts, in_file_, core::import_pose::PDB_file);
		} else {
			core::import_pose::pose_from_file(*pose_, in_file_, core::import_pose::PDB_file);
		}
		core::scoring::methods::EnergyMethodOptions opts;
		try {
			enmeth_ = core::scoring::ScoringManager::get_instance()->energy_method( score_type_, opts );
		} catch (utility::excn::EXCN_Base& excn){
			TR.Error
				<< "Unable to Test scoring with score type '"
				<< core::scoring::ScoreTypeManager::name_from_score_type(score_type_)
				<< "'" << std::endl
				<< "Fail with error:" << std::endl
				<< excn;
		}
		enmethtype_ = enmeth_->method_type();

		scorefxn_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
		try{
			// do this once in case there are one time setup requirements
			scorefxn_->set_weight(score_type_, 1);
			// Run for testing of score value *and* side-effects of calling
			// Exact equality here is intended. Close-to-zero implies substantial calculations
			if ( scorefxn_->score(*pose_) == 0 ) {
				TR.Warning << "WARNING: Performance test for score term '" << score_type_
					<< "' gives zero energy with test structure. Timing values may not be informative." << std::endl;
			}
			setup_successful_ = true;
		} catch (utility::excn::EXCN_Base& excn){
			TR.Error
				<< "Unable to Test scoring with score type '"
				<< core::scoring::ScoreTypeManager::name_from_score_type(score_type_)
				<< "'" << std::endl
				<< "Fail with error:" << std::endl
				<< excn;
		}
	}

	virtual void run(core::Real scaleFactor) {
		using namespace core::scoring::methods;

		if ( !setup_successful_ ) {
			TR.Error
				<< "Skipping running this test becuase the setup was not completed successfully." << std::endl;
			return;
		}

		if ( enmethtype_ == ci_1b || enmethtype_ == cd_1b ) {
			run_one_body_energy( scaleFactor );
		} else if ( enmethtype_ == ci_2b || enmethtype_ == cd_2b ) {
			run_short_ranged_two_body_energy( scaleFactor );
		} else if ( enmethtype_ == ci_lr_2b || enmethtype_ == cd_lr_2b ) {
			run_long_ranged_two_body_energy( scaleFactor );
		} else {
			// whole structure energies
			core::scoring::ScoreFunction scorefxn;
			try{
				core::Size reps( (core::Size)(base_scale_factor_*scaleFactor) );
				if ( reps == 0 ) { reps = 1; } // Do at least one repetion, regardless of scaling.
				for ( core::Size i=0; i < reps; i++ ) {
					scorefxn.set_weight(score_type_, 1);
					scorefxn.score(*pose_);
					pose_->energies().clear();
					scorefxn.set_weight(score_type_, 0);
				}
			} catch (utility::excn::EXCN_Base& excn){
				TR.Error
					<< "Unable to Test scoring with score type '"
					<< core::scoring::ScoreTypeManager::name_from_score_type(score_type_)
					<< "'" << std::endl
					<< "Fail with error:" << std::endl
					<< excn;
			}
		}
	}

	void run_one_body_energy( core::Real scaleFactor ) {
		using namespace core;
		using namespace core::scoring;
		using namespace core::scoring::methods;

		utility::vector1< core::conformation::Residue const * > res( pose_->size() );
		for ( Size ii = 1; ii <= pose_->size(); ++ii ) { res[ ii ] = & pose_->residue( ii ); }

		if ( ! utility::pointer::dynamic_pointer_cast< core::scoring::methods::OneBodyEnergy const > ( enmeth_ ) ) {
			TR.Error
				<< "Unable to Test scoring with score type '"
				<< core::scoring::ScoreTypeManager::name_from_score_type(score_type_)
				<< "'" << std::endl
				<< "Failed dynamic cast of energy method to OneBodyEnergy" << std::endl;
		}

		OneBodyEnergy const & e1b( static_cast< OneBodyEnergy const & > ( *enmeth_ ));
		EnergyMap emap;
		core::Size reps( (core::Size)(base_scale_factor_*scaleFactor) );
		if ( reps == 0 ) { reps = 1; } // Do at least one repetion, regardless of scaling.
		for ( core::Size ii = 1; ii <= reps; ++ii ) {
			emap.zero( e1b.score_types() );
			for ( core::Size jj = 1; jj <= pose_->size(); ++jj ) {
				e1b.residue_energy( *res[ jj ], *pose_, emap );
			}
		}
	}

	void run_short_ranged_two_body_energy( core::Real scaleFactor ) {
		using namespace core;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::graph;

		utility::vector1< core::conformation::Residue const * > res( pose_->size() );
		for ( Size ii = 1; ii <= pose_->size(); ++ii ) { res[ ii ] = & pose_->residue( ii ); }

		if ( ! utility::pointer::dynamic_pointer_cast< core::scoring::methods::ShortRangeTwoBodyEnergy const > ( enmeth_ ) ) {
			TR.Error
				<< "Unable to Test scoring with score type '"
				<< core::scoring::ScoreTypeManager::name_from_score_type(score_type_)
				<< "'" << std::endl
				<< "Failed dynamic cast of energy method to ShortRangeTwoBodyEnergy" << std::endl;
		}

		ShortRangeTwoBodyEnergy const & e2b( static_cast< ShortRangeTwoBodyEnergy const & > (*enmeth_) );
		EnergyMap emap;
		core::Size reps( (core::Size)(base_scale_factor_*scaleFactor) );
		if ( reps == 0 ) { reps = 1; } // Do at least one repetion, regardless of scaling.
		for ( Size ii = 1; ii <= reps; ++ii ) {
			emap.zero( e2b.score_types() );
			for ( Size jj = 1; jj <= pose_->size(); ++jj ) {
				for ( Node::EdgeListConstIter
						iter     = pose_->energies().energy_graph().get_node( jj )->const_upper_edge_list_begin(),
						iter_end = pose_->energies().energy_graph().get_node( jj )->const_upper_edge_list_end();
						iter != iter_end; ++iter ) {
					Size kk = (*iter)->get_second_node_ind();
					e2b.residue_pair_energy( *res[ jj ], *res[ kk ], *pose_, *scorefxn_, emap );
				}
			}
		}
	}

	void run_long_ranged_two_body_energy( core::Real scaleFactor ) {
		using namespace core;
		using namespace core::scoring;
		using namespace core::scoring::methods;
		using namespace core::graph;
		if ( ! utility::pointer::dynamic_pointer_cast< core::scoring::methods::LongRangeTwoBodyEnergy const > ( enmeth_ ) ) {
			TR.Error
				<< "Unable to Test scoring with score type '"
				<< core::scoring::ScoreTypeManager::name_from_score_type(score_type_)
				<< "'" << std::endl
				<< "Failed dynamic cast of energy method to LongRangeTwoBodyEnergy" << std::endl;
		}

		utility::vector1< core::conformation::Residue const * > res( pose_->size() );
		for ( Size ii = 1; ii <= pose_->size(); ++ii ) { res[ ii ] = & pose_->residue( ii ); }

		LongRangeTwoBodyEnergy const & e2b( static_cast< LongRangeTwoBodyEnergy const & > (*enmeth_) );
		EnergyMap emap;
		LREnergyContainerOP lrec = pose_->energies().nonconst_long_range_container( e2b.long_range_type() );

		core::Size reps( (core::Size)(base_scale_factor_*scaleFactor) );
		if ( reps == 0 ) { reps = 1; } // Do at least one repetion, regardless of scaling.
		for ( Size ii = 1; ii <= reps; ++ii ) {
			emap.zero( e2b.score_types() );
			for ( Size jj = 1; jj <= pose_->size(); ++jj ) {
				if ( ! lrec->any_upper_neighbors_for_residue( jj ) ) continue;
				for ( ResidueNeighborIteratorOP
						rni = lrec->upper_neighbor_iterator_begin( jj ),
						rniend = lrec->upper_neighbor_iterator_end( jj );
						(*rni) != (*rniend); ++(*rni) ) {
					Size kk = rni->upper_neighbor_id();
					e2b.residue_pair_energy( *res[ jj ], *res[ kk ], *pose_, *scorefxn_, emap );
				}
			}
		}

	}

	virtual void tearDown() {
		pose_.reset();
		scorefxn_.reset();
		enmeth_.reset();
	}

private:
	core::pose::PoseOP pose_;
	core::scoring::ScoreType score_type_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::methods::EnergyMethodOP enmeth_;
	core::scoring::methods::EnergyMethodType enmethtype_;
	core::Size base_scale_factor_;
	bool setup_successful_;
	std::string in_file_;
	bool centroid_;
};

using namespace core::scoring;

// These energy terms don't actually calculate anything in the profiled methods
// they use finalize_total_energy() to actually do the calculation
//ScoreEachBenchmark Score_surface_("core.scoring.Score_1000x_surface",surface,1000);
//ScoreEachBenchmark Score_hpatch_("core.scoring.Score_10x_hpatch",hpatch,10);

// fast zeros under current setup - reactivate with a structure that results in reasonable calculations
//ScoreEachBenchmark Score_dna_chi_("core.scoring.Score_1000x_dna_chi",dna_chi,1000);
//ScoreEachBenchmark Score_rna_bulge_("core.scoring.Score_1000x_rna_bulge",rna_bulge,1000);
//ScoreEachBenchmark Score_nmer_ref_("core.scoring.Score_1000x_nmer_ref",nmer_ref,1000);
//ScoreEachBenchmark Score_nmer_pssm_("core.scoring.Score_1000x_nmer_pssm",nmer_pssm,1000);
//ScoreEachBenchmark Score_mm_lj_inter_rep_("core.scoring.Score_100x_mm_lj_inter_rep",mm_lj_inter_rep,100);
//ScoreEachBenchmark Score_mm_lj_inter_atr_("core.scoring.Score_100x_mm_lj_inter_atr",mm_lj_inter_atr,100);
//ScoreEachBenchmark Score_dof_constraint_("core.scoring.Score_1000x_dof_constraint",dof_constraint,1000);
//ScoreEachBenchmark Score_h2o_intra_("core.scoring.Score_1000x_h2o_intra",h2o_intra,1000);
//ScoreEachBenchmark Score_special_rot_("core.scoring.Score_1000x_special_rot",special_rot,1000);
//ScoreEachBenchmark Score_fa_stack_("core.scoring.Score_1000x_fa_stack",fa_stack,1000);
//ScoreEachBenchmark Score_fa_stack_aro_("core.scoring.Score_1000x_fa_stack_aro",fa_stack_aro,1000);

// So fast is it really worth benchmarking them?

ScoreEachBenchmark Score_ref_("core.scoring.Score_100000x_ref",ref,100000);
ScoreEachBenchmark Score_yhh_planarity_("core.scoring.Score_10000x_yhh_planarity",yhh_planarity,10000);
//centroid
ScoreEachBenchmark Score_cbeta_("core.scoring.Score_10000x_cbeta",cbeta,10000,"test_in_cen.pdb",true);
ScoreEachBenchmark Score_env_("core.scoring.Score_5000x_env",env,5000,"test_in_cen.pdb",true);

// really fast
ScoreEachBenchmark Score_fa_pair_("core.scoring.Score_100x_fa_pair",fa_pair,100);
ScoreEachBenchmark Score_fa_pair_aro_aro_("core.scoring.Score_100x_fa_pair_aro_aro",fa_pair_aro_aro,100);
ScoreEachBenchmark Score_fa_pair_aro_pol_("core.scoring.Score_100x_fa_pair_aro_pol",fa_pair_aro_pol,100);
ScoreEachBenchmark Score_fa_pair_pol_pol_("core.scoring.Score_100x_fa_pair_pol_pol",fa_pair_pol_pol,100);
ScoreEachBenchmark Score_omega_("core.scoring.Score_10000x_omega",omega,10000);
ScoreEachBenchmark Score_rama_("core.scoring.Score_10000x_rama",rama,10000);
ScoreEachBenchmark Score_rama2b_("core.scoring.Score_10000x_rama2b",rama2b,10000);
// not reduced
ScoreEachBenchmark Score_p_aa_pp_("core.scoring.Score_1000x_p_aa_pp",p_aa_pp,1000);
// readjusted
ScoreEachBenchmark Score_p_aa_("core.scoring.Score_40000x_p_aa",p_aa,40000);
ScoreEachBenchmark Score_mm_lj_intra_rep_("core.scoring.Score_1000x_mm_lj_intra_rep",mm_lj_intra_rep,1000);
ScoreEachBenchmark Score_mm_lj_intra_atr_("core.scoring.Score_1000x_mm_lj_intra_atr",mm_lj_intra_atr,1000);
ScoreEachBenchmark Score_unfolded_("core.scoring.Score_1000x_unfolded",unfolded,1000);
ScoreEachBenchmark Score_e_pH_("core.scoring.Score_1000x_e_pH",e_pH,1000);
//disulfides
ScoreEachBenchmark Score_dslf_ss_dst_("core.scoring.Score_1000x_dslf_ss_dst",dslf_ss_dst,1000,"1bbi_disulf.pdb");
ScoreEachBenchmark Score_dslf_cs_ang_("core.scoring.Score_1000x_dslf_cs_ang",dslf_cs_ang,1000,"1bbi_disulf.pdb");
ScoreEachBenchmark Score_dslf_ss_dih_("core.scoring.Score_1000x_dslf_ss_dih",dslf_ss_dih,1000,"1bbi_disulf.pdb");
ScoreEachBenchmark Score_dslf_ca_dih_("core.scoring.Score_1000x_dslf_ca_dih",dslf_ca_dih,1000,"1bbi_disulf.pdb");
ScoreEachBenchmark Score_dslf_cbs_ds_("core.scoring.Score_1000x_dslf_cbs_ds",dslf_cbs_ds,1000,"1bbi_disulf.pdb");
ScoreEachBenchmark Score_dslf_fa13_("core.scoring.Score_10000x_dslf_fa13",dslf_fa13,10000,"1bbi_disulf.pdb");
//centroid
ScoreEachBenchmark Score_cen_env_smooth_("core.scoring.Score_10000x_cen_env_smooth",cen_env_smooth,10000,"test_in_cen.pdb",true);
ScoreEachBenchmark Score_cbeta_smooth_("core.scoring.Score_10000x_cbeta_smooth",cbeta_smooth,10000,"test_in_cen.pdb",true);

// fast
ScoreEachBenchmark Score_pro_close_("core.scoring.Score_1000x_pro_close",pro_close,1000);
ScoreEachBenchmark Score_fa_dun_("core.scoring.Score_100x_fa_dun",fa_dun,100);
ScoreEachBenchmark Score_fa_dun_dev_("core.scoring.Score_100x_fa_dun_dev",fa_dun_dev,100);
ScoreEachBenchmark Score_fa_dun_rot_("core.scoring.Score_100x_fa_dun_rot",fa_dun_rot,100);
ScoreEachBenchmark Score_fa_dun_semi_("core.scoring.Score_100x_fa_dun_semi",fa_dun_semi,100);
//ref_nc currently gives zero energies, though it takes a while to do so
ScoreEachBenchmark Score_ref_nc_("core.scoring.Score_100x_ref_nc",ref_nc,100);
ScoreEachBenchmark Score_co_("core.scoring.Score_100x_co",co,100);
//chainbreak terms currently give zero energies, though it takes quite a while to do so
ScoreEachBenchmark Score_chainbreak_("core.scoring.Score_100x_chainbreak",chainbreak,100);
// not reduced
ScoreEachBenchmark Score_linear_chainbreak_("core.scoring.Score_100x_linear_chainbreak",linear_chainbreak,100);
ScoreEachBenchmark Score_overlap_chainbreak_("core.scoring.Score_100x_overlap_chainbreak",overlap_chainbreak,100);
ScoreEachBenchmark Score_distance_chainbreak_("core.scoring.Score_100x_distance_chainbreak",distance_chainbreak,100);
ScoreEachBenchmark Score_cart_bonded_("core.scoring.Score_100x_cart_bonded",cart_bonded,100);
ScoreEachBenchmark Score_mm_stretch_("core.scoring.Score_1000x_mm_stretch",mm_stretch,1000);
ScoreEachBenchmark Score_ch_bond_("core.scoring.Score_100x_ch_bond",ch_bond,100);
ScoreEachBenchmark Score_neigh_vect_("core.scoring.Score_100x_neigh_vect",neigh_vect,100);
ScoreEachBenchmark Score_hbond_sr_bb_("core.scoring.Score_100x_hbond_sr_bb",hbond_sr_bb,100);
ScoreEachBenchmark Score_hbond_lr_bb_("core.scoring.Score_100x_hbond_lr_bb",hbond_lr_bb,100);
ScoreEachBenchmark Score_hbond_bb_sc_("core.scoring.Score_100x_hbond_bb_sc",hbond_bb_sc,100);
//hbond_sr_bb_sc currently gives zero energies, though it takes a while to do so
ScoreEachBenchmark Score_hbond_sr_bb_sc_("core.scoring.Score_100x_hbond_sr_bb_sc",hbond_sr_bb_sc,100);
ScoreEachBenchmark Score_hbond_lr_bb_sc_("core.scoring.Score_100x_hbond_lr_bb_sc",hbond_lr_bb_sc,100);
ScoreEachBenchmark Score_hbond_sc_("core.scoring.Score_100x_hbond_sc",hbond_sc,100);
ScoreEachBenchmark Score_envsmooth_("core.scoring.Score_100x_envsmooth",envsmooth,100);
ScoreEachBenchmark Score_fa_cust_pair_dist_("core.scoring.Score_100x_fa_cust_pair_dist",fa_cust_pair_dist,100);
//fast centroid terms
ScoreEachBenchmark Score_rg_("core.scoring.Score_100x_rg",rg,100,"test_in_cen.pdb",true);
ScoreEachBenchmark Score_hs_pair_("core.scoring.Score_100x_hs_pair",hs_pair,100,"test_in_cen.pdb",true);
ScoreEachBenchmark Score_ss_pair_("core.scoring.Score_100x_ss_pair",ss_pair,100,"test_in_cen.pdb",true);
ScoreEachBenchmark Score_rsigma_("core.scoring.Score_100x_rsigma",rsigma,100,"test_in_cen.pdb",true);
ScoreEachBenchmark Score_sheet_("core.scoring.Score_100x_sheet",sheet,100,"test_in_cen.pdb",true);
ScoreEachBenchmark Score_vdw_("core.scoring.Score_100x_vdw",vdw,100,"test_in_cen.pdb",true);
ScoreEachBenchmark Score_cen_pair_smooth_("core.scoring.Score_100x_cen_pair_smooth",cen_pair_smooth,100,"test_in_cen.pdb",true);
ScoreEachBenchmark Score_cenpack_smooth_("core.scoring.Score_100x_cenpack_smooth",cenpack_smooth,100,"test_in_cen.pdb",true);
ScoreEachBenchmark Score_cenpack_("core.scoring.Score_100x_cenpack",cenpack,100,"test_in_cen.pdb",true);
ScoreEachBenchmark Score_pair_("core.scoring.Score_100x_pair",pair,100,"test_in_cen.pdb",true);

// average
ScoreEachBenchmark Score_fa_atr_("core.scoring.Score_10x_fa_atr",fa_atr,10);
ScoreEachBenchmark Score_fa_rep_("core.scoring.Score_10x_fa_rep",fa_rep,10);
ScoreEachBenchmark Score_fa_sol_("core.scoring.Score_10x_fa_sol",fa_sol,10);
ScoreEachBenchmark Score_fa_intra_atr_("core.scoring.Score_1000x_fa_intra_atr",fa_intra_atr,1000);
ScoreEachBenchmark Score_fa_intra_rep_("core.scoring.Score_1000x_fa_intra_rep",fa_intra_rep,1000);
ScoreEachBenchmark Score_fa_intra_sol_("core.scoring.Score_1000x_fa_intra_sol",fa_intra_sol,1000);
// temp until lk hack works with analytic etable ScoreEachBenchmark Score_lk_hack_("core.scoring.Score_100x_lk_hack",lk_hack,100);
// temp until lk ball works with analytic etable ScoreEachBenchmark Score_lk_ball_("core.scoring.Score_100x_lk_ball",lk_ball,100);
// temp until lk ball works with analytic etable ScoreEachBenchmark Score_lk_ball_iso_("core.scoring.Score_100x_lk_ball_iso)",lk_ball_iso,100);
ScoreEachBenchmark Score_mm_bend_("core.scoring.Score_200x_mm_bend",mm_bend,200);
ScoreEachBenchmark Score_mm_twist_("core.scoring.Score_200x_mm_twist",mm_twist,200);
ScoreEachBenchmark Score_lk_costheta_("core.scoring.Score_10x_lk_costheta",lk_costheta,10);
ScoreEachBenchmark Score_lk_polar_("core.scoring.Score_10x_lk_polar",lk_polar,10);
ScoreEachBenchmark Score_lk_nonpolar_("core.scoring.Score_10x_lk_nonpolar",lk_nonpolar,10);
//h2o_hbond currently gives zero energies, though it takes a while to do so
ScoreEachBenchmark Score_h2o_hbond_("core.scoring.Score_100x_h2o_hbond",h2o_hbond,100);
ScoreEachBenchmark Score_gauss_("core.scoring.Score_10x_gauss",gauss,10);
// Seg-faulting in debug mode ScoreEachBenchmark Score_geom_sol_("core.scoring.Score_100x_geom_sol",geom_sol,100);
ScoreEachBenchmark Score_occ_sol_fitted_("core.scoring.Score_10x_occ_sol_fitted",occ_sol_fitted,10);
ScoreEachBenchmark Score_occ_sol_fitted_onebody_("core.scoring.Score_5x_occ_sol_fitted_onebody",occ_sol_fitted_onebody,5);
ScoreEachBenchmark Score_fa_elec_("core.scoring.Score_10x_fa_elec",fa_elec,10);

// slow 100x
ScoreEachBenchmark Score_sa_("core.scoring.Score_1x_sa",sa,1);

// very slow 10x
ScoreEachBenchmark Score_gb_elec_("core.scoring.Score_1x_gb_elec",gb_elec,1);
ScoreEachBenchmark Score_pack_stat_("core.scoring.Score_1x_pack_stat",pack_stat,1);
//occ_sol_exact currently gives zero energies, though it takes quite a while to do so
ScoreEachBenchmark Score_occ_sol_exact_("core.scoring.Score_1x_occ_sol_exact",occ_sol_exact,1);

#endif // include guard
