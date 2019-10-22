// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/scoring/mhc_epitope_energy/MHCEpitopePredictorSVM.cxxtest.hh
/// @brief  Test suite for core::scoring::mhc_epitope_energy::MHCEpitopePredictorSVM, the MHC epitope predictor using a NMer SVM-based predictor.
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

// Test headers
#include <cxxtest/TestSuite.h>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorSVM.hh>
#include <core/scoring/nmer/NMerSVMEnergy.hh>

// Package Headers
#include <test/core/init_util.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <utility/pointer/memory.hh>

static basic::Tracer TR("core.scoring.mhc_epitope_energy.MHCEpitopePredictorSVM.cxxtest");

// --------------- Test Class --------------- //

// using declarations
using namespace core;
//using namespace core::pose;
using namespace core::scoring;
using namespace core::scoring::mhc_epitope_energy;

class MHCEpitopePredictorSVMTests : public CxxTest::TestSuite {

public:

	// --------------- Fixtures --------------- //
	// Shared variables
	std::string sequence;

	// Shared initialization goes here.
	void setUp() {
		core_init();

		// Setup test sequences
		sequence = "YFCTRAFRILAWIGI";
	}

	// Shared finalization goes here.
	void tearDown() {
	}

	/// @brief Check that the MHCEpitopePredictorSVM machinery is working.
	/// @author Brahm Yachnin
	void test_mhc_predictor_svm() {
		//////Create NMer objects to configure Predictors with: everything, no rank, no PSSM, no rank or PSSM.

		//Create an NMer object for use with the Predictor, with rank and PSSM in use
		core::scoring::methods::NMerSVMEnergyOP nmer_svm_rank_pssm( utility::pointer::make_shared<core::scoring::methods::NMerSVMEnergy>
			(
			9, //nmer_length
			false, //gate_svm_scores
			3, //term_length
			true, //use_pssm_features
			true, //avg_rank_as_energy
			0, //nmer_svm_scorecut
			utility::vector1<std::string>({
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10101_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10301_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10401_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10701_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10802_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB11101_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB11302_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB11501_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" )
			}), //svm_fname_vec
			utility::vector1<std::string>({
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB10101.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB10301.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB10401.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB10701.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB10802.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB11101.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB11302.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB11501.libsvm.test.out.sort.gz" )
			}), //svm_rank_fname_vec
			utility::vector1<std::string>({
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB10101_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB10301_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB10401_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB10701_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB10802_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB11101_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB11302_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB11501_nooverlap.9mer.norm.pssm" )
			}) //pssm_fname_vec
			)
		);
		//Create an NMer object for use with the Predictor, no rank and with PSSM
		core::scoring::methods::NMerSVMEnergyOP nmer_svm_pssm( utility::pointer::make_shared<core::scoring::methods::NMerSVMEnergy>
			(
			9, //nmer_length
			false, //gate_svm_scores
			3, //term_length
			true, //use_pssm_features
			false, //avg_rank_as_energy
			0, //nmer_svm_scorecut
			utility::vector1<std::string>({
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10101_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10301_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10401_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10701_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10802_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB11101_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB11302_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB11501_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" )
			}), //svm_fname_vec
			utility::vector1<std::string>(), //svm_rank_fname_vec
			utility::vector1<std::string>({
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB10101_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB10301_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB10401_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB10701_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB10802_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB11101_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB11302_nooverlap.9mer.norm.pssm" ),
			basic::database::full_name( "sequence/mhc_pssms/HLA-DRB11501_nooverlap.9mer.norm.pssm" )
			}) //pssm_fname_vec
			)
		);
		//Create an NMer object for use with the Predictor, with rank and with no PSSM
		core::scoring::methods::NMerSVMEnergyOP nmer_svm_rank( utility::pointer::make_shared<core::scoring::methods::NMerSVMEnergy>
			(
			9, //nmer_length
			false, //gate_svm_scores
			3, //term_length
			false, //use_pssm_features
			true, //avg_rank_as_energy
			0, //nmer_svm_scorecut
			utility::vector1<std::string>({
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10101_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10301_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10401_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10701_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10802_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB11101_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB11302_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB11501_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" )
			}), //svm_fname_vec
			utility::vector1<std::string>({
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB10101.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB10301.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB10401.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB10701.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB10802.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB11101.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB11302.libsvm.test.out.sort.gz" ),
			basic::database::full_name( "sequence/mhc_rank_svm_scores/HLA-DRB11501.libsvm.test.out.sort.gz" )
			}), //svm_rank_fname_vec
			utility::vector1<std::string>() //pssm_fname_vec
			)
		);
		//Create an NMer object for use with the Predictor, with no rank or PSSM
		core::scoring::methods::NMerSVMEnergyOP nmer_svm_only( utility::pointer::make_shared<core::scoring::methods::NMerSVMEnergy>
			(
			9, //nmer_length
			false, //gate_svm_scores
			3, //term_length
			false, //use_pssm_features
			false, //avg_rank_as_energy
			0, //nmer_svm_scorecut
			utility::vector1<std::string>({
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10101_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10301_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10401_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10701_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB10802_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB11101_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB11302_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" ),
			basic::database::full_name( "sequence/mhc_svms/HLA-DRB11501_nooverlap.libsvm.dat.noscale.nu0.5.min_mse.model" )
			}), //svm_fname_vec
			utility::vector1<std::string>(), //svm_rank_fname_vec
			utility::vector1<std::string>() //pssm_fname_vec
			)
		);

		//////Create MHCEpitopePredictorSVM objects corresponding to each NMer object.
		//Create a MHCEpitopePredictorSVM object, using the rank-pssm NMer object declared above.
		MHCEpitopePredictorSVMOP pred_svm_rank_pssm( utility::pointer::make_shared<MHCEpitopePredictorSVM>(nmer_svm_rank_pssm) );
		//Create a MHCEpitopePredictorSVM object, using the rank NMer object declared above.
		MHCEpitopePredictorSVMOP pred_svm_rank( utility::pointer::make_shared<MHCEpitopePredictorSVM>(nmer_svm_rank) );
		//Create a MHCEpitopePredictorSVM object, using the pssm NMer object declared above.
		MHCEpitopePredictorSVMOP pred_svm_pssm( utility::pointer::make_shared<MHCEpitopePredictorSVM>(nmer_svm_pssm) );
		//Create a MHCEpitopePredictorSVM object, using the svm only NMer object declared above.
		MHCEpitopePredictorSVMOP pred_svm_only( utility::pointer::make_shared<MHCEpitopePredictorSVM>(nmer_svm_only) );

		//////Tests for the pred_svm_rank_pssm/nmer_svm_rank_pssm set

		//Check that the NMer object is the same as the one held by the Predictor
		TS_ASSERT_EQUALS( nmer_svm_rank_pssm, pred_svm_rank_pssm->get_svm() );

		//Check the score of the known sequence
		core::Real score_pred_rank_pssm = pred_svm_rank_pssm->score(sequence);
		TS_ASSERT_DELTA(score_pred_rank_pssm, 0.57222, 0.00001);

		//Check that we get the same value when predicting directly from the NMer object
		core::Real energy_nmer = 0;
		core::Real rank_nmer = 0;
		utility::vector1<core::Real>svm_energies(nmer_svm_rank_pssm->n_svms(), core::Real(0.));
		utility::vector1<core::Real>svm_ranks(nmer_svm_rank_pssm->n_svms(), core::Real(0.));
		nmer_svm_rank_pssm->get_residue_energy_from_sequence(
			sequence,
			1 + nmer_svm_rank_pssm->term_length(), //p1 position, excluding the termini
			energy_nmer, //Input energy average
			rank_nmer, //Input rank average
			svm_energies, //Vector of energies, per svm
			svm_ranks //Vector of ranks, per svm
		);
		TS_ASSERT_EQUALS(score_pred_rank_pssm, rank_nmer);

		//////Tests for the pred_svm_rank/nmer_svm_rank set
		//Check that the NMer object is the same as the one held by the Predictor
		TS_ASSERT_EQUALS( nmer_svm_rank, pred_svm_rank->get_svm() );

		//Check the score of the known sequence
		core::Real score_pred_rank = pred_svm_rank->score(sequence);
		TS_ASSERT_DELTA(score_pred_rank, 0.61451, 0.00001);

		//Check that we get the same value when predicting directly from the NMer object
		energy_nmer = 0;
		rank_nmer = 0;
		svm_energies.assign(nmer_svm_rank->n_svms(), core::Real(0.));
		svm_ranks.assign(nmer_svm_rank->n_svms(), core::Real(0.));
		nmer_svm_rank->get_residue_energy_from_sequence(
			sequence,
			1 + nmer_svm_rank->term_length(), //p1 position, excluding the termini
			energy_nmer, //Input energy average
			rank_nmer, //Input rank average
			svm_energies, //Vector of energies, per svm
			svm_ranks //Vector of ranks, per svm
		);
		TS_ASSERT_EQUALS(score_pred_rank, rank_nmer);

		//////Tests for the pred_svm_pssm/nmer_svm_pssm set
		//Check that the NMer object is the same as the one held by the Predictor
		TS_ASSERT_EQUALS( nmer_svm_pssm, pred_svm_pssm->get_svm() );

		//Check the score of the known sequence
		core::Real score_pred_pssm = pred_svm_pssm->score(sequence);
		TS_ASSERT_DELTA(score_pred_pssm, 0.20510, 0.00001);

		//Check that we get the same value when predicting directly from the NMer object
		energy_nmer = 0;
		rank_nmer = 0;
		svm_energies.assign(nmer_svm_pssm->n_svms(), core::Real(0.));
		svm_ranks.assign(nmer_svm_pssm->n_svms(), core::Real(0.));
		nmer_svm_pssm->get_residue_energy_from_sequence(
			sequence,
			1 + nmer_svm_pssm->term_length(), //p1 position, excluding the termini
			energy_nmer, //Input energy average
			rank_nmer, //Input rank average
			svm_energies, //Vector of energies, per svm
			svm_ranks //Vector of ranks, per svm
		);
		TS_ASSERT_EQUALS(score_pred_pssm, energy_nmer);

		//////Tests for the pred_svm_only/nmer_svm_only set
		//Check that the NMer object is the same as the one held by the Predictor
		TS_ASSERT_EQUALS( nmer_svm_only, pred_svm_only->get_svm() );

		//Check the score of the known sequence
		core::Real score_pred_svm_only = pred_svm_only->score(sequence);
		TS_ASSERT_DELTA(score_pred_svm_only, 0.22395, 0.00001);

		//Check that we get the same value when predicting directly from the NMer object
		energy_nmer = 0;
		rank_nmer = 0;
		svm_energies.assign(nmer_svm_only->n_svms(), core::Real(0.));
		svm_ranks.assign(nmer_svm_only->n_svms(), core::Real(0.));
		nmer_svm_only->get_residue_energy_from_sequence(
			sequence,
			1 + nmer_svm_only->term_length(), //p1 position, excluding the termini
			energy_nmer, //Input energy average
			rank_nmer, //Input rank average
			svm_energies, //Vector of energies, per svm
			svm_ranks //Vector of ranks, per svm
		);
		TS_ASSERT_EQUALS(score_pred_svm_only, energy_nmer);

		TR << "End of test_mhc_predictor_svm." << std::endl;
	}

};
