// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/RamaMutationSelectorTests.cxxtest.hh
/// @brief  Unit tests for the RamaMutationSelector residue selector.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <test/util/pdb1ubq.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/cyclic_peptide/RamaMutationSelector.hh>


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/RamaPrePro.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

// Protocol Headers
#include <protocols/simple_moves/MutateResidue.hh>

// Basic Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("RamaMutationSelectorTests");


class RamaMutationSelectorTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "" );
		testpose_ = pdb1ubq5to13_poseop();
	}

	void tearDown(){

	}


	void do_select_by_xxx_attempt( std::string const &res_name, core::Real const &multiplier, core::Real const &threshold ) {
		core::pose::Pose pose( *(testpose_) );
		utility::vector1< core::Real > pro_energies( pose.total_residue(), 999999 );
		utility::vector1< core::Real > pro_energies2( pose.total_residue(), 999999 );

		core::scoring::ScoringManager const &score_man( *(core::scoring::ScoringManager::get_instance()) );
		core::scoring::RamaPrePro const &rama( score_man.get_RamaPrePro() );

		core::scoring::ScoreFunction sfxn;
		sfxn.set_weight( core::scoring::rama_prepro, multiplier );

		for ( core::Size ir(2), irmax(pose.total_residue()); ir<irmax; ++ir ) {
			core::pose::Pose pose2( pose );
			if ( res_name != "" ) {
				protocols::simple_moves::MutateResidue mutres( ir, res_name );
				mutres.apply(pose2);
			}
			(sfxn)(pose2);
			pro_energies[ir] = pose2.energies().residue_total_energies(ir)[ core::scoring::rama_prepro ] * multiplier;

			utility::vector1<core::Real> gradient; //Dummy var.
			utility::vector1<core::Real> mainchain_tors(2);
			mainchain_tors[1] = pose2.phi(ir);
			mainchain_tors[2] = pose2.psi(ir);
			rama.eval_rpp_rama_score( pose2.conformation(), pose2.residue_type_ptr(ir), pose2.residue_type_ptr(ir+1), mainchain_tors, pro_energies2[ir], gradient, false);
			pro_energies2[ir] *= multiplier;

			//DELETE THE FOLLOWING:
			//std::stringstream tempname;
			//tempname << "temp_" << res_name << "_" << ir << ".pdb";
			//pose2.dump_scored_pdb( tempname.str(), sfxn );
		}

		protocols::cyclic_peptide::RamaMutationSelector selector;
		selector.set_target_type( res_name );
		selector.set_score_threshold( threshold );
		selector.set_rama_prepro_multiplier( multiplier );

		core::select::residue_selector::ResidueSubset const selection( selector.apply( pose ) );
		core::select::residue_selector::ResidueSubset expected( pose.total_residue(), false );
		for ( core::Size i=2, imax=expected.size(); i<imax; ++i ) {
			expected[i] = pro_energies2[i] <= threshold;
		}

		for ( core::Size i(1), imax(expected.size()); i<=imax; ++i ) {
			TS_ASSERT( selection[i] == expected[i] );
		}

		TR << "\nRES\t2_BOD_EN\t1_BOD_EN\tSEL'D\tEXP'D\n";
		for ( core::Size i(1), imax(pro_energies.size()); i<=imax; ++i ) {
			TR << i << "\t" << pro_energies[i] << "\t" << pro_energies2[i] << (selection[i] ? "\tTRUE" : "\tFALSE") << (expected[i] ? "\tTRUE" : "\tFALSE") << "\n";
		}

		TR.flush();
	}

	void test_select_by_valine(){
		do_select_by_xxx_attempt( "VAL", 0.8, -0.1 );
	}

	void test_select_by_proline(){
		do_select_by_xxx_attempt( "PRO", 0.8, 3.8 );
	}

	void test_select_by_dproline(){
		do_select_by_xxx_attempt( "DPRO", 0.8, 9.81 );
	}

	void test_select_by_aib(){
		do_select_by_xxx_attempt( "AIB", 0.8, 1.5 );
	}

	void test_select_by_native(){
		do_select_by_xxx_attempt( "", 0.8, -0.3 );
	}

	core::pose::PoseCOP testpose_;

};



