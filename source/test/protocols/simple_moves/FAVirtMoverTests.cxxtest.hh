// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/simple_moves/FAVirtMoverTests.cxxtest.hh
/// @brief  Tests for FA to Virt mover and vice versa.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/simple_moves/ConvertRealToVirtualMover.hh>
#include <protocols/simple_moves/ConvertVirtualToRealMover.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>

// Protocol Headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>



static THREAD_LOCAL basic::Tracer TR("FAVirtMoverTests");


class FAVirtMoverTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){

		using namespace core::pose;
		using namespace core::import_pose;
		using namespace basic::options;

		core_init_with_additional_options( "-include_sugars");

		pose_from_file( pose_, "protocols/carbohydrates/N-linked_14-mer_glycan.pdb" , core::import_pose::PDB_file);

	}

	void test_fa_to_virt(){
		using namespace protocols::simple_moves;

		ConvertRealToVirtualMover fa_to_virt = ConvertRealToVirtualMover();

		//Convert whole pose to virtual.  Score should be zero.
		fa_to_virt.apply(pose_);

		for ( core::Size i = 1; i <= pose_.total_residue(); ++i ) {
			TS_ASSERT_EQUALS( pose_.residue( i ).is_virtual_residue(), true );
			for ( core::Size xx = 1; xx <= pose_.residue( i ).natoms(); ++xx ) {
				TS_ASSERT_EQUALS(pose_.residue( i ).is_virtual( xx ), true);
			}
		}

	}

	void test_virt_to_fa() {
		using namespace protocols::simple_moves;

		ConvertRealToVirtualMover fa_to_virt = ConvertRealToVirtualMover();
		ConvertVirtualToRealMover virt_to_fa = ConvertVirtualToRealMover();

		//Convert whole pose to virtual.  Score should be zero.
		fa_to_virt.apply( pose_ );

		//Convert back, make sure all residues and atoms are now real.
		virt_to_fa.apply( pose_ );

		for ( core::Size i = 1; i <= pose_.total_residue(); ++i ) {
			TS_ASSERT_EQUALS( pose_.residue( i ).is_virtual_residue(), false );

			// Can't be tested, as some virtual atoms exist for carbohydrates, etc.
			//for ( core::Size xx = 1; xx <= pose_.residue( i ).natoms(); ++xx ){
			// TS_ASSERT_EQUALS(pose_.residue( i ).is_virtual( xx ), false);
			//}
		}
	}

	void test_virt_scoring() {
		using namespace protocols::simple_moves;
		using namespace core::select::residue_selector;

		core::Real tol = 1e-6;

		ConvertRealToVirtualMover fa_to_virt = ConvertRealToVirtualMover();
		ConvertVirtualToRealMover virt_to_fa = ConvertVirtualToRealMover();

		core::scoring::ScoreFunctionOP scorefxn  = core::scoring::get_score_function();
		//fa_to_virt.apply( pose_ );

		core::Real score1 = scorefxn->score(pose_);


		utility::vector1< bool > subset( pose_.total_residue(), false);
		for ( core::Size i = 6; i <= pose_.total_residue(); ++i ) {
			subset[ i ] = true;
		}

		ReturnResidueSubsetSelectorOP subset_selector = ReturnResidueSubsetSelectorOP( new ReturnResidueSubsetSelector(subset) );

		//pose_.real_to_virtual(3);

		fa_to_virt.set_residue_selector(subset_selector);

		fa_to_virt.apply(pose_);

		core::Real score2 = scorefxn->score(pose_);

		virt_to_fa.apply(pose_);

		core::Real score3 = scorefxn->score(pose_);

		//std::cout << "Score3: " << score3 << std::endl;

		//std::cout << "Score1: " << score1 << " Score2: "<< score2 << std::endl;
		std::cout << "Score1: " << score1 << " Score2: "<< score2 << " Score3: "<< score3 << std::endl;

		TS_ASSERT( score1 != score2 );
		TS_ASSERT( score3 != score2 );
		TS_ASSERT( score1 - score3 < tol );

	}
	void tearDown(){

	}

private:

	core::pose::Pose pose_;







};



