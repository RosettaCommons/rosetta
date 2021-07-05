// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/metal_interface/RemoveMetalConnectionsMoverTests.cxxtest.hh
/// @brief  Unit tests for the RemoveMetalConnectionsMover, which removes the bonds
/// that were set up by the SetupMetalsMover or the -auto_setup_metals flag.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/metal_interface/RemoveMetalConnectionsMover.hh>
#include <protocols/simple_moves/SetupMetalsMover.hh>
#include <protocols/constraint_movers/ClearConstraintsMover.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("RemoveMetalConnectionsMoverTests");


class RemoveMetalConnectionsMoverTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init();

	}

	void tearDown(){

	}

	/// @brief Test that this mover does what it is supposed to.
	/// @details Sets up metals on pose then removes automatic setup.
	void test_mover(){
		// Load SOD1, a metalloprotein:
		core::pose::PoseOP pose( core::import_pose::pose_from_file("protocols/metal_interface/1spd_cleaned.pdb") );

		// Set up metals:
		{
			protocols::simple_moves::SetupMetalsMover set_up_metals;
			set_up_metals.apply( *pose );
		}

		// Check that we've got bonds:
		TR << "After first metal setup:" << std::endl;
		TR << "\tHis46:\t" << pose->residue_type(46).name() << std::endl;
		TR << "\tHis48:\t" << pose->residue_type(48).name() << std::endl;
		TR << "\tHis63:\t" << pose->residue_type(63).name() << std::endl;
		TR << "\tHis71:\t" << pose->residue_type(71).name() << std::endl;
		TR << "\tHis80:\t" << pose->residue_type(80).name() << std::endl;
		TR << "\tAsp83:\t" << pose->residue_type(83).name() << std::endl;
		TR << "\tHis120:\t" << pose->residue_type(120).name() << std::endl;
		TR << "\tCu154:\t" << pose->residue_type(154).name() << std::endl;
		TR << "\tZn155:\t" << pose->residue_type(155).name() << std::endl;
		TS_ASSERT_EQUALS( pose->residue_type(46).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(48).name(), "HIS:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(63).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(71).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(80).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(83).name(), "ASP:MP-OD2-connect:MP-OD2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(120).name(), "HIS:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(154).name(), "CU:MP-CU-metal_connect:MP-CU-metal_connect:MP-CU-metal_connect:MP-CU-metal_connect" );
		TS_ASSERT_EQUALS( pose->residue_type(155).name(), "ZN:MP-ZN-metal_connect:MP-ZN-metal_connect:MP-ZN-metal_connect:MP-ZN-metal_connect" );

		// Remove metals:
		{
			protocols::metal_interface::RemoveMetalConnectionsMover remove_metals;
			protocols::constraint_movers::ClearConstraintsMover clear_csts;
			remove_metals.apply(*pose);
			clear_csts.apply(*pose);
		}

		// Check that we've removed bonds:
		TR << "After first metal setup removal:" << std::endl;
		TR << "\tHis46:\t" << pose->residue_type(46).name() << std::endl;
		TR << "\tHis48:\t" << pose->residue_type(48).name() << std::endl;
		TR << "\tHis63:\t" << pose->residue_type(63).name() << std::endl;
		TR << "\tHis71:\t" << pose->residue_type(71).name() << std::endl;
		TR << "\tHis80:\t" << pose->residue_type(80).name() << std::endl;
		TR << "\tAsp83:\t" << pose->residue_type(83).name() << std::endl;
		TR << "\tHis120:\t" << pose->residue_type(120).name() << std::endl;
		TR << "\tCu154:\t" << pose->residue_type(154).name() << std::endl;
		TR << "\tZn155:\t" << pose->residue_type(155).name() << std::endl;
		TS_ASSERT_EQUALS( pose->residue_type(46).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(48).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(63).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(71).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(80).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(83).name(), "ASP" );
		TS_ASSERT_EQUALS( pose->residue_type(120).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(154).name(), "CU" );
		TS_ASSERT_EQUALS( pose->residue_type(155).name(), "ZN" );

		// Set up metals again:
		{
			protocols::simple_moves::SetupMetalsMover set_up_metals;
			set_up_metals.apply( *pose );
		}

		// Check that we've got bonds again:
		TR << "After second metal setup:" << std::endl;
		TR << "\tHis46:\t" << pose->residue_type(46).name() << std::endl;
		TR << "\tHis48:\t" << pose->residue_type(48).name() << std::endl;
		TR << "\tHis63:\t" << pose->residue_type(63).name() << std::endl;
		TR << "\tHis71:\t" << pose->residue_type(71).name() << std::endl;
		TR << "\tHis80:\t" << pose->residue_type(80).name() << std::endl;
		TR << "\tAsp83:\t" << pose->residue_type(83).name() << std::endl;
		TR << "\tHis120:\t" << pose->residue_type(120).name() << std::endl;
		TR << "\tCu154:\t" << pose->residue_type(154).name() << std::endl;
		TR << "\tZn155:\t" << pose->residue_type(155).name() << std::endl;
		TS_ASSERT_EQUALS( pose->residue_type(46).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(48).name(), "HIS:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(63).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(71).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(80).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(83).name(), "ASP:MP-OD2-connect:MP-OD2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(120).name(), "HIS:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(154).name(), "CU:MP-CU-metal_connect:MP-CU-metal_connect:MP-CU-metal_connect:MP-CU-metal_connect" );
		TS_ASSERT_EQUALS( pose->residue_type(155).name(), "ZN:MP-ZN-metal_connect:MP-ZN-metal_connect:MP-ZN-metal_connect:MP-ZN-metal_connect" );

		// Remove metals again:
		{
			protocols::metal_interface::RemoveMetalConnectionsMover remove_metals;
			protocols::constraint_movers::ClearConstraintsMover clear_csts;
			remove_metals.apply(*pose);
			clear_csts.apply(*pose);
		}

		// Check that we've removed bonds again:
		TR << "After second metal setup removal:" << std::endl;
		TR << "\tHis46:\t" << pose->residue_type(46).name() << std::endl;
		TR << "\tHis48:\t" << pose->residue_type(48).name() << std::endl;
		TR << "\tHis63:\t" << pose->residue_type(63).name() << std::endl;
		TR << "\tHis71:\t" << pose->residue_type(71).name() << std::endl;
		TR << "\tHis80:\t" << pose->residue_type(80).name() << std::endl;
		TR << "\tAsp83:\t" << pose->residue_type(83).name() << std::endl;
		TR << "\tHis120:\t" << pose->residue_type(120).name() << std::endl;
		TR << "\tCu154:\t" << pose->residue_type(154).name() << std::endl;
		TR << "\tZn155:\t" << pose->residue_type(155).name() << std::endl;
		TS_ASSERT_EQUALS( pose->residue_type(46).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(48).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(63).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(71).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(80).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(83).name(), "ASP" );
		TS_ASSERT_EQUALS( pose->residue_type(120).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(154).name(), "CU" );
		TS_ASSERT_EQUALS( pose->residue_type(155).name(), "ZN" );

	}

	/// @brief Test that this mover does what it is supposed to.
	/// @details Sets up metals on pose then removes automatic setup for everything except the
	/// HIS63 bonds.
	void test_mover_not_his63() {
		// Load SOD1, a metalloprotein:
		core::pose::PoseOP pose( core::import_pose::pose_from_file("protocols/metal_interface/1spd_cleaned.pdb") );

		// Set up a residue selector for the HIS63 and Cu atoms:
		core::select::residue_selector::ResidueIndexSelectorOP his_63_selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
		his_63_selector->append_index(63); //Select His63.
		core::select::residue_selector::NotResidueSelectorOP selector(
			utility::pointer::make_shared< core::select::residue_selector::NotResidueSelector >( his_63_selector )
		);

		// Set up metals:
		{
			protocols::simple_moves::SetupMetalsMover set_up_metals;
			set_up_metals.apply( *pose );
		}

		// Check that we've got bonds:
		TR << "After first metal setup:" << std::endl;
		TR << "\tHis46:\t" << pose->residue_type(46).name() << std::endl;
		TR << "\tHis48:\t" << pose->residue_type(48).name() << std::endl;
		TR << "\tHis63:\t" << pose->residue_type(63).name() << std::endl;
		TR << "\tHis71:\t" << pose->residue_type(71).name() << std::endl;
		TR << "\tHis80:\t" << pose->residue_type(80).name() << std::endl;
		TR << "\tAsp83:\t" << pose->residue_type(83).name() << std::endl;
		TR << "\tHis120:\t" << pose->residue_type(120).name() << std::endl;
		TR << "\tCu154:\t" << pose->residue_type(154).name() << std::endl;
		TR << "\tZn155:\t" << pose->residue_type(155).name() << std::endl;
		TS_ASSERT_EQUALS( pose->residue_type(46).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(48).name(), "HIS:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(63).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(71).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(80).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(83).name(), "ASP:MP-OD2-connect:MP-OD2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(120).name(), "HIS:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(154).name(), "CU:MP-CU-metal_connect:MP-CU-metal_connect:MP-CU-metal_connect:MP-CU-metal_connect" );
		TS_ASSERT_EQUALS( pose->residue_type(155).name(), "ZN:MP-ZN-metal_connect:MP-ZN-metal_connect:MP-ZN-metal_connect:MP-ZN-metal_connect" );

		// Remove metals:
		{
			protocols::metal_interface::RemoveMetalConnectionsMover remove_metals;
			remove_metals.set_residue_selector( selector );
			protocols::constraint_movers::ClearConstraintsMover clear_csts;
			remove_metals.apply(*pose);
			clear_csts.apply(*pose);
		}

		// Check that we've got the right bonds:
		TR << "After first metal setup removal:" << std::endl;
		TR << "\tHis46:\t" << pose->residue_type(46).name() << std::endl;
		TR << "\tHis48:\t" << pose->residue_type(48).name() << std::endl;
		TR << "\tHis63:\t" << pose->residue_type(63).name() << std::endl;
		TR << "\tHis71:\t" << pose->residue_type(71).name() << std::endl;
		TR << "\tHis80:\t" << pose->residue_type(80).name() << std::endl;
		TR << "\tAsp83:\t" << pose->residue_type(83).name() << std::endl;
		TR << "\tHis120:\t" << pose->residue_type(120).name() << std::endl;
		TR << "\tCu154:\t" << pose->residue_type(154).name() << std::endl;
		TR << "\tZn155:\t" << pose->residue_type(155).name() << std::endl;
		TS_ASSERT_EQUALS( pose->residue_type(46).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(48).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(63).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(71).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(80).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(83).name(), "ASP" );
		TS_ASSERT_EQUALS( pose->residue_type(120).name(), "HIS" );
		TS_ASSERT_EQUALS( pose->residue_type(154).name(), "CU:MP-CU-metal_connect" );
		TS_ASSERT_EQUALS( pose->residue_type(155).name(), "ZN:MP-ZN-metal_connect" );


	}

	/// @brief Test that this mover does what it is supposed to.
	/// @details Sets up metals on pose then removes automatic setup for only the
	/// HIS63-Cu bond.
	void test_mover_his63_cu_only() {
		// Load SOD1, a metalloprotein:
		core::pose::PoseOP pose( core::import_pose::pose_from_file("protocols/metal_interface/1spd_cleaned.pdb") );

		// Set up a residue selector for the HIS63 and Cu atoms:
		core::select::residue_selector::ResidueIndexSelectorOP selector( utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >() );
		selector->append_index(63); //Select His63.
		selector->append_index(154); //Select Cu154.

		// Set up metals:
		{
			protocols::simple_moves::SetupMetalsMover set_up_metals;
			set_up_metals.apply( *pose );
		}

		// Check that we've got bonds:
		TR << "After first metal setup:" << std::endl;
		TR << "\tHis46:\t" << pose->residue_type(46).name() << std::endl;
		TR << "\tHis48:\t" << pose->residue_type(48).name() << std::endl;
		TR << "\tHis63:\t" << pose->residue_type(63).name() << std::endl;
		TR << "\tHis71:\t" << pose->residue_type(71).name() << std::endl;
		TR << "\tHis80:\t" << pose->residue_type(80).name() << std::endl;
		TR << "\tAsp83:\t" << pose->residue_type(83).name() << std::endl;
		TR << "\tHis120:\t" << pose->residue_type(120).name() << std::endl;
		TR << "\tCu154:\t" << pose->residue_type(154).name() << std::endl;
		TR << "\tZn155:\t" << pose->residue_type(155).name() << std::endl;
		TS_ASSERT_EQUALS( pose->residue_type(46).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(48).name(), "HIS:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(63).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(71).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(80).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(83).name(), "ASP:MP-OD2-connect:MP-OD2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(120).name(), "HIS:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(154).name(), "CU:MP-CU-metal_connect:MP-CU-metal_connect:MP-CU-metal_connect:MP-CU-metal_connect" );
		TS_ASSERT_EQUALS( pose->residue_type(155).name(), "ZN:MP-ZN-metal_connect:MP-ZN-metal_connect:MP-ZN-metal_connect:MP-ZN-metal_connect" );

		// Remove metals:
		{
			protocols::metal_interface::RemoveMetalConnectionsMover remove_metals;
			remove_metals.set_residue_selector( selector );
			protocols::constraint_movers::ClearConstraintsMover clear_csts;
			remove_metals.apply(*pose);
			clear_csts.apply(*pose);
		}

		// Check that we've got the right bonds:
		TR << "After first metal setup removal:" << std::endl;
		TR << "\tHis46:\t" << pose->residue_type(46).name() << std::endl;
		TR << "\tHis48:\t" << pose->residue_type(48).name() << std::endl;
		TR << "\tHis63:\t" << pose->residue_type(63).name() << std::endl;
		TR << "\tHis71:\t" << pose->residue_type(71).name() << std::endl;
		TR << "\tHis80:\t" << pose->residue_type(80).name() << std::endl;
		TR << "\tAsp83:\t" << pose->residue_type(83).name() << std::endl;
		TR << "\tHis120:\t" << pose->residue_type(120).name() << std::endl;
		TR << "\tCu154:\t" << pose->residue_type(154).name() << std::endl;
		TR << "\tZn155:\t" << pose->residue_type(155).name() << std::endl;
		TS_ASSERT_EQUALS( pose->residue_type(46).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(48).name(), "HIS:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(63).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(71).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(80).name(), "HIS:MP-ND1-connect:MP-ND1-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(83).name(), "ASP:MP-OD2-connect:MP-OD2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(120).name(), "HIS:MP-NE2-connect:MP-NE2-pruneH" );
		TS_ASSERT_EQUALS( pose->residue_type(154).name(), "CU:MP-CU-metal_connect:MP-CU-metal_connect:MP-CU-metal_connect" );
		TS_ASSERT_EQUALS( pose->residue_type(155).name(), "ZN:MP-ZN-metal_connect:MP-ZN-metal_connect:MP-ZN-metal_connect:MP-ZN-metal_connect" );



	}

};
