// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/helical_bundle/HelicalBundlePDBInfoTests.cxxtest.hh
/// @brief  Unit tests confirming that helical bundle poses work properly with PDBInfo objects.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/Remarks.hh>

// Protocol Headers
#include <protocols/helical_bundle/MakeBundle.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("HelicalBundlePDBInfoTests");


class HelicalBundlePDBInfoTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();
		testpose_ = core::import_pose::pose_from_file( "protocols/helical_bundle/helical_bundle_input.pdb" );
	}

	void tearDown(){

	}



	void test_pdbinfo_on_append(){
		core::pose::PoseOP pose( testpose_->clone() );
		protocols::helical_bundle::MakeBundle mkbundle;

		mkbundle.set_default_crick_params_file("alpha_helix.crick_params");
		mkbundle.set_default_helix_length(20);
		utility::vector1<std::string> resnames(1);
		resnames[1] = "ALA";
		mkbundle.set_default_residue_name(resnames);
		mkbundle.set_reset_pose(false);
		mkbundle.add_helix();
		mkbundle.helix(1)->set_r0(5.0);
		mkbundle.helix(1)->set_omega0(0.01);
		mkbundle.add_helix();
		mkbundle.helix(2)->set_r0(7.0);
		mkbundle.helix(2)->set_omega0(-0.01);
		mkbundle.set_default_delta_omega0(3.141592654);
		mkbundle.apply(*pose);

		pose->pdb_info()->show( TR );
		core::io::Remarks const &remarks( pose->pdb_info()->remarks() );
		for ( core::Size i(0), imax(remarks.size()); i<imax; ++i ) {
			TR << remarks[i].num << "\t" << remarks[i].value << std::endl;
		}
		for ( core::Size ir(1), irmax(pose->total_residue()); ir<=irmax; ++ir ) {
			utility::vector1 < std::string > reslabels(pose->pdb_info()->get_reslabels(ir));
			TR << pose->residue(ir).name3() << ir << "\t";
			for ( core::Size j(1), jmax(reslabels.size()); j<=jmax; ++j ) {
				TR << reslabels[j] << std::endl;
			}
			if ( reslabels.size() == 0 ) TR << std::endl;
		}

		TS_ASSERT_EQUALS( pose->pdb_info()->get_reslabels(15)[1], "HBNet" );
		TS_ASSERT_EQUALS( pose->pdb_info()->get_reslabels(16)[1], "HBNet" );
		TS_ASSERT_EQUALS( pose->pdb_info()->get_reslabels(19)[1], "HBNet" );
		TS_ASSERT_EQUALS( pose->pdb_info()->get_reslabels(52)[1], "HBNet" );
		TS_ASSERT_EQUALS( pose->pdb_info()->get_reslabels(86)[1], "HBNet" );
		TS_ASSERT_EQUALS( pose->pdb_info()->get_reslabels(121)[1], "HBNet" );
		TS_ASSERT_EQUALS( pose->pdb_info()->get_reslabels(124)[1], "HBNet" );
		TS_ASSERT_EQUALS( pose->pdb_info()->get_reslabels(125)[1], "HBNet" );
	}

private:

	core::pose::PoseCOP testpose_;

};



