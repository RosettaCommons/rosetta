// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/conformation/membrane/ImplicitLipidInfo.cxxtest.hh
/// @brief  Unit test for ImplicitLipidInfo class
/// @author Rebecca Alford (rfalford12@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.hh>
#include <core/conformation/Conformation.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <numeric/xyzVector.hh>
#include <core/types.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static basic::Tracer TR("ImplicitLipidInfo");


class ImplicitLipidInfo : public CxxTest::TestSuite {

public:

	void setUp() {

		core_init_with_additional_options( "-mute all");

		using namespace core::import_pose;
		using namespace core::pose;
		using namespace protocols::membrane;

		// Load alpha helical protein test case
		ahelical_pose_ = utility::pointer::make_shared< Pose >();
		pose_from_file( *ahelical_pose_, "core/conformation/membrane/1U19_tr_ignorechain.pdb" , core::import_pose::PDB_file );
		AddMembraneMoverOP add_memb( new AddMembraneMover( "from_structure" ) );
		add_memb->apply( *ahelical_pose_ );

	}

	void tearDown(){

	}

	void test_implicit_lipid_info_initialization() {

		using namespace core::conformation::membrane;

		TR << "Testing lipid data stored in ImplicitLipidInfo" << std::endl;

		ImplicitLipidInfoOP implicit_lipids( ahelical_pose_->conformation().membrane_info()->implicit_lipids() );

		// Test data initialization for DLPC lipid type
		TS_ASSERT_EQUALS( implicit_lipids->chain_type(), "12:0/12:0" );
		TS_ASSERT_EQUALS( implicit_lipids->headgroup_type(), "PC" );
		TS_ASSERT_EQUALS( implicit_lipids->lipid_composition_name(), "DLPC" );
		TS_ASSERT_EQUALS( implicit_lipids->lipid_composition_name_long(), "1,2-dilauroyl-sn-glycero-3-phosphocholine" );
		TS_ASSERT_EQUALS( implicit_lipids->degrees_of_saturation(), 0 );
		TS_ASSERT_EQUALS( implicit_lipids->temperature(), 37 );
		TS_ASSERT_EQUALS( implicit_lipids->is_helical(), true );

		// Check thickness and steepness values
		TS_ASSERT_EQUALS( implicit_lipids->water_thickness(), 15.351 );
		TS_ASSERT_EQUALS( implicit_lipids->water_steepness(), 0.343 );
		TS_ASSERT_EQUALS( implicit_lipids->water_pseudo_thickness(), 199.57 );

	}



private:

	// Starting pose
	core::pose::PoseOP ahelical_pose_;

};
