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
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
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

		// Setup some test coordinates
		// outside - water exposed (above the membrane)
		p1_ = ahelical_pose_->conformation().residue( 23 ).atom( "CA" ).xyz();

		// in pore, water exposed
		p2_ = ahelical_pose_->conformation().residue( 125 ).atom( "CA" ).xyz();

		// interface exposed
		p3_ = ahelical_pose_->conformation().residue( 201 ).atom( "CA" ).xyz();

		// lipid exposed
		p4_ = ahelical_pose_->conformation().residue( 209 ).atom( "CA" ).xyz();

		// outside - water exposed (below the membrane)
		p5_ = ahelical_pose_->conformation().residue( 245 ).atom( "CA" ).xyz();

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
		TS_ASSERT_EQUALS( implicit_lipids->has_pore(), true );
		TS_ASSERT_EQUALS( implicit_lipids->is_helical(), true );

		// Check thickness and steepness values
		TS_ASSERT_EQUALS( implicit_lipids->water_thickness(), 15.351 );
		TS_ASSERT_EQUALS( implicit_lipids->water_steepness(), 0.343 );
		TS_ASSERT_EQUALS( implicit_lipids->water_pseudo_thickness(), 199.57 );

	}

	void test_f_depth(){

		using namespace core::conformation::membrane;

		TR << "Testing ImplicitLipidInfo function f_depth" << std::endl;

		ImplicitLipidInfoOP implicit_lipids( ahelical_pose_->conformation().membrane_info()->implicit_lipids() );

		TS_ASSERT_DELTA( implicit_lipids->f_depth( p1_.z() ), 0.9639, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_depth( p2_.z() ), 0.014, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_depth( p3_.z() ), 0.572,0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_depth( p4_.z() ), 0.018,0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_depth( p5_.z() ), 0.9868,0.001 );

	}

	void test_g_radius(){

		using namespace core::conformation::membrane;

		TR << "Testing ImplicitLipidInfo function g_radius" << std::endl;

		ImplicitLipidInfoOP implicit_lipids( ahelical_pose_->conformation().membrane_info()->implicit_lipids() );

		TS_ASSERT_DELTA( implicit_lipids->g_radius( p1_ ), 0.5487, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->g_radius( p2_ ), 0.1333, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->g_radius( p3_ ), 8.5096, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->g_radius( p4_ ), 280.8989,0.001 );
		TS_ASSERT_DELTA( implicit_lipids->g_radius( p5_ ), 1.7862,0.001 );

	}

	void test_f_cavity(){

		using namespace core::conformation::membrane;

		TR << "Testing ImplicitLipidInfo function f_cavity" << std::endl;

		ImplicitLipidInfoOP implicit_lipids( ahelical_pose_->conformation().membrane_info()->implicit_lipids() );

		TS_ASSERT_DELTA( implicit_lipids->f_cavity( p1_ ), 0.9975, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_cavity( p2_ ), 0.9999, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_cavity( p3_ ), 0.000, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_cavity( p4_ ), 0.00, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_cavity( p5_ ), 0.0030, 0.001 );

	}

	void test_f_hydration(){

		using namespace core::conformation::membrane;

		TR << "Testing ImplicitLipidInfo f_hydration function" << std::endl;

		ImplicitLipidInfoOP implicit_lipids( ahelical_pose_->conformation().membrane_info()->implicit_lipids() );

		TS_ASSERT_DELTA( implicit_lipids->f_hydration( p1_ ), 0.999, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_hydration( p2_ ), 0.999, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_hydration( p3_ ), 0.572, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_hydration( p4_ ), 0.0181, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_hydration( p5_ ), 0.986, 0.001 );

	}

	void test_f_depth_gradient(){

		using namespace core::conformation::membrane;

		TR << "Testing ImplicitLipidInfo f_depth_gradient function" << std::endl;

		ImplicitLipidInfoOP implicit_lipids( ahelical_pose_->conformation().membrane_info()->implicit_lipids() );

		TS_ASSERT_DELTA( implicit_lipids->f_depth_gradient( p1_.z() ), 0.0119, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_depth_gradient( p2_.z() ), 0.00489, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_depth_gradient( p3_.z() ), 0.0839,0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_depth_gradient( p4_.z() ), 0.00612,0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_depth_gradient( p5_.z() ), 0.0044,0.001 );

	}

	void test_df_cavity(){

		using namespace core::conformation::membrane;

		TR << "Testing ImplicitLipidInfo df_cavity function" << std::endl;

		ImplicitLipidInfoOP implicit_lipids( ahelical_pose_->conformation().membrane_info()->implicit_lipids() );

		TS_ASSERT( true );

		TS_ASSERT_DELTA( implicit_lipids->g_radius_gradient( p1_ ), -0.0265, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->g_radius_gradient( p2_ ), -0.0001, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->g_radius_gradient( p3_ ), 0, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->g_radius_gradient( p4_ ), 0.00 ,0.001 );
		TS_ASSERT_DELTA( implicit_lipids->g_radius_gradient( p5_ ), -0.0517,0.001 );

	}

	void test_f_cavity_gradient(){

		using namespace core::conformation::membrane;

		TR << "Testing ImplicitLipidInfo f_cavity_gradient function" << std::endl;

		ImplicitLipidInfoOP implicit_lipids( ahelical_pose_->conformation().membrane_info()->implicit_lipids() );

		TS_ASSERT_DELTA( implicit_lipids->f_cavity_gradient( implicit_lipids->g_radius( p1_ ) ), 0.0448, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_cavity_gradient( implicit_lipids->g_radius( p2_ ) ), 0.000, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_cavity_gradient( implicit_lipids->g_radius( p3_ ) ), 0, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_cavity_gradient( implicit_lipids->g_radius( p4_ ) ), 0.0, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_cavity_gradient( implicit_lipids->g_radius( p5_ ) ), 0.0168, 0.001 );
	}

	void test_f_hydration_gradient(){

		using namespace core::conformation::membrane;

		TR << "Testing ImplicitLipidInfo f_hydration_gradient function" << std::endl;

		ImplicitLipidInfoOP implicit_lipids( ahelical_pose_->conformation().membrane_info()->implicit_lipids() );

		TS_ASSERT_DELTA( implicit_lipids->f_hydration_gradient( p1_ ), -0.0009, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_hydration_gradient( p2_ ), -0.0001, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_hydration_gradient( p3_ ), 0.08398, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_hydration_gradient( p4_ ), 0.0061, 0.001 );
		TS_ASSERT_DELTA( implicit_lipids->f_hydration_gradient( p5_ ), 0.0043, 0.001 );

	}

private:

	// Starting pose
	core::pose::PoseOP ahelical_pose_;

	// Test coordinates
	numeric::xyzVector< core::Real > p1_;
	numeric::xyzVector< core::Real > p2_;
	numeric::xyzVector< core::Real > p3_;
	numeric::xyzVector< core::Real > p4_;
	numeric::xyzVector< core::Real > p5_;
	numeric::xyzVector< core::Real > p6_;


};
