// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/membrane/MembraneEnergyLandscapeSamplerTest.cxxtest.hh
/// @brief  This is an unit test for the MembraneEnergyLandscapeSampler mover
/// @author Rituparnasamanta (rituparna@utexas.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/AddMembraneMoverCreator.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane_benchmark/MembraneEnergyLandscapeSampler.hh>


static basic::Tracer TR("MembraneEnergyLandscapeSamplerTest");


class MembraneEnergyLandscapeSamplerTest : public CxxTest::TestSuite {
	//Define Variables
	protocols::membrane_benchmark::MembraneEnergyLandscapeSampler test_memb;
public:

	void setUp() {
		core_init();

	}

	void tearDown() {

	}
	/// @brief Calculate the axis and center of an alpha_helix
	void test_axis_center_of_alpha_helix(){

		using namespace core;
		using namespace protocols::membrane;
		using namespace core::import_pose;
		using namespace core::pose;

		// 1. TM domain of the M2 proton channel (single helix): 1mp6
		PoseOP m2_pose ( utility::pointer::make_shared< Pose >() );
		pose_from_file( *m2_pose, "protocols/membrane/1mp6_transformed.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb1 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1mp6.span" );
		add_memb1->apply( *m2_pose );

		//std::string filename5( "/Users/rsamant2/posetomembrane_1mp6.pdb" );
		//m2_pose->dump_pdb( filename5 );

		TR << "This is a unit test for: "<< test_memb.get_name()<< std::endl;
		Vector m2_axis( test_memb.getaxis(*m2_pose, 0) );
		Vector expected_axis_0( -0.124216,0.0455334,0.99121 );

		TR <<"The axis is:" << m2_axis.x() << "i+ " << m2_axis.y() << "j+" << m2_axis.z() << "k" << std::endl;
		//TR << "testing the weird abs function" << std::abs(m2_axis.z()) << std::endl;
		//the axis after add membrane should be close to z-axis
		TS_ASSERT_DELTA( std::abs(m2_axis.z()), 1.00, 0.01 );
		TS_ASSERT_DELTA( expected_axis_0, m2_axis, 0.001 );
		//TS_ASSERT(true)

		Vector m2_center( test_memb.getcenter(*m2_pose, 0) );
		Vector expected_center_0( 0.10425, -0.11875, 0.91825 );
		TR <<"The center is:" << m2_center.x() << "i+ " << m2_center.y() << "j+" << m2_center.z() << "k" << std::endl;

		TS_ASSERT_DELTA( expected_center_0, m2_center, 0.001 );
	}

	/// @brief Calculate the axis and center for a multipass protein
	void test_axis_center_of_protein(){

		using namespace core;
		using namespace protocols::membrane;
		using namespace core::import_pose;
		using namespace core::pose;

		// 1. TM domain of the 1U19__tr: Multipass protein
		PoseOP mp_pose ( utility::pointer::make_shared< Pose >() );
		pose_from_file( *mp_pose, "protocols/membrane/1U19__tr.pdb" , core::import_pose::PDB_file);
		AddMembraneMoverOP add_memb1 = utility::pointer::make_shared< AddMembraneMover >( "protocols/membrane/1U19__tr.span");
		add_memb1->apply( *mp_pose );


		Vector mp_axis( test_memb.getaxis(*mp_pose, 2) );
		Vector expected_axis_2( -0.060021, -0.0647234, 0.996097 );

		TR <<"The axis is:" << mp_axis.x() << "i+ " << mp_axis.y() << "j+" << mp_axis.z() << "k" << std::endl;

		//the axis after add membrane should be close to z-axis
		TS_ASSERT_DELTA( std::abs(mp_axis.z()), 1.00, 0.01 );
		TS_ASSERT_DELTA( expected_axis_2, mp_axis, 0.001 );

		Vector mp_center( test_memb.getcenter(*mp_pose, 2) );
		Vector expected_center_2( -0.11102, 0.278718, -0.56644 );
		TR <<"The center is:" << mp_center.x() << "i+ " << mp_center.y() << "j+" << mp_center.z() << "k" << std::endl;

		TS_ASSERT_DELTA( expected_center_2, mp_center, 0.001 );
		//TS_ASSERT(true)
	}




};
