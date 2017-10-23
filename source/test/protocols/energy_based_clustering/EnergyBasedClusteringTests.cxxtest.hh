// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/energy_based_clustering/EnergyBasedClusteringTests.cxxtest.hh
/// @brief  Unit tests for the EnergyBasedClusteringProtocol.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>

// Protocols Headers
#include <protocols/energy_based_clustering/EnergyBasedClusteringProtocol.hh>
#include <protocols/energy_based_clustering/EnergyBasedClusteringOptions.hh>

// Utility, etc Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("EnergyBasedClusteringTests");

class EnergyBasedClusteringTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp(){
		core_init_with_additional_options( "-energy_based_clustering:use_CB true" );

	}

	void tearDown(){

	}

	/// @brief Test the EnergyBasedClusteringProtocol::storeposedata() function.
	void test_storeposedata(){
		protocols::energy_based_clustering::EnergyBasedClusteringProtocol clusterer;
		core::pose::Pose pose;
		utility::vector1< core::Real > posedata;
		utility::vector1< numeric::xyzVector <core::Real> > alignmentdata;
		utility::vector1< core::Real > dihedral_mode_resconstruction_data;
		utility::vector1 < core::id::NamedAtomID > extra_atom_list;

		core::pose::make_pose_from_sequence( pose, "AX[OU3_ALA]X[B3A]", "fa_standard", true );
		for ( core::Size i(1); i<=pose.total_residue(); ++i ) {
			pose.set_omega(i, 180);
			pose.set_phi( i, -61 );
			pose.set_psi( i, -41 );
			if ( i > 1 ) pose.set_theta( i, 173.0);
			if ( i == 2 ) pose.set_mu( i, -64.1 );
		}
		pose.update_residue_neighbors();

		clusterer.storeposedata( pose, posedata, alignmentdata, dihedral_mode_resconstruction_data, protocols::energy_based_clustering::EBC_bb_cartesian, extra_atom_list );

		TS_ASSERT_EQUALS( alignmentdata.size(), 17 );

		TR << "\nAlignmentdata:\n";
		for ( core::Size i(1); i<=alignmentdata.size(); ++i ) {
			TR << alignmentdata[i].x() << "\t" << alignmentdata[i].y() << "\t" << alignmentdata[i].z() << "\n";
		}
		TR << std::endl;

		TS_ASSERT_DELTA( alignmentdata[1].x(), pose.xyz( core::id::NamedAtomID( "N", 1 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[1].y(), pose.xyz( core::id::NamedAtomID( "N", 1 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[1].z(), pose.xyz( core::id::NamedAtomID( "N", 1 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[2].x(), pose.xyz( core::id::NamedAtomID( "CA", 1 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[2].y(), pose.xyz( core::id::NamedAtomID( "CA", 1 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[2].z(), pose.xyz( core::id::NamedAtomID( "CA", 1 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[3].x(), pose.xyz( core::id::NamedAtomID( "C", 1 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[3].y(), pose.xyz( core::id::NamedAtomID( "C", 1 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[3].z(), pose.xyz( core::id::NamedAtomID( "C", 1 ) ).z(), 1e-6 );\
			TS_ASSERT_DELTA( alignmentdata[4].x(), pose.xyz( core::id::NamedAtomID( "O", 1 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[4].y(), pose.xyz( core::id::NamedAtomID( "O", 1 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[4].z(), pose.xyz( core::id::NamedAtomID( "O", 1 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[5].x(), pose.xyz( core::id::NamedAtomID( "CB", 1 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[5].y(), pose.xyz( core::id::NamedAtomID( "CB", 1 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[5].z(), pose.xyz( core::id::NamedAtomID( "CB", 1 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[6].x(), pose.xyz( core::id::NamedAtomID( "N", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[6].y(), pose.xyz( core::id::NamedAtomID( "N", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[6].z(), pose.xyz( core::id::NamedAtomID( "N", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[7].x(), pose.xyz( core::id::NamedAtomID( "CA", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[7].y(), pose.xyz( core::id::NamedAtomID( "CA", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[7].z(), pose.xyz( core::id::NamedAtomID( "CA", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[8].x(), pose.xyz( core::id::NamedAtomID( "CM", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[8].y(), pose.xyz( core::id::NamedAtomID( "CM", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[8].z(), pose.xyz( core::id::NamedAtomID( "CM", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[9].x(), pose.xyz( core::id::NamedAtomID( "NU", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[9].y(), pose.xyz( core::id::NamedAtomID( "NU", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[9].z(), pose.xyz( core::id::NamedAtomID( "NU", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[10].x(), pose.xyz( core::id::NamedAtomID( "C", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[10].y(), pose.xyz( core::id::NamedAtomID( "C", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[10].z(), pose.xyz( core::id::NamedAtomID( "C", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[11].x(), pose.xyz( core::id::NamedAtomID( "O", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[11].y(), pose.xyz( core::id::NamedAtomID( "O", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[11].z(), pose.xyz( core::id::NamedAtomID( "O", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[12].x(), pose.xyz( core::id::NamedAtomID( "CB", 2 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[12].y(), pose.xyz( core::id::NamedAtomID( "CB", 2 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[12].z(), pose.xyz( core::id::NamedAtomID( "CB", 2 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[13].x(), pose.xyz( core::id::NamedAtomID( "N", 3 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[13].y(), pose.xyz( core::id::NamedAtomID( "N", 3 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[13].z(), pose.xyz( core::id::NamedAtomID( "N", 3 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[14].x(), pose.xyz( core::id::NamedAtomID( "CA", 3 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[14].y(), pose.xyz( core::id::NamedAtomID( "CA", 3 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[14].z(), pose.xyz( core::id::NamedAtomID( "CA", 3 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[15].x(), pose.xyz( core::id::NamedAtomID( "CM", 3 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[15].y(), pose.xyz( core::id::NamedAtomID( "CM", 3 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[15].z(), pose.xyz( core::id::NamedAtomID( "CM", 3 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[16].x(), pose.xyz( core::id::NamedAtomID( "C", 3 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[16].y(), pose.xyz( core::id::NamedAtomID( "C", 3 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[16].z(), pose.xyz( core::id::NamedAtomID( "C", 3 ) ).z(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[17].x(), pose.xyz( core::id::NamedAtomID( "CB", 3 ) ).x(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[17].y(), pose.xyz( core::id::NamedAtomID( "CB", 3 ) ).y(), 1e-6 );
		TS_ASSERT_DELTA( alignmentdata[17].z(), pose.xyz( core::id::NamedAtomID( "CB", 3 ) ).z(), 1e-6 );
	}



};
