// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test/protocols/match/LigandConformer.cxxtest.hh
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


// Test headers
#include <cxxtest/TestSuite.h>

#include <protocols/match/Hit.hh>
#include <protocols/toolbox/match_enzdes_util/LigandConformer.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>

#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Utility headers
#include <utility/exit.hh>

/// Project headers
#include <core/types.hh>

// C++ headers
#include <string>
#include <iostream>
#include <sstream>

//Auto Headers
#include <core/id/AtomID.hh>
#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>


using namespace protocols::match;
//using namespace protocols::match::downstream;
using namespace protocols::toolbox::match_enzdes_util;


// --------------- Test Class --------------- //

class LigandConformerTests : public CxxTest::TestSuite {

public:

	typedef core::Real Real;
	typedef core::Size Size;

	// --------------- Fixtures --------------- //

	// Define a test fixture (some initial state that several tests share)
	// In CxxTest, setUp()/tearDown() are executed around each test case. If you need a fixture on the test
	// suite level, i.e. something that gets constructed once before all the tests in the test suite are run,
	// suites have to be dynamically created. See CxxTest sample directory for example.

	utility::vector1< Real > d_1, ang_U2D1_1, tor_U3D1_1, ang_U1D2_1, tor_U1D3_1, tor_U2D2_1;
	utility::vector1< Real > d_2, ang_U2D1_2, tor_U3D1_2, ang_U1D2_2, tor_U1D3_2, tor_U2D2_2;
	utility::vector1< Real > d_3, ang_U2D1_3, tor_U3D1_3, ang_U1D2_3, tor_U1D3_3, tor_U2D2_3;

	// Shared initialization goes here.
	void setUp() {
		using namespace core::chemical;

		core_init();

		ResidueTypeSet & restype_set(
			ChemicalManager::get_instance()->nonconst_residue_type_set( FA_STANDARD ));
		if ( ! restype_set.has_name( "CARBARYL_LG1" ) ) {
			utility::vector1< std::string > carbaryl_list;
			carbaryl_list.push_back( "protocols/match/carbaryl_LG1.params" );
			restype_set.read_files_for_custom_residue_types(carbaryl_list);
		}

		if ( ! restype_set.has_name( "MBH_LG1" ) ) {
			utility::vector1< std::string > carbaryl_list;
			carbaryl_list.push_back( "protocols/match/MBH_LG.params" );
			restype_set.read_files_for_custom_residue_types(carbaryl_list);
		}


		d_1.clear(); ang_U2D1_1.clear(); tor_U3D1_1.clear(); ang_U1D2_1.clear(); tor_U1D3_1.clear(); tor_U2D2_1.clear();
		d_2.clear(); ang_U2D1_2.clear(); tor_U3D1_2.clear(); ang_U1D2_2.clear(); tor_U1D3_2.clear(); tor_U2D2_2.clear();
		d_3.clear(); ang_U2D1_3.clear(); tor_U3D1_3.clear(); ang_U1D2_3.clear(); tor_U1D3_3.clear(); tor_U2D2_3.clear();

	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_ligand_conformer_recover_input_coords() {
		using namespace core;
		using namespace core::chemical;
		using namespace core::io::pdb;
		using namespace core::pose;

		Pose carbaryl_pose;
		core::import_pose::pose_from_pdb( carbaryl_pose,
			*(ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )),
			"protocols/match/E1cb_carbaryl_1his_oxy_1bb_10_2.pdb" );

		std::string const  at4( "C12" ),  at5( "C9" ),   at6( "C13" );
		std::string const oat1( "C17" ), oat2( "C19" ), oat3( "C22" ); // orientation atoms, oats

		Size const at4id( carbaryl_pose.residue(1).atom_index( at4 ) );
		Size const at5id( carbaryl_pose.residue(1).atom_index( at5 ) );
		Size const at6id( carbaryl_pose.residue(1).atom_index( at6 ) );

		Size const oat1id( carbaryl_pose.residue(1).atom_index( oat1 ) );
		Size const oat2id( carbaryl_pose.residue(1).atom_index( oat2 ) );
		Size const oat3id( carbaryl_pose.residue(1).atom_index( oat3 ) );

		LigandConformerOP ligconf( new LigandConformer );
		ligconf->initialize_from_residue(
			at4id, at5id, at6id, oat1id, oat2id, oat3id, carbaryl_pose.residue(1) );

		numeric::HomogeneousTransform< Real > oat3frame(
			carbaryl_pose.residue(1).xyz( oat1id ),
			carbaryl_pose.residue(1).xyz( oat2id ),
			carbaryl_pose.residue(1).xyz( oat3id ) );

		Vector euler = oat3frame.euler_angles_deg();

		Hit hit;
		hit.first()[ 1 ] = 1;
		hit.first()[ 2 ] = 1;
		hit.first()[ 3 ] = 1;
		hit.second()[ 1 ] = carbaryl_pose.residue(1).xyz( oat3id ).x();
		hit.second()[ 2 ] = carbaryl_pose.residue(1).xyz( oat3id ).y();
		hit.second()[ 3 ] = carbaryl_pose.residue(1).xyz( oat3id ).z();
		hit.second()[ 4 ] = euler( 1 );
		hit.second()[ 5 ] = euler( 2 );
		hit.second()[ 6 ] = euler( 3 );

		utility::vector1< core::id::AtomID > atomids( carbaryl_pose.residue(1).natoms() );
		for ( Size ii = 1; ii <= carbaryl_pose.residue(1).natoms(); ++ii ) atomids[ ii ] = core::id::AtomID( ii, 1 );
		utility::vector1< Vector > coords( carbaryl_pose.residue(1).natoms()  );

		ligconf->coordinates_from_orientation( hit.second(), atomids, coords );
		for ( Size ii = 1; ii <= carbaryl_pose.residue(1).natoms(); ++ii ) {
			TS_ASSERT( coords[ ii ].distance( carbaryl_pose.residue(1).xyz( ii ) ) < 1e-6 );
			//std::cout << coords[ ii ].x() << " " << coords[ ii ].y() << " " << coords[ ii ].z() << " vs ";
			//std::cout << carbaryl_pose.residue(1).atom_name( ii ) << " " << carbaryl_pose.residue(1).xyz( ii ).x() << " " << carbaryl_pose.residue(1).xyz( ii ).y() << " " << carbaryl_pose.residue(1).xyz( ii ).z() << std::endl;
		}

		//SingleDownstreamResidueWriter writer;
		//writer.set_restype( & (carbaryl_pose.residue(1).type()) );
		//writer.set_downstream_builder( rigid_builder );
		//writer.write_downstream_coordinates( hit, std::cout );

		utility::vector1< Vector > rotated_coords( carbaryl_pose.residue(1).natoms() );
		numeric::xyzMatrix< double > zrot = numeric::z_rotation_matrix_degrees( 30 );
		for ( Size ii = 1; ii <= carbaryl_pose.residue(1).natoms(); ++ii ) rotated_coords[ ii ] = zrot * carbaryl_pose.residue(1).xyz( ii );

		numeric::HomogeneousTransform< Real > oat3frame2(
			rotated_coords[ oat1id ],
			rotated_coords[ oat2id ],
			rotated_coords[ oat3id ] );

		Vector euler2 = oat3frame2.euler_angles_deg();

		Hit hit2;
		hit2.first()[ 1 ] = 1;
		hit2.first()[ 2 ] = 1;
		hit2.first()[ 3 ] = 1;
		hit2.second()[ 1 ] = rotated_coords[ oat3id ].x();
		hit2.second()[ 2 ] = rotated_coords[ oat3id ].y();
		hit2.second()[ 3 ] = rotated_coords[ oat3id ].z();
		hit2.second()[ 4 ] = euler2( 1 );
		hit2.second()[ 5 ] = euler2( 2 );
		hit2.second()[ 6 ] = euler2( 3 );

		ligconf->coordinates_from_orientation( hit2.second(), atomids, coords );
		for ( Size ii = 1; ii <= carbaryl_pose.residue(1).natoms(); ++ii ) {
			TS_ASSERT( coords[ ii ].distance( rotated_coords[ ii ] ) < 1e-6 );
			//std::cout << coords[ ii ].x() << " " << coords[ ii ].y() << " " << coords[ ii ].z() << " vs ";
			//std::cout << carbaryl_pose.residue(1).atom_name( ii ) << " " << rotated_coords[ ii ].x() << " " << rotated_coords[ ii ].y() << " " << rotated_coords[ ii ].z() << std::endl;
		}

	}

	void test_ligand_conformer_measure_global_coords() {
		using namespace core;
		using namespace core::chemical;
		using namespace core::io::pdb;
		using namespace core::pose;

		Pose carbaryl_pose;
		core::import_pose::pose_from_pdb( carbaryl_pose,
			*(ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )),
			"protocols/match/E1cb_carbaryl_1his_oxy_1bb_10_2.pdb" );

		std::string const dat1( "C12" ), dat2( "C9" ),  dat3( "C13" );
		std::string const oat1( "C17" ), oat2( "C19" ), oat3( "C22" ); // orientation atoms, oats

		Size const dat1id( carbaryl_pose.residue(1).atom_index( dat1 ) );
		Size const dat2id( carbaryl_pose.residue(1).atom_index( dat2 ) );
		Size const dat3id( carbaryl_pose.residue(1).atom_index( dat3 ) );

		Size const oat1id( carbaryl_pose.residue(1).atom_index( oat1 ) );
		Size const oat2id( carbaryl_pose.residue(1).atom_index( oat2 ) );
		Size const oat3id( carbaryl_pose.residue(1).atom_index( oat3 ) );

		LigandConformerOP ligconf( new LigandConformer );
		ligconf->initialize_from_residue(
			dat1id, dat2id, dat3id, oat1id, oat2id, oat3id, carbaryl_pose.residue(1) );

		numeric::HomogeneousTransform< Real > oat3frame(
			carbaryl_pose.residue(1).xyz( oat1id ),
			carbaryl_pose.residue(1).xyz( oat2id ),
			carbaryl_pose.residue(1).xyz( oat3id ) );

		numeric::HomogeneousTransform< Real > at3frame(
			carbaryl_pose.residue(1).xyz( dat1id ),
			carbaryl_pose.residue(1).xyz( dat2id ),
			carbaryl_pose.residue(1).xyz( dat3id ) );

		Vector euler = oat3frame.euler_angles_deg();

		for ( Size ii = 1; ii <= 3; ++ii ) if ( euler( ii ) < 0 ) euler( ii ) += 360.0;

		Hit hit;
		hit.first()[ 1 ] = 1;
		hit.first()[ 2 ] = 1;
		hit.first()[ 3 ] = 1;
		hit.second()[ 1 ] = carbaryl_pose.residue(1).xyz( oat3id ).x();
		hit.second()[ 2 ] = carbaryl_pose.residue(1).xyz( oat3id ).y();
		hit.second()[ 3 ] = carbaryl_pose.residue(1).xyz( oat3id ).z();
		hit.second()[ 4 ] = euler( 1 );
		hit.second()[ 5 ] = euler( 2 );
		hit.second()[ 6 ] = euler( 3 );

		Real6 global_coords = ligconf->global_orientation_from_frame3( at3frame );

		for ( Size ii = 1; ii <= 6; ++ii ) {
			TS_ASSERT_DELTA( global_coords[ ii ], hit.second()[ ii ], 1e-6 );
		}
	}

};


