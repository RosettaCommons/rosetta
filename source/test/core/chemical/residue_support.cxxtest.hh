// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/chemical/residue_support.cxxtest.hh
/// @brief unit tests for the residue_support file
/// @author Rocco Moretti (rmorettiase@gmail.com)


// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/residue_support.hh>
#include <core/chemical/Atom.hh>

// Project Headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/MMAtomType.hh>
#include <core/chemical/residue_io.hh>

#include <core/types.hh>

// Platform Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>
#include <set>
#include <cmath>

static basic::Tracer TR("core.chemical.restype_support.cxxtest");

class residue_support_Tests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	void test_renaming() {
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		ResidueTypeSetCOP rsd_types = cm->residue_type_set(tag);

		ResidueTypeCOP rsd_ref( rsd_types->name_map("LYS").get_self_ptr() );
		ResidueTypeOP rsd;
		std::set< std::string > names;

		// Already named - should be no changes.
		rsd = ResidueTypeOP( new ResidueType( *rsd_ref ) );
		rename_atoms( *rsd, /*preserve=*/true);
		for ( core::Size ii(1); ii <= rsd->natoms(); ++ii ) {
			TS_ASSERT_EQUALS( rsd->atom(ii).name(), rsd_ref->atom(ii).name() );
		}

		// Force renaming.
		rsd = ResidueTypeOP( new ResidueType( *rsd_ref ) );
		names.clear();
		rename_atoms( *rsd, /*preserve=*/false);
		for ( core::Size ii(1); ii <= rsd->natoms(); ++ii ) {
			TS_ASSERT_EQUALS( names.count( rsd->atom(ii).name() ), 0 ); // Names should be unique
			names.insert( rsd->atom(ii).name() );
		}
		//non-terminal Lysine: 6 carbons, 1 oxygen, 2 nitrogen, 13 hydrogens.
		TS_ASSERT_EQUALS( names.size(), 22 );
		TS_ASSERT_EQUALS( names.count( " O1 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " O2 " ), 0 );
		TS_ASSERT_EQUALS( names.count( " N1 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " N2 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " C1 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " C6 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " C7 " ), 0 );
		TS_ASSERT_EQUALS( names.count( " H5 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " H13" ), 1 );
		TS_ASSERT_EQUALS( names.count( " H14" ), 0 );

		// Partial renaming.
		rsd = ResidueTypeOP( new ResidueType( *rsd_ref ) );
		rsd->atom(" NZ ").name(" N  ");
		rsd->atom("1HB ").name("");
		rsd->atom("2HB ").name("");
		rsd->atom(" HA ").name(" H  ");
		rename_atoms( *rsd, /*preserve=*/true);
		names.clear();
		for ( core::Size ii(1); ii <= rsd->natoms(); ++ii ) {
			TS_ASSERT_EQUALS( names.count( rsd->atom(ii).name() ), 0 ); // Names should be unique
			names.insert( rsd->atom(ii).name() );
		}
		TS_ASSERT_EQUALS( names.size(), 22 );
		TS_ASSERT_EQUALS( names.count( " O  " ), 1 );
		TS_ASSERT_EQUALS( names.count( " CA " ), 1 );
		TS_ASSERT_EQUALS( names.count( " N  " ), 0 );
		TS_ASSERT_EQUALS( names.count( " NZ " ), 0 );
		TS_ASSERT_EQUALS( names.count( " N1 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " N2 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " H1 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " H2 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " H3 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " H4 " ), 1 );
		TS_ASSERT_EQUALS( names.count( " H5 " ), 0 );
		TS_ASSERT_EQUALS( names.count( " H  " ), 0 );
		TS_ASSERT_EQUALS( names.count( " HA " ), 0 );
		TS_ASSERT_EQUALS( names.count( ""     ), 0 );
		TS_ASSERT_EQUALS( names.count( "1HZ " ), 1 );

	}

	void test_rigid_matrix() {
		using namespace core;
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		// minirosett_database/chemical/atom_type_sets/<tag>
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		// minirosetta_database/chemical/element_sets/<tag>
		ElementSetCOP element_types = cm->element_set("default");
		// minirosetta_database/chemical/mm_atom_type_sets/<tag>
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);

		ResidueTypeOP rsd( new ResidueType( atom_types, element_types, mm_atom_types, NULL ) );

		rsd->add_atom( "C1", "aroC", "VIRT", 0 );
		rsd->add_atom( "C2", "aroC", "VIRT", 0 );
		rsd->add_atom( "C3", "aroC", "VIRT", 0 );
		rsd->add_atom( "C4", "aroC", "VIRT", 0 );
		rsd->add_atom( "C5", "aroC", "VIRT", 0 );
		rsd->add_atom( "C6", "aroC", "VIRT", 0 );
		rsd->add_atom( "C7", "aroC", "VIRT", 0 );
		rsd->add_atom( "C8", "aroC", "VIRT", 0 );
		rsd->add_atom( "C9", "aroC", "VIRT", 0 );
		rsd->add_atom( "H10", "Haro", "VIRT", 0 );
		rsd->add_atom( "H11", "Haro", "VIRT", 0 );
		rsd->add_atom( "H12", "Haro", "VIRT", 0 );
		rsd->add_atom( "O13", "OOC", "VIRT", 0 );
		// shape :
		//    8-9-10
		//    #
		// 11 7=6
		//  |   |
		//  2-4-5-13
		//  | |   |
		//  1-3   12
		rsd->add_bond( "C1", "C2" );
		rsd->add_bond( "C1", "C3" );
		rsd->add_bond( "C3", "C4" );
		rsd->add_bond( "C4", "C2" );
		rsd->add_bond( "C4", "C5" );
		rsd->add_bond( "C5", "C6" );
		rsd->add_bond( "C6", "C7", DoubleBond );
		rsd->add_bond( "C7", "C8", TripleBond );
		rsd->add_bond( "C8", "C9" );
		rsd->add_bond( "C9", "H10" );
		rsd->add_bond( "C2", "H11" );
		rsd->add_bond( "C5", "O13" );
		rsd->add_bond( "O13", "H12" );

		rsd->atom("C1").ideal_xyz( Vector(0,0,0) );
		rsd->atom("C2").ideal_xyz( Vector(0,1,0) );
		rsd->atom("C3").ideal_xyz( Vector(1,0,0) );
		rsd->atom("C4").ideal_xyz( Vector(1,1,0) );
		rsd->atom("C5").ideal_xyz( Vector(2.1,1,0) );
		rsd->atom("C6").ideal_xyz( Vector(2.1,2.1,0) );
		rsd->atom("C7").ideal_xyz( Vector(1,2.1,0) );
		rsd->atom("C8").ideal_xyz( Vector(1,3.3,0) );
		rsd->atom("C9").ideal_xyz( Vector(2.5,3.3,0) );
		rsd->atom("H10").ideal_xyz( Vector(3.0,3.3,0) );
		rsd->atom("H11").ideal_xyz( Vector(0,1.5,0) );
		rsd->atom("H12").ideal_xyz( Vector(3.3,0,0) );
		rsd->atom("O13").ideal_xyz( Vector(3.3,1,0) );

		rsd->bond("C1","C2").ringness(BondInRing);
		rsd->bond("C1","C3").ringness(BondInRing);
		rsd->bond("C3","C4").ringness(BondInRing);
		rsd->bond("C2","C4").ringness(BondInRing);

		core::Size natoms( rsd->natoms() );
		utility::vector1< utility::vector1< core::Real > > distances(natoms, utility::vector1< core::Real >( natoms, 1e9 ) );

		calculate_rigid_matrix( *rsd, distances );

		core::Real const delta( 0.001 );
		// Distances should be symetrical.
		for ( core::Size ii(1); ii <= natoms; ++ii ) {
			for ( core::Size jj(ii+1); jj <= natoms; ++jj ) {
				TS_ASSERT_DELTA( distances[ii][jj], distances[jj][ii], delta );
			}
		}
		//Atom ordering should be the same as numerical ordering.
		TS_ASSERT_DELTA( distances[1][2], 1.0, delta );
		TS_ASSERT_DELTA( distances[2][1], 1.0, delta );
		TS_ASSERT_DELTA( distances[1][4], sqrt( 1.0*1.0 + 1.0*1.0), delta );
		TS_ASSERT_DELTA( distances[1][4], distances[4][1], delta );
		TS_ASSERT_DELTA( distances[1][5], sqrt(2) + 1.1, delta );
		TS_ASSERT_DELTA( distances[5][3], 1 + 1.1, delta );
		TS_ASSERT_DELTA( distances[1][6], sqrt(2) + 1.1 + 1.1, delta );
		TS_ASSERT_DELTA( distances[6][4], 1.1 + 1.1, delta );
		//Bonding issues
		TS_ASSERT_DELTA( distances[7][8],1.2, delta );
		TS_ASSERT_DELTA( distances[7][6],1.1, delta );
		TS_ASSERT_DELTA( distances[6][8], sqrt( 1.1*1.1 + 1.2*1.2 ), delta );
		TS_ASSERT_DELTA( distances[4][8], 1.1 + 1.1 + sqrt( 1.1*1.1 + 1.2*1.2 ), delta );
		//C9 is a Non-rotameric stub, so should be non-rotameric
		TS_ASSERT_DELTA( distances[9][6],sqrt( 0.4*0.4 + 1.2*1.2), delta );
		//Across everything
		TS_ASSERT_DELTA( distances[1][9], sqrt(2) + 1.1 + 1.1 + sqrt( 0.4*0.4 + 1.2*1.2 ), delta );
		//Hydrogens add to rigid unit.
		TS_ASSERT_DELTA( distances[7][10],sqrt(2*2 + 1.2*1.2), delta );
		TS_ASSERT_DELTA( distances[3][11],sqrt(1*1+1.5*1.5), delta );
		//But polar hydrogens don't
		TS_ASSERT_DELTA( distances[5][12],1.2+1, delta );
	}


	void test_possible_nbr() {
		using namespace core;
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		// minirosett_database/chemical/atom_type_sets/<tag>
		AtomTypeSetCOP atom_types = cm->atom_type_set(tag);
		// minirosetta_database/chemical/element_sets/<tag>
		ElementSetCOP element_types = cm->element_set("default");
		// minirosetta_database/chemical/mm_atom_type_sets/<tag>
		MMAtomTypeSetCOP mm_atom_types = cm->mm_atom_type_set(tag);

		ResidueTypeOP rsd( new ResidueType( atom_types, element_types, mm_atom_types, NULL ) );

		rsd->add_atom( "C1", "aroC", "VIRT", 0 );
		rsd->add_atom( "C2", "aroC", "VIRT", 0 );
		rsd->add_atom( "C3", "aroC", "VIRT", 0 );
		rsd->add_atom( "C4", "aroC", "VIRT", 0 );
		rsd->add_atom( "C5", "aroC", "VIRT", 0 );
		rsd->add_atom( "C6", "aroC", "VIRT", 0 );
		// "E" shape, and all rigid
		rsd->add_bond( "C1", "C2", DoubleBond );
		rsd->add_bond( "C2", "C3", DoubleBond );
		rsd->add_bond( "C3", "C4", DoubleBond );
		rsd->add_bond( "C4", "C5", DoubleBond );
		rsd->add_bond( "C3", "C6", DoubleBond );

		rsd->atom("C1").ideal_xyz( Vector(0,0,0) );
		rsd->atom("C2").ideal_xyz( Vector(2,2,0) );
		rsd->atom("C3").ideal_xyz( Vector(3,1,0) );
		rsd->atom("C4").ideal_xyz( Vector(4,0,0) );
		rsd->atom("C5").ideal_xyz( Vector(2,-2,0) );
		rsd->atom("C6").ideal_xyz( Vector(2,0,0) );

		core::Real maxdist(0);
		VD nbr;

		//  // C6 is omitted because it's singly bonded.
		//  nbr = ResidueType::null_vertex;
		//  maxdist = find_nbr_dist(*rsd, nbr);
		//  TS_ASSERT_EQUALS( nbr, rsd->vd_from_name("C3") );
		//  TS_ASSERT_DELTA( maxdist, 3 * sqrt( 2 ), 1e-6 ); //To C1 & C5
		//
		//  // Doubly bonded H is still omitted.
		//  rsd->add_bond( "C6", "C2" );
		//  rsd->atom("C6").element_type( element_types->element("H") );
		//  nbr = ResidueType::null_vertex;
		//  maxdist = find_nbr_dist(*rsd, nbr);
		//  TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C3");
		//  TS_ASSERT_DELTA( maxdist, 3 * sqrt(2) , 1e-6 ); //To C1 & C5
		//
		//  // But a ring through an H will count
		//  rsd->add_bond( "C6", "C2" );
		//  nbr = ResidueType::null_vertex;
		//  maxdist = find_nbr_dist(*rsd, nbr);
		//  TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C6");
		//  TS_ASSERT_DELTA( maxdist, 2.0 , 1e-6 );
		//
		//
		//  rsd->atom("C4").ideal_xyz( Vector(3,2,0) );
		//  rsd->atom("C5").ideal_xyz( Vector(1,3,0) );
		//  nbr = ResidueType::null_vertex;
		//  maxdist = find_nbr_dist(*rsd, nbr);
		//  TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C2" );
		//  TS_ASSERT_DELTA( maxdist, sqrt( 2.0*2.0+2.0*2.0 ), 1e-6 ); //To C1
		//
		//  ////Hydrogens aren't included in furthest distance calculation,
		//  //// but they still count as bonded groups.
		//  //rsd->atom("C1").element_type( element_types->element("H") );
		//  //nbr = ResidueType::null_vertex;
		//  //maxdist = find_nbr_dist(*rsd, nbr);
		//  //TS_ASSERT_EQUALS( nbr, rsd->vd_from_name("C2") );
		//  //TS_ASSERT_DELTA( maxdist, sqrt( 2.0 ), 1e-6 ); //To C3 & C5 - As hydrogen, C6 = sqrt(8) doesn't count.

		// C6 is omitted because it's singly bonded.
		nbr = ResidueType::null_vertex;
		maxdist = find_nbr_dist(*rsd, nbr);
		TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C3" );
		TS_ASSERT_DELTA( maxdist, sqrt( 1.0*1.0 + 3.0*3.0 ), 1e-6 ); //To C1 & C5

		rsd->add_bond( "C6", "C2" );
		nbr = ResidueType::null_vertex;
		maxdist = find_nbr_dist(*rsd, nbr);
		TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C6" );
		TS_ASSERT_DELTA( maxdist, 2.0 , 1e-6 );

		// Doubly bonded H is still omitted.
		rsd->atom("C6").element_type( element_types->element("H") );
		nbr = ResidueType::null_vertex;
		maxdist = find_nbr_dist(*rsd, nbr);
		TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C3" );
		TS_ASSERT_DELTA( maxdist, sqrt( 1.0*1.0 + 3.0*3.0 ) , 1e-6 ); //To C1 & C5

		rsd->atom("C4").ideal_xyz( Vector(3,2,0) );
		rsd->atom("C5").ideal_xyz( Vector(1,3,0) );
		nbr = ResidueType::null_vertex;
		maxdist = find_nbr_dist(*rsd, nbr);
		TS_ASSERT_EQUALS( rsd->atom_name(nbr), "C2" );
		TS_ASSERT_DELTA( maxdist, sqrt( 2.0*2.0+2.0*2.0 ), 1e-6 ); //To C1
	}

	void test_recharging() {
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		std::string const tag(FA_STANDARD);
		ResidueTypeSetCOP rsd_types = cm->residue_type_set(tag);

		ResidueTypeOP rsd( new ResidueType( rsd_types->name_map("TYR") ) );

		TR << "Testing standard charging." << std::endl;

		for ( core::Size ii(1); ii <= rsd->natoms(); ++ii ) {
			rsd->atom(ii).formal_charge( 0.0 );
			rsd->atom(ii).charge( -1.234 ); // Will be reset
		}

		rosetta_recharge_fullatom( *rsd );

		// 21 atoms in TYR -- Nbb + Cabb + CObb + OCbb + CH2 + aroC * 6 + OH + HNbb + Hapo * 3 + Haro * 4 + Hpol
		core::Real Y_naive_charge = -0.47 + 0.07 + 0.51 + -0.51 + -0.18 + -0.115*6 + -0.66 + 0.31 + 0.115*4 + 0.095*3 + 0.43;
		TS_ASSERT_DELTA( rsd->atom(" N  ").charge(), -0.470 - Y_naive_charge/21, 1e-4 ); //"Nbb "
		TS_ASSERT_DELTA( rsd->atom(" CA ").charge(),  0.070 - Y_naive_charge/21, 1e-4 ); //"CAbb"
		TS_ASSERT_DELTA( rsd->atom(" OH ").charge(), -0.660 - Y_naive_charge/21, 1e-4 ); //"OH  "

		core::Real net_charge(0);
		for ( core::Size ii(1); ii <= rsd->natoms(); ++ii ) {
			net_charge += rsd->atom(ii).charge();
		}
		TS_ASSERT_DELTA( net_charge, 0, 1e-4 );

		TR << "Testing charging with formal charges" << std::endl;

		rsd->atom(" OH ").formal_charge( -1 );
		rsd->atom(" O  ").formal_charge( -1 );

		rosetta_recharge_fullatom( *rsd );

		TS_ASSERT_DELTA( rsd->atom(" N  ").charge(), -0.470 - Y_naive_charge/21 + -2.0/21, 1e-4 ); //"Nbb "
		TS_ASSERT_DELTA( rsd->atom(" CA ").charge(),  0.070 - Y_naive_charge/21 + -2.0/21, 1e-4 ); //"CAbb"
		TS_ASSERT_DELTA( rsd->atom(" OH ").charge(), -0.660 - Y_naive_charge/21 + -2.0/21, 1e-4 ); //"OH  "

		net_charge = 0;
		for ( core::Size ii(1); ii <= rsd->natoms(); ++ii ) {
			net_charge += rsd->atom(ii).charge();
		}
		TS_ASSERT_DELTA( net_charge, -2.0, 1e-4 );

		TR << "Testing that virtual atoms aren't recharged and aren't counted during recharge" << std::endl;

		rsd = ResidueTypeOP( new ResidueType( rsd_types->name_map("PRO") ) );

		for ( core::Size ii(1); ii <= rsd->natoms(); ++ii ) {
			rsd->atom(ii).formal_charge( 0.0 );
			rsd->atom(ii).charge( -1.234 ); // Will be reset
		}

		rsd->atom(" CB ").formal_charge( -1 );
		rsd->atom(" CG ").formal_charge( -1 );

		rosetta_recharge_fullatom( *rsd );

		// 14 non-virtual in PRO:  Npro + CAbb + CObb + OCbb + CH2 * 3 + Hapo * 7
		core::Real P_naive_charge = -0.37 + 0.07 + 0.51 + -0.51 + -0.18 * 3 + 0.095 * 7;
		TS_ASSERT_DELTA( rsd->atom(" N  ").charge(), -0.370 - P_naive_charge/14 + -2.0/14, 1e-4 ); //"Nbb "
		TS_ASSERT_DELTA( rsd->atom(" CA ").charge(),  0.070 - P_naive_charge/14 + -2.0/14, 1e-4 ); //"CAbb"
		TS_ASSERT_EQUALS( rsd->atom(" NV ").charge(), 0 ); // VIRT
	}

	void test_make_centroid() {
		using namespace core::chemical;

		ChemicalManager * cm(ChemicalManager::get_instance());
		ResidueTypeSetCOP fa_rts = cm->residue_type_set( FA_STANDARD );
		ResidueTypeSetCOP cen_rts = cm->residue_type_set( CENTROID );

		utility::io::izstream paramslist("core/chemical/params/cen_types_list.txt");
		std::string cenfile, fafile;
		paramslist >> cenfile >> fafile;
		while ( paramslist ) {
			TR << "Comparing converted " << fafile << " with " << cenfile << std::endl;

			core::chemical::ResidueTypeOP fa_rsd = read_topology_file("core/chemical/params/"+fafile, fa_rts );
			core::chemical::ResidueTypeOP cen_rsd = read_topology_file("core/chemical/params/"+cenfile, cen_rts );

			core::chemical::ResidueTypeOP converted = make_centroid( *fa_rsd );

			TS_ASSERT_EQUALS( converted->atom_type_set().name(), cen_rsd->atom_type_set().name() );
			TS_ASSERT_EQUALS( converted->natoms(), cen_rsd->natoms() );

			for ( core::Size ii(1); ii <= cen_rsd->natoms(); ++ii ) {
				core::Size jj( converted->atom_index( cen_rsd->atom_name(ii) ) );
				TSM_ASSERT_EQUALS( cen_rsd->atom_name(ii), converted->atom_type( jj ).name(), cen_rsd->atom_type( ii ).name() );
			}
			paramslist >> cenfile >> fafile;
		}

	}

};
