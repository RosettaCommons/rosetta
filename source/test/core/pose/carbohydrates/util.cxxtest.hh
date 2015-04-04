// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file	 test/core/pose/carbohydrates/util.cxxtest.hh
/// @brief   Test suite for utility functions for carbohydrate-containing poses
/// @author  Labonte <JWLabonte@jhu.edu>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/pose/carbohydrates/util.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Project headers
#include <core/types.hh>
#include <core/id/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/vector1.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>


class CarbohydratePoseUtilityFunctionTests : public CxxTest::TestSuite {
public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		using namespace core::pose;
		using namespace core::import_pose;
		using namespace basic::options;

		core_init_with_additional_options( "-override_rsd_type_limit" );

		option[ OptionKeys::in::include_sugars ]( true );
		option[ OptionKeys::in::file::read_pdb_link_records ]( true );

		// Test branched oligosaccharide.
		pose_from_pdb( Lex_, "core/chemical/carbohydrates/Lex.pdb" );

		// Test oligosaccharide with exocyclic linkage.
		pose_from_pdb( isomaltose_, "core/chemical/carbohydrates/isomaltose.pdb" );
	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////
	void test_find_seqpos_of_saccharides_parent_residue()
	{
		using namespace core::pose::carbohydrates;

		TS_TRACE( "Testing find_seqpos_of_saccharides_parent_residue() function." );

		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_parent_residue( Lex_.residue( 1 ) ), 0 );  // 1st has no parent.
		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_parent_residue( Lex_.residue( 2 ) ), 1 );
		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_parent_residue( Lex_.residue( 3 ) ), 1 );  // 3rd is on a branch.

		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_parent_residue( isomaltose_.residue( 2 ) ), 1 );  // a ->6 linkage
	}

	void test_get_glycosidic_bond_residues()
	{
		using namespace std;
		using namespace core::conformation;
		using namespace core::pose::carbohydrates;

		TS_TRACE( "Testing get_glycosidic_bond_residues() function." );

		string first_residue;
		pair< ResidueCOP, ResidueCOP > residues;

		first_residue = Lex_.residue( 1 ).name();
		residues = get_glycosidic_bond_residues( Lex_, 1 );  // meaningless; should return residue 1 twice
		TS_ASSERT_EQUALS( residues.first, residues.second );

		residues = get_glycosidic_bond_residues( Lex_, 2 );  // a main-chain linkage
		TS_ASSERT_EQUALS( residues.second->name(), first_residue );

		residues = get_glycosidic_bond_residues( Lex_, 3 );  // a branch-point linkage
		TS_ASSERT_EQUALS( residues.second->name(), first_residue );


		first_residue = isomaltose_.residue( 1 ).name();
		residues = get_glycosidic_bond_residues( isomaltose_, 2 );  // a ->6 linkage
		TS_ASSERT_EQUALS( residues.second->name(), first_residue );
	}

	void test_get_reference_atoms()
	{
		using namespace core::id;
		using namespace core::chemical::carbohydrates;
		using namespace core::pose::carbohydrates;

		TS_TRACE( "Testing get_reference_atoms() function." );

		utility::vector1< AtomID > atoms;

		atoms = get_reference_atoms( phi_torsion, Lex_, 2 );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 1 ].rsd() ).atom_name( atoms[ 1 ].atomno() ), " VO5" );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 2 ].rsd() ).atom_name( atoms[ 2 ].atomno() ), " C1 " );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 3 ].rsd() ).atom_name( atoms[ 3 ].atomno() ), " O4 " );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 4 ].rsd() ).atom_name( atoms[ 4 ].atomno() ), " C4 " );

		atoms = get_reference_atoms( psi_torsion, Lex_, 2 );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 1 ].rsd() ).atom_name( atoms[ 1 ].atomno() ), " C1 " );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 2 ].rsd() ).atom_name( atoms[ 2 ].atomno() ), " O4 " );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 3 ].rsd() ).atom_name( atoms[ 3 ].atomno() ), " C4 " );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 4 ].rsd() ).atom_name( atoms[ 4 ].atomno() ), " C3 " );


		atoms = get_reference_atoms( phi_torsion, Lex_, 3 );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 1 ].rsd() ).atom_name( atoms[ 1 ].atomno() ), " VO5" );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 2 ].rsd() ).atom_name( atoms[ 2 ].atomno() ), " C1 " );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 3 ].rsd() ).atom_name( atoms[ 3 ].atomno() ), " O3 " );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 4 ].rsd() ).atom_name( atoms[ 4 ].atomno() ), " C3 " );

		atoms = get_reference_atoms( psi_torsion, Lex_, 3 );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 1 ].rsd() ).atom_name( atoms[ 1 ].atomno() ), " C1 " );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 2 ].rsd() ).atom_name( atoms[ 2 ].atomno() ), " O3 " );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 3 ].rsd() ).atom_name( atoms[ 3 ].atomno() ), " C3 " );
		TS_ASSERT_EQUALS( Lex_.residue( atoms[ 4 ].rsd() ).atom_name( atoms[ 4 ].atomno() ), " C2 " );


		atoms = get_reference_atoms( phi_torsion, isomaltose_, 2 );
		TS_ASSERT_EQUALS( isomaltose_.residue( atoms[ 1 ].rsd() ).atom_name( atoms[ 1 ].atomno() ), " VO5" );
		TS_ASSERT_EQUALS( isomaltose_.residue( atoms[ 2 ].rsd() ).atom_name( atoms[ 2 ].atomno() ), " C1 " );
		TS_ASSERT_EQUALS( isomaltose_.residue( atoms[ 3 ].rsd() ).atom_name( atoms[ 3 ].atomno() ), " O6 " );
		TS_ASSERT_EQUALS( isomaltose_.residue( atoms[ 4 ].rsd() ).atom_name( atoms[ 4 ].atomno() ), " C6 " );

		atoms = get_reference_atoms( psi_torsion, isomaltose_, 2 );
		TS_ASSERT_EQUALS( isomaltose_.residue( atoms[ 1 ].rsd() ).atom_name( atoms[ 1 ].atomno() ), " C1 " );
		TS_ASSERT_EQUALS( isomaltose_.residue( atoms[ 2 ].rsd() ).atom_name( atoms[ 2 ].atomno() ), " O6 " );
		TS_ASSERT_EQUALS( isomaltose_.residue( atoms[ 3 ].rsd() ).atom_name( atoms[ 3 ].atomno() ), " C6 " );
		TS_ASSERT_EQUALS( isomaltose_.residue( atoms[ 4 ].rsd() ).atom_name( atoms[ 4 ].atomno() ), " C5 " );

		atoms = get_reference_atoms( omega_torsion, isomaltose_, 2 );
		TS_ASSERT_EQUALS( isomaltose_.residue( atoms[ 1 ].rsd() ).atom_name( atoms[ 1 ].atomno() ), " O6 " );
		TS_ASSERT_EQUALS( isomaltose_.residue( atoms[ 2 ].rsd() ).atom_name( atoms[ 2 ].atomno() ), " C6 " );
		TS_ASSERT_EQUALS( isomaltose_.residue( atoms[ 3 ].rsd() ).atom_name( atoms[ 3 ].atomno() ), " C5 " );
		TS_ASSERT_EQUALS( isomaltose_.residue( atoms[ 4 ].rsd() ).atom_name( atoms[ 4 ].atomno() ), " C4 " );
	}

	void test_TorsionID_query_functions()
	{
		using namespace core::id;
		using namespace core::pose::carbohydrates;

		TS_TRACE( "Testing functions that query if the given TorsionID is of a glycosidic torsion." );

		// Are they phis?
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 1, BB, 3 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 1, BB, 4 ) ) );  // psi(2)
		TS_ASSERT( is_phi_torsion( Lex_, TorsionID( 1, BB, 5 ) ) );  // phi(2)

		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 1, CHI, 3 ) ) );  // psi(3)
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 1, CHI, 4 ) ) );  // psi(2)
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 1, CHI, 5 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 1, CHI, 6 ) ) );

		// Nus can never be phis!
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 1, NU, 1 ) ) );

		// Termini cannot have a phi(n+1)!
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 2, BB, 1 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 2, BB, 2 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 2, BB, 3 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 2, BB, 4 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 2, BB, 5 ) ) );

		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 2, CHI, 1 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 2, CHI, 3 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 2, CHI, 4 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 2, CHI, 5 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 2, CHI, 6 ) ) );

		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 3, BB, 1 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 3, BB, 2 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 3, BB, 3 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 3, BB, 4 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 3, BB, 5 ) ) );

		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 3, CHI, 1 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 3, CHI, 3 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 3, CHI, 4 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 3, CHI, 5 ) ) );
		TS_ASSERT( ! is_phi_torsion( Lex_, TorsionID( 3, CHI, 6 ) ) );

		TS_ASSERT( ! is_phi_torsion( isomaltose_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_phi_torsion( isomaltose_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_phi_torsion( isomaltose_, TorsionID( 1, BB, 3 ) ) );
		TS_ASSERT( ! is_phi_torsion( isomaltose_, TorsionID( 1, BB, 5 ) ) );  // omega(2)
		TS_ASSERT( ! is_phi_torsion( isomaltose_, TorsionID( 1, BB, 6 ) ) );  // psi(2)
		TS_ASSERT( is_phi_torsion( isomaltose_, TorsionID( 1, BB, 7 ) ) );  // phi(2)

		TS_ASSERT( ! is_phi_torsion( isomaltose_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( ! is_phi_torsion( isomaltose_, TorsionID( 1, CHI, 3 ) ) );
		TS_ASSERT( ! is_phi_torsion( isomaltose_, TorsionID( 1, CHI, 4 ) ) );
		TS_ASSERT( ! is_phi_torsion( isomaltose_, TorsionID( 1, CHI, 5 ) ) );  // omega(2)
		TS_ASSERT( ! is_phi_torsion( isomaltose_, TorsionID( 1, CHI, 6 ) ) );  // psi(2)


		// Are they psis?
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 1, BB, 3 ) ) );
		TS_ASSERT( is_psi_torsion( Lex_, TorsionID( 1, BB, 4 ) ) );  // psi(2)
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 1, BB, 5 ) ) );  // phi(2)

		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( is_psi_torsion( Lex_, TorsionID( 1, CHI, 3 ) ) );  // psi(3)
		//TS_ASSERT( is_psi_torsion( Lex_, TorsionID( 1, CHI, 4 ) ) );  // psi(2)
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 1, CHI, 5 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 1, CHI, 6 ) ) );

		// Nus can never be psis!
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 1, NU, 1 ) ) );

		// Termini cannot have a psi(n+1)!
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 2, BB, 1 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 2, BB, 2 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 2, BB, 3 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 2, BB, 4 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 2, BB, 5 ) ) );

		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 2, CHI, 1 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 2, CHI, 3 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 2, CHI, 4 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 2, CHI, 5 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 2, CHI, 6 ) ) );

		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 3, BB, 1 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 3, BB, 2 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 3, BB, 3 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 3, BB, 4 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 3, BB, 5 ) ) );

		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 3, CHI, 1 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 3, CHI, 3 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 3, CHI, 4 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 3, CHI, 5 ) ) );
		TS_ASSERT( ! is_psi_torsion( Lex_, TorsionID( 3, CHI, 6 ) ) );

		TS_ASSERT( ! is_psi_torsion( isomaltose_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_psi_torsion( isomaltose_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_psi_torsion( isomaltose_, TorsionID( 1, BB, 3 ) ) );
		TS_ASSERT( ! is_psi_torsion( isomaltose_, TorsionID( 1, BB, 5 ) ) );  // omega(2)
		TS_ASSERT( is_psi_torsion( isomaltose_, TorsionID( 1, BB, 6 ) ) );  // psi(2)
		TS_ASSERT( ! is_psi_torsion( isomaltose_, TorsionID( 1, BB, 7 ) ) );  // phi(2)

		TS_ASSERT( ! is_psi_torsion( isomaltose_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( ! is_psi_torsion( isomaltose_, TorsionID( 1, CHI, 3 ) ) );
		TS_ASSERT( ! is_psi_torsion( isomaltose_, TorsionID( 1, CHI, 4 ) ) );
		TS_ASSERT( ! is_psi_torsion( isomaltose_, TorsionID( 1, CHI, 5 ) ) );  // omega(2)
		//TS_ASSERT( is_psi_torsion( isomaltose_, TorsionID( 1, CHI, 6 ) ) );  // psi(2)


		// Are they omegas?
		TS_ASSERT( ! is_omega_torsion( Lex_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_omega_torsion( Lex_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_omega_torsion( Lex_, TorsionID( 1, BB, 3 ) ) );
		TS_ASSERT( ! is_omega_torsion( Lex_, TorsionID( 1, BB, 4 ) ) );  // psi(2)
		TS_ASSERT( ! is_omega_torsion( Lex_, TorsionID( 1, BB, 5 ) ) );  // phi(2)

		TS_ASSERT( ! is_omega_torsion( Lex_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( ! is_omega_torsion( Lex_, TorsionID( 1, CHI, 3 ) ) );  // psi(3)
		TS_ASSERT( ! is_omega_torsion( Lex_, TorsionID( 1, CHI, 4 ) ) );  // psi(2)
		TS_ASSERT( ! is_omega_torsion( Lex_, TorsionID( 1, CHI, 5 ) ) );
		TS_ASSERT( ! is_omega_torsion( Lex_, TorsionID( 1, CHI, 6 ) ) );

		// Nus can never be omegas!
		TS_ASSERT( ! is_omega_torsion( Lex_, TorsionID( 1, NU, 1 ) ) );

		// Termini cannot have an omega(n+1)!
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 2, BB, 1 ) ) );
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 2, BB, 2 ) ) );
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 2, BB, 3 ) ) );
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 2, BB, 6 ) ) );
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 2, BB, 7 ) ) );

		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 2, CHI, 1 ) ) );
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 2, CHI, 3 ) ) );
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 2, CHI, 4 ) ) );
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 2, CHI, 5 ) ) );
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 2, CHI, 6 ) ) );

		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 1, BB, 3 ) ) );
		TS_ASSERT( is_omega_torsion( isomaltose_, TorsionID( 1, BB, 5 ) ) );  // omega(2)
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 1, BB, 6 ) ) );  // psi(2)
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 1, BB, 7 ) ) );  // phi(2)

		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 1, CHI, 3 ) ) );
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 1, CHI, 4 ) ) );
		//TS_ASSERT( is_omega_torsion( isomaltose_, TorsionID( 1, CHI, 5 ) ) );  // omega(2)
		TS_ASSERT( ! is_omega_torsion( isomaltose_, TorsionID( 1, CHI, 6 ) ) );  // psi(2)
	}


private:  // Private datum ////////////////////////////////////////////////////
	core::pose::Pose Lex_;  // Lewisx: beta-D-Galp-(1->4)-[alpha-D-Fucp-(1->3)]-D-GlcpNAc
	core::pose::Pose isomaltose_;  // a (1alpha->6) disaccharide of D-glucose

};  // class CarbohydratePoseUtilityFunctionTests
