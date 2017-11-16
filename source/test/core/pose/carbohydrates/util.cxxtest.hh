// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/pose/carbohydrates/util.cxxtest.hh
/// @brief   Test suite for utility functions for carbohydrate-containing poses
/// @author  Labonte <JWLabonte@jhu.edu>
/// @author  Jared Adolf-Bryfogle <jadolfbr@gmail.com>


// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/pose/carbohydrates/util.hh>
#include <core/conformation/carbohydrates/util.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>

// Package headers
#include <core/conformation/Residue.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>

// Project headers
#include <core/types.hh>
#include <core/id/types.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// Basic headers
#include <basic/Tracer.hh>


static basic::Tracer TR( "core.pose.carbohydrates.util.cxxtest" );


class CarbohydratePoseUtilityFunctionTests : public CxxTest::TestSuite {
public:  // Standard methods //////////////////////////////////////////////////
	// Initialization
	void setUp()
	{
		using namespace core::pose;
		using namespace core::import_pose;


		core_init_with_additional_options( "-include_sugars" );

		// Test branched oligosaccharide.
		pose_from_file( Lex_, "core/chemical/carbohydrates/Lex.pdb" , PDB_file);

		// Test oligosaccharide with exocyclic linkage.
		pose_from_file( isomaltose_, "core/chemical/carbohydrates/isomaltose.pdb", PDB_file );

		// Test oligosaccharide with multiple branches off a single residue.
		make_pose_from_saccharide_sequence( bisected_man_,
			"a-D-Manp-(1->3)-[a-D-Manp-(1->6)]-[b-d-GlcpNAc-(1->4)]-b-D-Manp" );

		// Test exocyclic carbon in linkage.
		pose_from_file( exo_test_,
			"core/chemical/carbohydrates/alpha-L-Fucp-_1-6_-D-GlcpNAc-_1-4_-D-GlcpNAc.pdb", PDB_file);

		std::string const man9_s( "a-D-Manp-(1->2)-a-D-Manp-(1->2)-a-D-Manp-(1->3)-[a-D-Manp-(1->2)-a-D-Manp-(1->3)-"
			"[a-D-Manp-(1->2)-a-D-Manp-(1->6)]-a-D-Manp-(1->6)]-b-D-Manp-(1->4)-b-D-GlcpNAc-(1->4)-b-D-GlcpNAc" );
		man9_op_ = pose_from_saccharide_sequence( man9_s, "fa_standard", true, false ); //No need to idealize.

		TR << *man9_op_ << std::endl;
	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////
	void test_find_seqpos_of_saccharides_relatives()
	{
		using namespace std;
		using namespace core::pose::carbohydrates;
		using namespace core::conformation::carbohydrates;

		TR << "Testing find_seqpos_of_saccharides_parent_residue() function." << endl;

		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_parent_residue( Lex_.residue( 1 ) ), 0 );  // 1st has no parent.
		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_parent_residue( Lex_.residue( 2 ) ), 1 );
		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_parent_residue( Lex_.residue( 3 ) ), 1 );  // 3rd is on a branch.

		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_parent_residue( isomaltose_.residue( 2 ) ), 1 );  // a ->6 linkage


		TR << "Testing find_seqpos_of_saccharides_child_residue_at() function." << endl;

		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_child_residue_at( bisected_man_.residue( 1 ), 1 ), 0 );
		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_child_residue_at( bisected_man_.residue( 1 ), 2 ), 0 );
		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_child_residue_at( bisected_man_.residue( 1 ), 3 ), 2 );  // main
		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_child_residue_at( bisected_man_.residue( 1 ), 4 ), 3 );
		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_child_residue_at( bisected_man_.residue( 1 ), 5 ), 4 );  // as O6
		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_child_residue_at( bisected_man_.residue( 1 ), 6 ), 4 );
		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_child_residue_at( bisected_man_.residue( 1 ), 7 ), 0 );  // OoB

		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_child_residue_at( bisected_man_.residue( 2 ), 4 ), 0 );  // term.
		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_child_residue_at( bisected_man_.residue( 3 ), 4 ), 0 );  // term.
		TS_ASSERT_EQUALS( find_seqpos_of_saccharides_child_residue_at( bisected_man_.residue( 4 ), 4 ), 0 );  // term.
	}

	void test_get_glycosidic_bond_residues()
	{
		using namespace std;
		using namespace core::conformation;
		using namespace core::pose::carbohydrates;
		using namespace core::conformation::carbohydrates;

		TR <<  "Testing get_glycosidic_bond_residues() function."  << std::endl;

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

	void test_get_linkage_position_of_saccharide_residue()
	{
		using namespace core::pose::carbohydrates;
		using namespace core::conformation::carbohydrates;

		TR <<  "Testing get_linkage_position_of_saccharide_residue() function."  << std::endl;
		TS_ASSERT_EQUALS( Lex_.glycan_tree_set()->get_linkage_position(1), 0 );
		TS_ASSERT_EQUALS( Lex_.glycan_tree_set()->get_linkage_position(2), 4 );
		TS_ASSERT_EQUALS( Lex_.glycan_tree_set()->get_linkage_position(3), 3 );
	}

	void test_get_reference_atoms()
	{
		using namespace core::id;
		using namespace core::chemical::carbohydrates;
		using namespace core::pose::carbohydrates;
		using namespace core::conformation::carbohydrates;

		TR <<  "Testing get_reference_atoms() function."  << std::endl;

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
		using namespace std;
		using namespace core::id;
		using namespace core::pose::carbohydrates;
		using namespace core::conformation::carbohydrates;

		TR << "Testing functions that query if the given TorsionID is of a glycosidic torsion." << endl;

		TR << " Testing recognition of phi..." << endl;

		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 1, BB, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 1, BB, 4 ) ) );  // psi(2)
		TS_ASSERT( is_glycosidic_phi_torsion( Lex_, TorsionID( 1, BB, 5 ) ) );  // phi(2)

		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 1, BB, 6 ) ) );  // out of bounds

		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 1, CHI, 3 ) ) );  // psi(3)
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 1, CHI, 4 ) ) );  // psi(2)
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 1, CHI, 5 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 1, CHI, 6 ) ) );

		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 1, CHI, 7 ) ) );  // out of bounds

		// Nus can never be phis!
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 1, NU, 1 ) ) );

		TS_ASSERT( is_glycosidic_phi_torsion( Lex_, TorsionID( 1, BRANCH, 1 ) ) );  // phi(3)
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 1, BRANCH, 2 ) ) );  // out of bounds

		// Termini cannot have a phi(n+1)!
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 2, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 2, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 2, BB, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 2, BB, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 2, BB, 5 ) ) );

		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 2, CHI, 1 ) ) );  // virtual-atom-moving torsion
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 2, CHI, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 2, CHI, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 2, CHI, 5 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 2, CHI, 6 ) ) );

		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 3, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 3, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 3, BB, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 3, BB, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 3, BB, 5 ) ) );

		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 3, CHI, 1 ) ) );  // virtual-atom-moving torsion
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 3, CHI, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 3, CHI, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 3, CHI, 5 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( Lex_, TorsionID( 3, CHI, 6 ) ) );

		// Test a system with an exocyclic linkage.
		TS_ASSERT( ! is_glycosidic_phi_torsion( isomaltose_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( isomaltose_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( isomaltose_, TorsionID( 1, BB, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( isomaltose_, TorsionID( 1, BB, 5 ) ) );  // omega(2)
		TS_ASSERT( ! is_glycosidic_phi_torsion( isomaltose_, TorsionID( 1, BB, 6 ) ) );  // psi(2)
		TS_ASSERT( is_glycosidic_phi_torsion( isomaltose_, TorsionID( 1, BB, 7 ) ) );  // phi(2)

		TS_ASSERT( ! is_glycosidic_phi_torsion( isomaltose_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( isomaltose_, TorsionID( 1, CHI, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( isomaltose_, TorsionID( 1, CHI, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( isomaltose_, TorsionID( 1, CHI, 5 ) ) );  // omega(2)
		TS_ASSERT( ! is_glycosidic_phi_torsion( isomaltose_, TorsionID( 1, CHI, 6 ) ) );  // psi(2)

		// Test a system with multiple branches.
		TS_ASSERT( ! is_glycosidic_phi_torsion( bisected_man_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( bisected_man_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( bisected_man_, TorsionID( 1, BB, 3 ) ) );  // psi(2)
		TS_ASSERT( is_glycosidic_phi_torsion( bisected_man_, TorsionID( 1, BB, 4 ) ) );  // phi(2)

		TS_ASSERT( ! is_glycosidic_phi_torsion( bisected_man_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_phi_torsion( bisected_man_, TorsionID( 1, CHI, 3 ) ) );  // virtual-atom-moving
		TS_ASSERT( ! is_glycosidic_phi_torsion( bisected_man_, TorsionID( 1, CHI, 4 ) ) );  // psi(3)
		TS_ASSERT( ! is_glycosidic_phi_torsion( bisected_man_, TorsionID( 1, CHI, 5 ) ) );  // omega(4)
		TS_ASSERT( ! is_glycosidic_phi_torsion( bisected_man_, TorsionID( 1, CHI, 6 ) ) );  // psi(4)

		TS_ASSERT( is_glycosidic_phi_torsion( bisected_man_, TorsionID( 1, BRANCH, 1 ) ) );  // phi(3)
		TS_ASSERT( is_glycosidic_phi_torsion( bisected_man_, TorsionID( 1, BRANCH, 2 ) ) );  // phi(4)


		TR << " Testing recognition of psi..." << endl;

		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 1, BB, 3 ) ) );
		TS_ASSERT( is_glycosidic_psi_torsion( Lex_, TorsionID( 1, BB, 4 ) ) );  // psi(2)
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 1, BB, 5 ) ) );  // phi(2)

		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 1, BB, 6 ) ) );  // out of bounds

		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( is_glycosidic_psi_torsion( Lex_, TorsionID( 1, CHI, 3 ) ) );  // psi(3)
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 1, CHI, 4 ) ) );  // virtual-atom-moving torsion
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 1, CHI, 5 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 1, CHI, 6 ) ) );

		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 1, CHI, 7 ) ) );  // out of bounds

		// Nus can never be psis!
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 1, NU, 1 ) ) );

		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 1, BRANCH, 1 ) ) );  // phi(3)
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 1, BRANCH, 2 ) ) );  // out of bounds

		// Termini cannot have a psi(n+1)!
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 2, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 2, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 2, BB, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 2, BB, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 2, BB, 5 ) ) );

		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 2, CHI, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 2, CHI, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 2, CHI, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 2, CHI, 5 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 2, CHI, 6 ) ) );

		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 3, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 3, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 3, BB, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 3, BB, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 3, BB, 5 ) ) );

		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 3, CHI, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 3, CHI, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 3, CHI, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 3, CHI, 5 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( Lex_, TorsionID( 3, CHI, 6 ) ) );

		// Test a system with an exocyclic linkage.
		TS_ASSERT( ! is_glycosidic_psi_torsion( isomaltose_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( isomaltose_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( isomaltose_, TorsionID( 1, BB, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( isomaltose_, TorsionID( 1, BB, 5 ) ) );  // omega(2)
		TS_ASSERT( is_glycosidic_psi_torsion( isomaltose_, TorsionID( 1, BB, 6 ) ) );  // psi(2)
		TS_ASSERT( ! is_glycosidic_psi_torsion( isomaltose_, TorsionID( 1, BB, 7 ) ) );  // phi(2)

		TS_ASSERT( ! is_glycosidic_psi_torsion( isomaltose_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( isomaltose_, TorsionID( 1, CHI, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( isomaltose_, TorsionID( 1, CHI, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( isomaltose_, TorsionID( 1, CHI, 5 ) ) );  // omega(2)
		TS_ASSERT( ! is_glycosidic_psi_torsion( isomaltose_, TorsionID( 1, CHI, 6 ) ) );  // virtual-atom-moving torsion

		// Test a system with multiple branches.
		TS_ASSERT( ! is_glycosidic_psi_torsion( bisected_man_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( bisected_man_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( is_glycosidic_psi_torsion( bisected_man_, TorsionID( 1, BB, 3 ) ) );  // psi(2)
		TS_ASSERT( ! is_glycosidic_psi_torsion( bisected_man_, TorsionID( 1, BB, 4 ) ) );  // phi(2)

		TS_ASSERT( ! is_glycosidic_psi_torsion( bisected_man_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_psi_torsion( bisected_man_, TorsionID( 1, CHI, 3 ) ) );  // virtual-atom-moving
		TS_ASSERT( is_glycosidic_psi_torsion( bisected_man_, TorsionID( 1, CHI, 4 ) ) );  // psi(3)
		TS_ASSERT( ! is_glycosidic_psi_torsion( bisected_man_, TorsionID( 1, CHI, 5 ) ) );  // omega(4)
		TS_ASSERT( is_glycosidic_psi_torsion( bisected_man_, TorsionID( 1, CHI, 6 ) ) );  // psi(4)

		TS_ASSERT( ! is_glycosidic_psi_torsion( bisected_man_, TorsionID( 1, BRANCH, 1 ) ) );  // phi(3)
		TS_ASSERT( ! is_glycosidic_psi_torsion( bisected_man_, TorsionID( 1, BRANCH, 2 ) ) );  // phi(4)


		TR << " Testing recognition of omega1..." << endl;

		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, BB, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, BB, 4 ) ) );  // psi(2)
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, BB, 5 ) ) );  // phi(2)

		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, BB, 6 ) ) );  // out of bounds

		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, CHI, 3 ) ) );  // psi(3)
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, CHI, 4 ) ) );  // psi(2)
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, CHI, 5 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, CHI, 6 ) ) );

		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, CHI, 7 ) ) );  // out of bounds

		// Nus can never be omegas!
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, NU, 1 ) ) );

		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, BRANCH, 1 ) ) );  // phi(3)
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 1, BRANCH, 2 ) ) );  // out of bounds

		// Termini cannot have an omega(n+1)!
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 2, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 2, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 2, BB, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 2, BB, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 2, BB, 5 ) ) );

		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 2, CHI, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 2, CHI, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 2, CHI, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 2, CHI, 5 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 2, CHI, 6 ) ) );

		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 3, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 3, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 3, BB, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 3, BB, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 3, BB, 5 ) ) );

		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 3, CHI, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 3, CHI, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 3, CHI, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 3, CHI, 5 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( Lex_, TorsionID( 3, CHI, 6 ) ) );

		// Test a system with an exocyclic linkage.
		TS_ASSERT( ! is_glycosidic_omega_torsion( isomaltose_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( isomaltose_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( isomaltose_, TorsionID( 1, BB, 3 ) ) );
		TS_ASSERT( is_glycosidic_omega_torsion( isomaltose_, TorsionID( 1, BB, 5 ) ) );  // omega(2)
		TS_ASSERT( ! is_glycosidic_omega_torsion( isomaltose_, TorsionID( 1, BB, 6 ) ) );  // psi(2)
		TS_ASSERT( ! is_glycosidic_omega_torsion( isomaltose_, TorsionID( 1, BB, 7 ) ) );  // phi(2)

		TS_ASSERT( ! is_glycosidic_omega_torsion( isomaltose_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( isomaltose_, TorsionID( 1, CHI, 3 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( isomaltose_, TorsionID( 1, CHI, 4 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( isomaltose_, TorsionID( 1, CHI, 5 ) ) );  // omega(2), already covered
		TS_ASSERT( ! is_glycosidic_omega_torsion( isomaltose_, TorsionID( 1, CHI, 6 ) ) );  // virtual-atom-moving

		// Test a system with multiple branches.
		TS_ASSERT( ! is_glycosidic_omega_torsion( bisected_man_, TorsionID( 1, BB, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( bisected_man_, TorsionID( 1, BB, 2 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( bisected_man_, TorsionID( 1, BB, 3 ) ) );  // psi(2)
		TS_ASSERT( ! is_glycosidic_omega_torsion( bisected_man_, TorsionID( 1, BB, 4 ) ) );  // phi(2)

		TS_ASSERT( ! is_glycosidic_omega_torsion( bisected_man_, TorsionID( 1, CHI, 1 ) ) );
		TS_ASSERT( ! is_glycosidic_omega_torsion( bisected_man_, TorsionID( 1, CHI, 3 ) ) );  // virtual-atom-moving
		TS_ASSERT( ! is_glycosidic_omega_torsion( bisected_man_, TorsionID( 1, CHI, 4 ) ) );  // psi(3)
		TS_ASSERT( is_glycosidic_omega_torsion( bisected_man_, TorsionID( 1, CHI, 5 ) ) );  // omega(4)
		TS_ASSERT( ! is_glycosidic_omega_torsion( bisected_man_, TorsionID( 1, CHI, 6 ) ) );  // psi(4)

		TS_ASSERT( ! is_glycosidic_omega_torsion( bisected_man_, TorsionID( 1, BRANCH, 1 ) ) );  // phi(3)
		TS_ASSERT( ! is_glycosidic_omega_torsion( bisected_man_, TorsionID( 1, BRANCH, 2 ) ) );  // phi(4)
	}

	void test_get_downstream_residue_that_this_torsion_moves()
	{
		using namespace core::id;
		using namespace core::pose::carbohydrates;
		using namespace core::conformation::carbohydrates;


		TR << "Testing get_downstream_residue_that_this_torsion_moves() function." << std::endl;

		// All main-chain torsions move the next residue in the main chain.
		TS_ASSERT_EQUALS( get_downstream_residue_that_this_torsion_moves( bisected_man_, TorsionID( 1, BB, 1 ) ), 2 );
		TS_ASSERT_EQUALS( get_downstream_residue_that_this_torsion_moves( bisected_man_, TorsionID( 1, BB, 2 ) ), 2 );
		TS_ASSERT_EQUALS( get_downstream_residue_that_this_torsion_moves( bisected_man_, TorsionID( 1, BB, 3 ) ), 2 );
		TS_ASSERT_EQUALS( get_downstream_residue_that_this_torsion_moves( bisected_man_, TorsionID( 1, BB, 4 ) ), 2 );

		TS_ASSERT_EQUALS( get_downstream_residue_that_this_torsion_moves( bisected_man_, TorsionID( 1, CHI, 1 ) ), 0 );
		TS_ASSERT_EQUALS( get_downstream_residue_that_this_torsion_moves( bisected_man_, TorsionID( 1, CHI, 2 ) ), 0 );
		TS_ASSERT_EQUALS( get_downstream_residue_that_this_torsion_moves( bisected_man_, TorsionID( 1, CHI, 3 ) ), 2 );
		TS_ASSERT_EQUALS( get_downstream_residue_that_this_torsion_moves( bisected_man_, TorsionID( 1, CHI, 4 ) ), 3 );
		TS_ASSERT_EQUALS( get_downstream_residue_that_this_torsion_moves( bisected_man_, TorsionID( 1, CHI, 5 ) ), 4 );
		TS_ASSERT_EQUALS( get_downstream_residue_that_this_torsion_moves( bisected_man_, TorsionID( 1, CHI, 6 ) ), 4 );

		TS_ASSERT_EQUALS(
			get_downstream_residue_that_this_torsion_moves( bisected_man_, TorsionID( 1, BRANCH, 1 ) ), 3 );
		TS_ASSERT_EQUALS(
			get_downstream_residue_that_this_torsion_moves( bisected_man_, TorsionID( 1, BRANCH, 2 ) ), 4 );

		TS_ASSERT_EQUALS( get_downstream_residue_that_this_torsion_moves( bisected_man_, TorsionID( 1, NU, 1 ) ), 0 );
	}

	void test_glycan_movemap_setup()
	{
		using namespace core::select::residue_selector;
		using namespace core::pose::carbohydrates;
		using namespace core::kinematics;
		using namespace core::id;

		//Main Chain 1->2
		utility::vector1< bool > subset( bisected_man_.size(), false);
		subset[2] = true; //Single Glycan residue.  Make sure all are enabled.

		ReturnResidueSubsetSelectorOP return_subset = ReturnResidueSubsetSelectorOP( new ReturnResidueSubsetSelector( subset));

		bool move_ring = true;
		bool move_bb = true;
		bool move_chi = true;

		MoveMapOP mm = create_glycan_movemap_from_residue_selector( bisected_man_, return_subset, move_chi, move_ring, move_bb );

		//mm->show(TR);

		TS_ASSERT( mm->get( TorsionID( 1, BB, 4) ));
		TS_ASSERT( mm->get( TorsionID( 1, BB, 3) ));

		//TS_ASSERT( mm->get( TorsionID( 1, CHI, 1) ) );
		//TS_ASSERT( mm->get( TorsionID( 1, CHI, 2) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, CHI, 3) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, CHI, 4) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, CHI, 5) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, CHI, 6) ) );

		TS_ASSERT( ! mm->get( TorsionID( 1, BRANCH, 1) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, BRANCH, 2) ) );


		//Branch Connection 1->3

		subset[2] = false;
		subset[3] = true;

		return_subset->set_residue_subset(subset);
		mm = create_glycan_movemap_from_residue_selector( bisected_man_, return_subset, move_chi, move_ring, move_bb );
		//mm->show(TR);

		TS_ASSERT( ! mm->get( TorsionID( 1, BB, 4) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, BB, 3) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, BB, 2) ) );
		TS_ASSERT( ! mm->get(   TorsionID( 1, CHI, 1) ) );
		TS_ASSERT( ! mm->get(   TorsionID( 1, CHI, 2) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, CHI, 3) ) );
		TS_ASSERT( mm->get(   TorsionID( 1, CHI, 4) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, CHI, 5) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, CHI, 6) ) );

		TS_ASSERT( mm->get(   TorsionID( 1, BRANCH, 1) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, BRANCH, 2) ) );
		TS_ASSERT( mm->get_nu( 3 ) );


		// Testing CHI
		subset[1] = true;
		subset[2] = false;
		subset[3] = false;

		move_ring = false;
		move_bb = false;
		move_chi = true;

		return_subset->set_residue_subset(subset);
		mm = create_glycan_movemap_from_residue_selector( bisected_man_, return_subset, move_chi, move_ring, move_bb );

		TS_ASSERT( ! mm->get( TorsionID( 1, BB, 4) ));
		TS_ASSERT( ! mm->get( TorsionID( 1, BB, 3) ));

		TS_ASSERT( mm->get( TorsionID( 1, CHI, 1) ) );
		TS_ASSERT( mm->get( TorsionID( 1, CHI, 2) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, CHI, 3) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, CHI, 4) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, CHI, 5) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, CHI, 6) ) );

		TS_ASSERT( ! mm->get( TorsionID( 1, BRANCH, 1) ) );
		TS_ASSERT( ! mm->get( TorsionID( 1, BRANCH, 2) ) );

	}
	void test_exocyclic_detection()
	{
		using namespace core::id;
		using namespace core::pose::carbohydrates;
		using namespace core::conformation::carbohydrates;

		TR << "Testing exocyclic linkage detection " << std::endl;
		//Can be improved to test branches as well.
		TS_ASSERT_THROWS_NOTHING( has_exocyclic_glycosidic_linkage( exo_test_.conformation(), 1) ); //Make sure we don't crash
		TS_ASSERT( ! has_exocyclic_glycosidic_linkage(exo_test_.conformation(), 1) ); //Make sure we get false.

		TS_ASSERT( exo_test_.glycan_tree_set()->has_exocyclic_glycosidic_linkage( 3 ) );
		TS_ASSERT( !  exo_test_.glycan_tree_set()->has_exocyclic_glycosidic_linkage(  2 ) );

	}

	void test_glycan_leafs()
	{
		using core::Size;
		using namespace utility;
		using namespace core::pose::carbohydrates;
		using namespace core::conformation::carbohydrates;

		TR << "MAN9 res: " << man9_op_->size() << std::endl;
		TS_ASSERT(man9_op_->size() == 11);

		Size tip_9[] = {9, 8};
		Size tip_11[] = {11, 10};
		Size tip_6[] = {6, 5, 4};

		//Test Leaf from Tip:
		utility::vector1< Size > branch_tip_9 (tip_9, tip_9 + sizeof(tip_9) / sizeof(Size) );
		utility::vector1< Size > branch_tip_11 (tip_11, tip_11 + sizeof(tip_11) / sizeof(Size) );
		utility::vector1< Size > branch_tip_6 (tip_6, tip_6 + sizeof(tip_6) / sizeof(Size) );

		utility::vector1< Size > branch_tip_9_test = get_resnums_in_leaf(*man9_op_, 9 /* tip */, 3 /* stop at */);
		utility::vector1< Size > branch_tip_11_test = get_resnums_in_leaf(*man9_op_, 11 /* tip */, 3 /* stop at */);
		utility::vector1< Size > branch_tip_6_test = get_resnums_in_leaf(*man9_op_, 6 /* tip */, 3 /* stop at */);

		TR << "Tip 9: " << utility::to_string (branch_tip_9_test) << std::endl;
		TR << "Tip 11: " << utility::to_string (branch_tip_11_test) << std::endl;
		TR << "Tip 6: " << utility::to_string (branch_tip_6_test) << std::endl;

		TS_ASSERT_EQUALS( get_resnums_in_leaf(*man9_op_, 9 /* tip */, 3 /* stop at */), branch_tip_9);
		TS_ASSERT_EQUALS( get_resnums_in_leaf(*man9_op_, 11 /* tip */, 3 /* stop at */), branch_tip_11);
		TS_ASSERT_EQUALS( get_resnums_in_leaf(*man9_op_, 6 /* tip */, 3 /* stop at */), branch_tip_6);

		TR.flush();
	}

	void test_glycan_branch()
	{
		using core::Size;
		using namespace utility;
		using namespace core::pose::carbohydrates;
		using namespace core::conformation::carbohydrates;

		TR << "Testing glycan branch" << std::endl;
		TS_ASSERT(man9_op_->size() == 11);

		//Test getting all residues and tips from specific positions.
		std::pair< vector1< Size >, vector1< Size > > res_and_tips;

		res_and_tips = get_carbohydrate_residues_and_tips_of_branch(man9_op_->conformation(), 3);
		TR << "Tips up of 3: " << res_and_tips.second << std::endl;
		TR << "Resn up of 3: " << res_and_tips.first << std::endl;
		Size tips3[] = {6, 9, 11};
		Size res3[] = {4, 7, 5, 6, 8, 10, 9, 11};
		utility::vector1< Size > tips3_v (tips3, tips3 + sizeof(tips3) / sizeof(Size) );
		utility::vector1< Size > res3_v (res3, res3 + sizeof(res3) / sizeof(Size) );
		TS_ASSERT_EQUALS( tips3_v, res_and_tips.second);
		TS_ASSERT_EQUALS( res3_v, res_and_tips.first);

		res_and_tips = get_carbohydrate_residues_and_tips_of_branch(man9_op_->conformation(), 7);
		TR << "Tips up of 7: " << res_and_tips.second << std::endl;
		TR << "Resn up of 7: " << res_and_tips.first << std::endl;
		Size tips7[] = {9, 11};
		Size res7[] = {8, 10, 9, 11};
		utility::vector1< Size > tips7_v (tips7, tips7 + sizeof(tips7) / sizeof(Size) );
		utility::vector1< Size > res7_v (res7, res7 + sizeof(res7) / sizeof(Size) );
		TS_ASSERT_EQUALS( tips7_v, res_and_tips.second);
		TS_ASSERT_EQUALS( res7_v, res_and_tips.first);

		res_and_tips = get_carbohydrate_residues_and_tips_of_branch(man9_op_->conformation(), 8);
		TR << "Tips up of 8: " << res_and_tips.second << std::endl;
		TR << "Resn up of 8: " << res_and_tips.first << std::endl;
		Size tips8[] = {9};
		Size res8[] = {9};
		utility::vector1< Size > tips8_v (tips8, tips8 + sizeof(tips8) / sizeof(Size) );
		utility::vector1< Size > res8_v (res8, res8 + sizeof(res8) / sizeof(Size) );
		TS_ASSERT_EQUALS( tips8_v, res_and_tips.second);
		TS_ASSERT_EQUALS( res8_v, res_and_tips.first);

		res_and_tips = get_carbohydrate_residues_and_tips_of_branch(man9_op_->conformation(), 4);
		TR << "Tips up of 4: " << res_and_tips.second << std::endl;
		TR << "Resn up of 4: " << res_and_tips.first << std::endl;
		Size tips4[] = {6};
		Size res4[] = {5, 6};
		utility::vector1< Size > tips4_v (tips4, tips4 + sizeof(tips4) / sizeof(Size) );
		utility::vector1< Size > res4_v (res4, res4 + sizeof(res4) / sizeof(Size) );
		TS_ASSERT_EQUALS( tips4_v, res_and_tips.second);
		TS_ASSERT_EQUALS( res4_v, res_and_tips.first);

		res_and_tips = get_carbohydrate_residues_and_tips_of_branch(man9_op_->conformation(), 2);
		TR << "Tips up of 2: " << res_and_tips.second << std::endl;
		TR << "Resn up of 2: " << res_and_tips.first << std::endl;
		Size tips2[] = {6, 9, 11};
		Size res2[] = {3, 4, 7, 5, 6, 8, 10, 9, 11};
		utility::vector1< Size > tips2_v (tips2, tips2 + sizeof(tips2) / sizeof(Size) );
		utility::vector1< Size > res2_v (res2, res2 + sizeof(res2) / sizeof(Size) );
		TS_ASSERT_EQUALS( tips2_v, res_and_tips.second);
		TS_ASSERT_EQUALS( res2_v, res_and_tips.first);


		/* Output:
		core.pose.carbohydrates.util: Children: [4, 7]
		Tips up of 3: [6, 9, 11]
		Resn up of 3: [4, 7, 5, 6, 8, 10, 9, 11]
		core.pose.carbohydrates.util: Children: [8, 10]
		Tips up of 7: [9, 11]
		Resn up of 7: [8, 10, 9, 11]
		core.pose.carbohydrates.util: Children: []
		Tips up of 8: [9]
		Resn up of 8: [9]
		core.pose.carbohydrates.util: Children: [5]
		Tips up of 4: [6]
		Resn up of 4: [5, 6]
		core.pose.carbohydrates.util: Children: [3]
		Tips up of 2: [6, 9, 11]
		Resn up of 2: [3, 4, 7, 5, 6, 8, 10, 9, 11]
		*/

		TR.flush();
	}

	void test_glycan_trimming()
	{
		using core::Size;
		using namespace utility;
		using namespace core::pose::carbohydrates;
		using namespace core::conformation::carbohydrates;

		man9_op_->dump_pdb("man9_pose.pdb");
		TR << "testing glycan trimming!" << std::endl;
		// Test actual trimming of glycan.

		TR << std::endl << "Deleting from 8" << std::endl;
		core::pose::PoseOP man9_copy = man9_op_->clone();
		TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 8); //Deletes 9
		TS_ASSERT_EQUALS(man9_copy->size(), 10);
		man9_copy->dump_pdb("man9_trim_at_8.pdb");

		TR << std::endl << "Deleting from 10" << std::endl;
		man9_copy = man9_op_->clone();
		TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 10); //Deletes 11
		TS_ASSERT_EQUALS(man9_copy->size(), 10);
		man9_copy->dump_pdb("man9_trim_at_10.pdb");

		TR << std::endl << "Deleting from 4" << std::endl;
		man9_copy = man9_op_->clone();
		TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 4); //Deletes 6,5
		TS_ASSERT_EQUALS(man9_copy->size(), 9);
		man9_copy->dump_pdb("man9_trim_at_4.pdb");

		TR << std::endl << "Deleting from 7" << std::endl;
		man9_copy = man9_op_->clone();
		TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 7); //Deletes 9,8,11,10
		TS_ASSERT_EQUALS(man9_copy->size(), 7);
		man9_copy->dump_pdb("man9_trim_at_7.pdb");

		TR << std::endl << "Deleting from 3" << std::endl;
		man9_copy = man9_op_->clone();
		TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 3); //Deletes 9,8,11,10,6,5,4
		TS_ASSERT_EQUALS(man9_copy->size(), 3);
		man9_copy->dump_pdb("man9_trim_at_3.pdb");

		TR << std::endl << "Deleting from 2" << std::endl;
		man9_copy = man9_op_->clone();
		TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 2); ///Deletes 9,8,11,10,6,5,4,3
		TS_ASSERT_EQUALS(man9_copy->size(), 2);
		man9_copy->dump_pdb("man9_trim_at_2.pdb");

		TR << std::endl << "Deleting from 1" << std::endl;
		man9_copy = man9_op_->clone();
		TS_ASSERT(man9_copy->size() ==11);
		delete_carbohydrate_branch( *man9_copy, 1); ///Deletes 9,8,11,10,6,5,4,3, 2
		TS_ASSERT_EQUALS(man9_copy->size(), 1);
		man9_copy->dump_pdb("man9_trim_at_1.pdb");

		TR.flush();
	}

	void test_glycosylate_pose() {
		core::pose::Pose N_linked_14_mer;
		core::pose::make_pose_from_sequence( N_linked_14_mer, "ANASA", "fa_standard" );
		core::pose::carbohydrates::glycosylate_pose( N_linked_14_mer, 2,
			"a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Glcp-(1->3)-a-D-Manp-(1->2)-a-D-Manp-(1->2)-a-D-Manp-(1->3)-"
			"[a-D-Manp-(1->2)-a-D-Manp-(1->3)-[a-D-Manp-(1->2)-a-D-Manp-(1->6)]-a-D-Manp-(1->6)]-b-D-Manp-(1->4)-"
			"b-D-GlcpNAc-(1->4)-b-D-GlcpNAc-" );

		// Saccharide gets reversed and reordered, such that it's
		//  14-13-12-11-10-9-
		//  [17-16-[19-18]-15]-8-
		// 7-6-
		TS_ASSERT_EQUALS( N_linked_14_mer.size(), 19 );
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(2).name(), "ASN:N-glycosylated");
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(6).name(), "->4)-beta-D-Glcp:2-AcNH");
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(8).name(), "->3)-beta-D-Manp:->6)-branch");
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(9).name(), "->2)-alpha-D-Manp");
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(14).name(), "->4)-alpha-D-Glcp:non-reducing_end");
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(15).name(), "->3)-alpha-D-Manp:->6)-branch");
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(17).name(), "->4)-alpha-D-Manp:non-reducing_end");
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(19).name(), "->4)-alpha-D-Manp:non-reducing_end");

		TS_ASSERT_EQUALS( N_linked_14_mer.residue(6).connected_residue_at_lower(), 2 );
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(9).connected_residue_at_lower(), 8 );
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(15).connected_residue_at_lower(), 8 );
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(16).connected_residue_at_lower(), 15 );
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(18).connected_residue_at_lower(), 15 );
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(8).connected_residue_at_upper(), 9 );
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(15).connected_residue_at_upper(), 16 );
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(2).connected_residue_at_resconn(3), 6 );
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(8).connected_residue_at_resconn(3), 15 );
		TS_ASSERT_EQUALS( N_linked_14_mer.residue(15).connected_residue_at_resconn(3), 18 );

		TS_ASSERT_EQUALS( N_linked_14_mer.fold_tree().to_string(),
			"FOLD_TREE  EDGE 1 5 -1  EDGE 2 6 -2  ND2  C1   EDGE 6 14 -1  EDGE 8 15 -2  O6   C1   EDGE 15 17 -1  EDGE 15 18 -2  O6   C1   EDGE 18 19 -1 ");

	}


private:  // Private data /////////////////////////////////////////////////////
	core::pose::Pose Lex_;  // Lewisx: beta-D-Galp-(1->4)-[alpha-D-Fucp-(1->3)]-D-GlcpNAc
	core::pose::Pose isomaltose_;  // a (1alpha->6) disaccharide of D-glucose
	core::pose::Pose bisected_man_;  // a-D-Manp-(1->3)-[b-d-GlcpNAc-(1->4)]-[a-D-Manp-(1->6)]-b-D-Manp
	core::pose::Pose exo_test_; // alpha-L-Fucp-(1->6)-D-GlcpNAc-(1->4)-D-GlcpNAc
	core::pose::PoseOP man9_op_; // a-D-Manp-(1->2)-a-D-Manp-(1->2)-a-D-Manp-(1->3)-[a-D-Manp-(1->2)-a-D-Manp-(1->3)-[a-D-Manp-(1->2)-a-D-Manp-(1->6)]-a-D-Manp-(1->6)]-b-D-Manp-(1->4)-b-D-GlcpNAc-(1->4)-b-D-GlcpNAc
};  // class CarbohydratePoseUtilityFunctionTests
