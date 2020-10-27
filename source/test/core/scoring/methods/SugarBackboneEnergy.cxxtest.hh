// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    test/core/scoring/SugarBackboneEnergy.cxxtest.hh
/// @brief   Test suite for the carbohydrate SugarBackboneEnergy term
/// @author  Andrew Leaver-Fay <aleaverfay@gmail.com>

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Unit header
#include <core/energy_methods/SugarBackboneEnergy.hh>

// Package headers
#include <core/id/PartialAtomID.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/carbohydrates/util.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/carbohydrates/GlycanNode.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

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

// Numeric headers
#include <numeric/xyz.functions.hh>

// Basic headers
#include <basic/Tracer.hh>


static basic::Tracer TR( "core.pose.carbohydrates.util.cxxtest" );


class SugarBackboneEnergyTests : public CxxTest::TestSuite {
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

		//TR << *man9_op_ << std::endl;

	}

	// Destruction
	void tearDown()
	{}


public:  // Tests /////////////////////////////////////////////////////////////
	void test_atoms_w_dof_derivatives()
	{
		using namespace std;
		using namespace core::pose::carbohydrates;
		using namespace core::conformation::carbohydrates;
		using namespace core::scoring::methods::carbohydrates;
		using core::id::AtomID;
		using core::id::PartialAtomID;
		using core::scoring::EnergyMap;
		using core::Size;

		typedef utility::vector1< PartialAtomID > PartialAtomIDs;
		typedef utility::vector1< AtomID > AtomIDs;
		typedef utility::vector1< core::Vector > Vectors;

		TS_ASSERT( true );
		for ( core::Size ii = 1; ii <= bisected_man_.total_residue(); ++ii ) {
			std::cout << "Residue " << ii << " is " << bisected_man_.residue(ii).name() << std::endl;
		}
		SugarBackboneEnergy sbbe;
		EnergyMap emap;

		for ( Size ii = 1; ii <= bisected_man_.total_residue(); ++ii ) {
			sbbe.residue_energy(bisected_man_.residue(ii), bisected_man_, emap);
		}

		{
			PartialAtomIDs pids = sbbe.atoms_with_dof_derivatives(
				bisected_man_.residue(2),
				bisected_man_ );

			TS_ASSERT_EQUALS( pids.size(), 6 );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(8, 2), pids[1] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(1, 2), pids[2] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(1, 1, 0), pids[3] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(1, 1, 1), pids[4] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(1, 1, 2), pids[5] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(1, 1, 3), pids[6] );

			AtomIDs ids; ids.reserve( pids.size() );
			Vectors coords; coords.reserve( pids.size() );
			for ( core::Size ii = 1; ii <= pids.size(); ++ii ) {
				ids.push_back( bisected_man_.conformation().resolve_partial_atom_id( pids[ ii ] ));
				coords.push_back( bisected_man_.xyz( ids[ ii ] ));
			}

			utility::vector1< core::Real > dihedrals_actual(2);
			dihedrals_actual[1] = bisected_man_.phi(2);
			dihedrals_actual[2] = bisected_man_.psi(2);
			// residue 2's linkage is not exocyclic
			// dihedrals_actual[3] = bisected_man_.omega(2);

			for ( core::Size ii = 1; ii <= 2; ++ii ) {
				core::Real dihedral_from_pids = numeric::dihedral_degrees(
					coords[ii], coords[ii+1], coords[ii+2], coords[ii+3] );
				TS_ASSERT_DELTA( dihedrals_actual[ii], dihedral_from_pids, 1e-6 );
			}
		}

		{
			PartialAtomIDs pids = sbbe.atoms_with_dof_derivatives(
				bisected_man_.residue(3),
				bisected_man_ );
			// for ( auto const & id : ids ) {
			//  std::cout << "3 id " << id << std::endl;
			// }

			TS_ASSERT_EQUALS( pids.size(), 6 );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(8, 3), pids[1] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(1, 3), pids[2] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(2, 1, 0), pids[3] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(2, 1, 1), pids[4] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(2, 1, 2), pids[5] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(2, 1, 3), pids[6] );

			AtomIDs ids; ids.reserve( pids.size() );
			Vectors coords; coords.reserve( pids.size() );
			for ( core::Size ii = 1; ii <= pids.size(); ++ii ) {
				ids.push_back( bisected_man_.conformation().resolve_partial_atom_id( pids[ ii ] ));
				coords.push_back( bisected_man_.xyz( ids[ ii ] ));
			}

			utility::vector1< core::Real > dihedrals_actual(2);
			dihedrals_actual[1] = bisected_man_.phi(3);
			dihedrals_actual[2] = bisected_man_.psi(3);
			// residue 3's linkage is not exocyclic
			// dihedrals_actual[3] = bisected_man_.omega(3);

			for ( core::Size ii = 1; ii <= 2; ++ii ) {
				core::Real dihedral_from_pids = numeric::dihedral_degrees(
					coords[ii], coords[ii+1], coords[ii+2], coords[ii+3] );
				TS_ASSERT_DELTA( dihedrals_actual[ii], dihedral_from_pids, 1e-6 );
			}
		}

		{
			PartialAtomIDs pids = sbbe.atoms_with_dof_derivatives(
				bisected_man_.residue(4),
				bisected_man_ );
			// for ( auto const & id : ids ) {
			//  std::cout << "4 id " << id << std::endl;
			// }

			TS_ASSERT_EQUALS( pids.size(), 6 );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(8, 4), pids[1] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(1, 4), pids[2] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(3, 1, 0), pids[3] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(3, 1, 1), pids[4] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(3, 1, 2), pids[5] );
			TS_ASSERT_EQUALS( core::id::PartialAtomID(3, 1, 3), pids[6] );

			AtomIDs ids; ids.reserve( pids.size() );
			Vectors coords; coords.reserve( pids.size() );
			for ( core::Size ii = 1; ii <= pids.size(); ++ii ) {
				ids.push_back( bisected_man_.conformation().resolve_partial_atom_id( pids[ ii ] ));
				coords.push_back( bisected_man_.xyz( ids[ ii ] ));
			}

			utility::vector1< core::Real > dihedrals_actual(3);
			dihedrals_actual[1] = bisected_man_.phi(4);
			dihedrals_actual[2] = bisected_man_.psi(4);
			dihedrals_actual[3] = bisected_man_.omega(4);

			for ( core::Size ii = 1; ii <= 3; ++ii ) {
				core::Real dihedral_from_pids = numeric::dihedral_degrees(
					coords[ii], coords[ii+1], coords[ii+2], coords[ii+3] );
				TS_ASSERT_DELTA( dihedrals_actual[ii], dihedral_from_pids, 1e-6 );
			}
		}
	}


private:  // Private data /////////////////////////////////////////////////////
	core::pose::Pose Lex_;  // Lewisx: beta-D-Galp-(1->4)-[alpha-D-Fucp-(1->3)]-D-GlcpNAc
	core::pose::Pose isomaltose_;  // a (1alpha->6) disaccharide of D-glucose
	core::pose::Pose bisected_man_;  // a-D-Manp-(1->3)-[b-d-GlcpNAc-(1->4)]-[a-D-Manp-(1->6)]-b-D-Manp
	core::pose::Pose exo_test_; // alpha-L-Fucp-(1->6)-D-GlcpNAc-(1->4)-D-GlcpNAc
	core::pose::PoseOP man9_op_; // a-D-Manp-(1->2)-a-D-Manp-(1->2)-a-D-Manp-(1->3)-[a-D-Manp-(1->2)-a-D-Manp-(1->3)-[a-D-Manp-(1->2)-a-D-Manp-(1->6)]-a-D-Manp-(1->6)]-b-D-Manp-(1->4)-b-D-GlcpNAc-(1->4)-b-D-GlcpNAc

	core::pose::Pose two_glycans_; //Simple PDB with two man5s

};  // class CarbohydratePoseUtilityFunctionTests
