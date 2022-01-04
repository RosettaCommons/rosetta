// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/CycpepRigidBodyPermutationMoverTests.cxxtest.hh
/// @brief  Unit tests for the CycpepRigidBodyPermutationMover.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>
#include <test/util/pose_funcs.hh>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/cyclic_peptide/CycpepRigidBodyPermutationMover.hh>

// Protocols Headers
#include <protocols/cyclic_peptide/DeclareBond.hh>

// Core Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/scoring/rms_util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID.hh>

// Utility Headers
#include <utility/vector1.hh>

// Basic Headers
#include <basic/Tracer.hh>

// STL Headers
#include <map>

static basic::Tracer TR("CycpepRigidBodyPermutationMoverTests");


class CycpepRigidBodyPermutationMoverTests : public CxxTest::TestSuite {
	//Define Variables

public:

	void setUp() {
		core_init();

	}

	void tearDown() {

	}

	/// @brief Compute RMSDs for every permutation to the original pose, and confirm that
	/// the applied permutation's RMSD is the lowest.
	bool min_perturbation_is(
		core::Size const permutation,
		bool const inverted,
		core::pose::Pose const & master_pose,
		core::pose::Pose const & pose
	) {
		bool best_inversion(false);
		core::Size best_permut(0);

		utility::vector1< utility::vector1< core::Size > > forward_permutations {
			{ 232, 233, 234, 235, 236, 237, 238, 239 },
			{ 233, 234, 235, 236, 237, 238, 239, 232 },
			{ 234, 235, 236, 237, 238, 239, 232, 233 },
			{ 235, 236, 237, 238, 239, 232, 233, 234 },
			{ 236, 237, 238, 239, 232, 233, 234, 235 },
			{ 237, 238, 239, 232, 233, 234, 235, 236 },
			{ 238, 239, 232, 233, 234, 235, 236, 237 },
			{ 239, 232, 233, 234, 235, 236, 237, 238 }
			};

		utility::vector1< utility::vector1< core::Size > > reverse_permutations {
			{ 232, 239, 238, 237, 236, 235, 234, 233 },
			{ 233, 232, 239, 238, 237, 236, 235, 234 },
			{ 234, 233, 232, 239, 238, 237, 236, 235 },
			{ 235, 234, 233, 232, 239, 238, 237, 236 },
			{ 236, 235, 234, 233, 232, 239, 238, 237 },
			{ 237, 236, 235, 234, 233, 232, 239, 238 },
			{ 238, 237, 236, 235, 234, 233, 232, 239 },
			{ 239, 238, 237, 236, 235, 234, 233, 232 },
			};

		core::Real lowest_rmsd(99999.9);
		utility::vector1< std::string > alignment_atoms{ "N", "CA", "C", "O", "CB" };

		for ( core::Size i(1); i<=8; ++i ) {
			std::map< core::id::AtomID, core::id::AtomID > atom_map;
			for ( core::Size j(1); j<=8; ++j ) {
				for ( std::string const & atname: alignment_atoms ) {
					atom_map[ core::id::AtomID( master_pose.residue_type( forward_permutations[1][j] ).atom_index(atname), forward_permutations[1][j] ) ] =
						core::id::AtomID( pose.residue_type( forward_permutations[i][j] ).atom_index( atname ), forward_permutations[i][j] );
				}
			}

			core::Real const rmsd( core::scoring::rms_at_corresponding_atoms_no_super( pose, master_pose, atom_map ) );
			if ( i == 0 || rmsd < lowest_rmsd ) {
				lowest_rmsd = rmsd;
				best_inversion = false;
				best_permut = i;
			}
			TR << "\tWithout inversion, offset " << i-1 << " yields RMSD " << rmsd << "." << std::endl;
		}

		for ( core::Size i(1); i<=8; ++i ) {
			std::map< core::id::AtomID, core::id::AtomID > atom_map;
			for ( core::Size j(1); j<=8; ++j ) {
				for ( std::string const & atname: alignment_atoms ) {
					atom_map[ core::id::AtomID( master_pose.residue_type( forward_permutations[1][j] ).atom_index(atname), forward_permutations[1][j] ) ] =
						core::id::AtomID( pose.residue_type( reverse_permutations[i][j] ).atom_index( atname ), reverse_permutations[i][j] );
				}
			}

			core::Real const rmsd( core::scoring::rms_at_corresponding_atoms_no_super( pose, master_pose, atom_map ) );
			if ( i == 0 || rmsd < lowest_rmsd ) {
				lowest_rmsd = rmsd;
				best_inversion = true;
				best_permut = i;
			}
			TR << "\tWith inversion, offset " << i-1 << " yields RMSD " << rmsd << "." << std::endl;
		}

		return best_permut == permutation + 1 && best_inversion == inverted;
	}

	/// @brief Perform every possible permutation of a cylic peptide.  For each one, confirm that
	/// the RMSD to the original pose with the appropriate permutation is less than the RMSD for every
	/// other permutation.
	void test_apply_setting() {
		using namespace protocols::cyclic_peptide;

		core::pose::PoseOP master_pose( core::import_pose::pose_from_file("protocols/cyclic_peptide/NDM1_bound_example.pdb") );
		DeclareBond decbond;
		decbond.set( 239, "C", 232, "N", false );
		decbond.apply(*master_pose);

		bool truefalse(false);
		for ( core::Size j(1); j<=2; ++j ) {
			for ( core::Size permut(0); permut <= 7; ++permut ) {
				TR << std::endl;
				TR << "Trying permutation=" << permut << " inversion=" << ( truefalse ? "TRUE" : "FALSE" ) << std::endl;

				core::select::residue_selector::ResidueIndexSelectorCOP selector(
					utility::pointer::make_shared< core::select::residue_selector::ResidueIndexSelector >("232-239")
				);

				CycpepRigidBodyPermutationMover my_mover;
				my_mover.set_mover_mode( "set_permutation" );
				my_mover.set_permutation_setting( permut );
				my_mover.set_inversion_setting( truefalse );
				my_mover.set_random_position_offset(0.0);
				my_mover.set_random_orientation_perturbation(0.0);
				my_mover.set_residue_selector(selector);

				core::pose::PoseOP pose( master_pose->clone() );
				my_mover.apply(*pose);
				TS_ASSERT( min_perturbation_is( permut, truefalse, *master_pose, *pose ) );
			}
			truefalse = true;
		}


	}


};
