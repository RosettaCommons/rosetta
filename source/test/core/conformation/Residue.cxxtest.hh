// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/Residue.cxxtest.hh
/// @brief  test suite for core::conformation::Residue
/// @author Christopher Miles (cmiles@uw.edu)

// Test Headers
#include <cxxtest/TestSuite.h>

// Unit Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/PseudoBond.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>

// Project headers
#include <test/core/init_util.hh>
#include <test/util/pose_funcs.hh>

// Utility Headers
#include <utility/vector1.hh>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


using core::pose::Pose;
using core::conformation::PseudoBond;
using core::conformation::PseudoBondCollection;
using core::conformation::PseudoBondCollectionOP;
using core::conformation::PseudoBondCollectionCOP;
using core::conformation::Residue;
using core::conformation::ResidueOP;

class ResidueTest : public CxxTest::TestSuite {
public:
	void setUp() {
		core_init();
	}

	void test_isDNA() {
		Pose pose;
		core::import_pose::pose_from_file(pose, "core/conformation/4gatA.pdb", core::import_pose::PDB_file);


		unsigned dna_start = 68;
		unsigned dna_end = 93;

		for ( unsigned i = 1; i < dna_start; ++i ) {
			TS_ASSERT(!pose.residue(i).is_DNA());
		}

		for ( unsigned i = dna_start; i <= dna_end; ++i ) {
			TS_ASSERT(pose.residue(i).is_DNA());
		}
	}

	void test_isLigand() {
		Pose pose;
		core::import_pose::pose_from_file(pose, "core/conformation/4gatA.pdb", core::import_pose::PDB_file);
		TS_ASSERT(pose.residue(67).is_ligand());
	}

	void test_residue_serialization() {
		TS_ASSERT( true );
#ifdef SERIALIZATION
		Pose trpcage = create_trpcage_ideal_pose();
		ResidueOP trp6 = trpcage.residue( 6 ).clone();

		// Let's pretend that residue 7 is only a single atom large so we
		// can add some pseudobonds between residue 6 and residue 8, and
		// then test whether pseudobonds are properly serialized.
		PseudoBond pb68;
		pb68.lr( 6 ); pb68.ur( 8 );
		pb68.lr_conn_id( 2 ); // connecting at C on residue 6
		pb68.ur_conn_id( 1 ); // connecting at N on residue 8
		PseudoBondCollectionOP pbc( new PseudoBondCollection );
		pbc->push_back( pb68 );
		trp6->set_pseudobonds_to_residue( 8, pbc );

		TS_ASSERT( trp6->aa() == core::chemical::aa_trp ); // this should be trp in the input pose; not a test, so much as an assertion.

		// Now serialize the coordinates
		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arch( oss );
			arch( trp6 );
		}

		ResidueOP trp6_copy;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arch( iss );
			arch( trp6_copy );
		}

		// trp6_copy should be exactly the same as trp6
		// and it should point to the same (global) instance of
		// the fullatom tryptophan residue type
		TS_ASSERT( & trp6->type() == & trp6_copy->type() );

		// The code below walks through the data members of class Residue in order

		// Coordinates and dihedrals should match exactly
		for ( core::Size ii = 1; ii <= trp6->natoms(); ++ii ) {
			for ( core::Size jj = 1; jj <= 3; ++jj ) {
				TS_ASSERT_EQUALS( trp6->xyz( ii )( jj ), trp6_copy->xyz( ii )( jj ) );
			}
			TS_ASSERT_EQUALS( trp6->atom( ii ).type(),    trp6_copy->atom( ii ).type() );
			TS_ASSERT_EQUALS( trp6->atom( ii ).mm_type(), trp6_copy->atom( ii ).mm_type() );
		}

		// Sadly, skipping the orbitals since I'm unsure how to activate them (I'm sure it's easy).
		// ADD CODE FOR ORBITALS

		TS_ASSERT_EQUALS( trp6->seqpos(), trp6_copy->seqpos() );
		TS_ASSERT_EQUALS( trp6->chain(),  trp6_copy->chain() );

		TS_ASSERT_EQUALS( trp6->chi().size(), trp6_copy->chi().size() );
		for ( core::Size ii = 1; ii <= trp6->chi().size(); ++ii ) {
			TS_ASSERT_EQUALS( trp6->chi()[ ii ], trp6_copy->chi()[ ii ] );
		}

		// now, this comparison isn't going to be all that useful because
		// trp doesn't have any nu dihedrals, but let's include the code
		// anyways so it could be compared in the future with some other
		// nu-containing residue
		TS_ASSERT_EQUALS( trp6->nus().size(), trp6_copy->nus().size() );
		for ( core::Size ii = 1; ii <= trp6->nus().size(); ++ii ) {
			TS_ASSERT_EQUALS( trp6->nus()[ ii ], trp6_copy->nus()[ ii ] );
		}

		TS_ASSERT_EQUALS( trp6->mainchain_torsions().size(), trp6_copy->mainchain_torsions().size() );
		for ( core::Size ii = 1; ii <= trp6->mainchain_torsions().size(); ++ii ) {
			TS_ASSERT_EQUALS( trp6->mainchain_torsions()[ ii ], trp6_copy->mainchain_torsions()[ ii ] );
		}

		for ( core::Size ii = 1; ii <= 3; ++ii ) {
			TS_ASSERT_EQUALS( trp6->actcoord()( ii ), trp6_copy->actcoord()( ii ) );
		}

		// ADD CODE HERE TO COMPARE DATA CACHE!

		// huh... can't directly access the nonstandard_polymer_ data member.

		// Make sure that the inter-residue connection information has been perfectly transfered
		// Connect map
		for ( core::Size ii = 1; ii <= trp6->n_residue_connections(); ++ii ) {
			TS_ASSERT_EQUALS( trp6->actual_residue_connection( ii ).resid(), trp6_copy->actual_residue_connection( ii ).resid() );
			TS_ASSERT_EQUALS( trp6->actual_residue_connection( ii ).connid(), trp6_copy->actual_residue_connection( ii ).connid() );
		}

		TS_ASSERT_EQUALS( trp6->connections_to_residue( 5 ), trp6_copy->connections_to_residue( 5 ) );
		TS_ASSERT_EQUALS( trp6->connections_to_residue( 7 ), trp6_copy->connections_to_residue( 7 ) );
		for ( core::Size ii = 1; ii <= trpcage.total_residue(); ++ii ) {
			TS_ASSERT_EQUALS( trp6_copy->is_bonded( ii ), ii == 5 || ii == 7 ); // the only inter-residue bonds are to 5 and 7
			TS_ASSERT_EQUALS( trp6_copy->is_pseudo_bonded( ii ), ii == 8 ); // the only pseudobond is to residue 8
		}

		PseudoBondCollectionCOP trp6_pbs      = trp6->get_pseudobonds_to_residue( 8 );
		PseudoBondCollectionCOP trp6_copy_pbs = trp6_copy->get_pseudobonds_to_residue( 8 );

		TS_ASSERT_EQUALS( trp6_pbs->size(), trp6_copy_pbs->size() );
		for ( PseudoBondCollection::PBIter trp6_pb = trp6_pbs->iter_begin(),
				trp6_copy_pb = trp6_copy_pbs->iter_begin();
				trp6_pb != trp6_pbs->iter_end(); ++trp6_pb, ++trp6_copy_pb ) {
			TS_ASSERT_EQUALS( trp6_pb->lr(), trp6_copy_pb->lr() );
			TS_ASSERT_EQUALS( trp6_pb->ur(), trp6_copy_pb->ur() );
			TS_ASSERT_EQUALS( trp6_pb->lr_conn_id(), trp6_copy_pb->lr_conn_id() );
			TS_ASSERT_EQUALS( trp6_pb->ur_conn_id(), trp6_copy_pb->ur_conn_id() );
		}

#endif // SERIALIZATION
	}


};


