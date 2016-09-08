// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  test/protocols/membrane/benchmark/MakeCanonicalHelix.cxxtest.hh
/// @brief Unit Test for creating ideal helical peptides or custom helical peptides
/// @author Rebecca Faye Alford (rfalford12@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

// Project Headers
#include <protocols/membrane/benchmark/MakeCanonicalHelix.hh>
#include <protocols/moves/Mover.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/sequence/util.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static THREAD_LOCAL basic::Tracer TR("protocols.membrane.benchmark.MakeCanonicalHelixTest.cxxtest");

using namespace core;
using namespace core::pose;
using namespace core::chemical;

class MakeCanonicalHelixTest : public CxxTest::TestSuite {

public: // test functions

	void setUp() {

		// Initialize core & options system
		core_init();

		// Setup current FASTA sequence
		std::string sequence( "ALALALALALALALA" );

		TR <<  "Constructing pose from sequence: ALALALALALALALA"  << std::endl;
		// Grab the current residue typeset
		ResidueTypeSetCOP const & residue_set(
			ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
		);
		pose_ = PoseOP( new Pose() );
		make_pose_from_sequence( *pose_, sequence, *residue_set );

		// Make pose into an ideal helix
		TR <<  "Creating a helical peptide from the 'pose from sequence'"  << std::endl;
		using namespace protocols::membrane::benchmark;
		MakeCanonicalHelixOP helical_pept( new MakeCanonicalHelix() );
		helical_pept->apply( *pose_ );

	}

	void tearDown() {
		// Nothing to see here
	}

	void test_default_dihedral_angles() {

		TR <<  "Testing correct setup of ideal peptide from sequence with anticipated dihedral angles"  << std::endl;
		for ( core::Size i = 1; i <= pose_->size(); ++i ) {

			TS_ASSERT_EQUALS( pose_->phi( i ), -57.0 );
			TS_ASSERT_EQUALS( pose_->psi( i ), -47.0 );
			TS_ASSERT_EQUALS( pose_->omega( i ), 175.0 );
		}
	}

private:

	PoseOP pose_;

};

