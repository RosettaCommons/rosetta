// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constaints/AmbiguousNMRDistanceConstraint.cxxtest.hh
/// @brief  test suite for AmbiguousNMRDistanceConstraint
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>

#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh> // DO NOT AUTO-REMOVE

//protocol headers
#include <protocols/cyclic_peptide/PeptideStubMover.hh>

//core headers
#include <core/pose/Pose.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>

// basic headers

// unit headers

#ifdef SERIALIZATION
#include <core/id/AtomID.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif


class AmbiguousNMRDistanceConstraintTests : public CxxTest::TestSuite
{
private:
    std::string const test_label_ = "PEPTIDE_STUB_MOVER_TEST_LABEL";
public:

    void setUp(){
    }

    void tearDown(){

	}

    void test_ala() {
        typedef utility::vector1< core::id::AtomID > Atoms;
        core::pose::Pose pose( (core::pose::Pose()) );

        protocols::cyclic_peptide::PeptideStubMover stubmover;
        stubmover.set_reset_mode(false);
        stubmover.add_residue( "Append", "ALA", 0, false, "", 0, pose.total_residue(), nullptr, "" );
        stubmover.set_residue_label( test_label_ );
        stubmover.apply(pose);

        //utility::vector1< core::id::AtomID > atoms_old;
        //utility::vector1< core::id::AtomID > atoms_new;
        Atoms atoms_old, atoms_new;

        core::Size res (1);
        core::chemical::AA const aa( pose.residue_type( res ).aa() );

        core::scoring::constraints::parse_NMR_name("HA", res, aa, atoms_old);
        core::scoring::constraints::parse_NMR_name("HA", res, aa, atoms_new);

        TS_ASSERT( atoms_new == atoms_old );
    }

	void test_serialize_AmbiguousNMRDistanceConstraint() {
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		using namespace core::scoring::constraints;
		using namespace core::scoring::func;
		using namespace core::id;

		FuncOP some_func( new HarmonicFunc( 1, 2 ));
		AmbiguousNMRDistanceConstraint::Atoms atoms1, atoms2;
		atoms1.push_back( AtomID( 1, 2 )); atoms1.push_back( AtomID( 1, 10 ) );
		atoms2.push_back( AtomID( 2, 3 )); atoms1.push_back( AtomID( 4, 16 ) );

		ConstraintOP instance( new AmbiguousNMRDistanceConstraint( atoms1, atoms2, some_func ) );

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( instance );
		}

		ConstraintOP instance2;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}

		TS_ASSERT( utility::pointer::dynamic_pointer_cast< AmbiguousNMRDistanceConstraint > ( instance2 ) );
		TS_ASSERT( *instance == *instance2 );
#endif // SERIALIZATION
	}

};

