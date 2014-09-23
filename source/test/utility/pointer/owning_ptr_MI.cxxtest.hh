// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/access_ptr.cxxtest.hh
/// @brief  Unit tests for utility::pointer::owning_ptr with multiple inheritance
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Headers
//#include <iostream> // Put this first to allow debug print-outs in project headers
#include <set>

// Project headers
#include <cxxtest/TestSuite.h>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCountMI.hh>
#include <utility/pointer/ReferenceCountMI_.hh>


// --- set up some helper classes/functions for these tests
// --- (these go in a special namespace to be used from this file only!)
namespace owning_ptr_MI_test {

class Atom :
	virtual public utility::pointer::ReferenceCountMI
{
public:

	Atom()
	{
		++nAtoms;
	}

	virtual
	~Atom()
	{
		--nAtoms;
	}

	static int nAtoms; // For tracking destructions
};
int Atom::nAtoms = 0;


class Atom_ :
	public Atom,
	public utility::pointer::ReferenceCountMI_
{
public:

	Atom_( double const weight, double const charge ) :
		weight_( weight ),
		charge_( charge )
	{}

	virtual
	~Atom_()
	{}

	double weight_;
	double charge_;
};

} // namespace owning_ptr_MI_test


// The tests live in the testutility::pointer namespace, but
// grant them access to our local namespace for classes/functions.

using namespace owning_ptr_MI_test;


class OwningPtrMITests : public CxxTest::TestSuite {

	public:

	/// @brief Explicit owning pointer
	void test_shared_ptr_MI_explicit() {

		typedef  utility::pointer::shared_ptr< Atom >  AtomP;
		Atom * Cp = new Atom_( 12.0, 0.0 ); // Raw heap pointer

		TS_ASSERT( Atom::nAtoms == 1 );
		TS_ASSERT( Cp->ref_count() == 0 );
		{
			AtomP C( Cp ); // First owning pointer has sole ownership
			TS_ASSERT( Cp->ref_count() == 1 );
			{
				AtomP C2( Cp ); // Second owning pointer shares ownership
				TS_ASSERT( Cp->ref_count() == 2 );
				{
					AtomP C3( C2 ); // Third owning pointer created from 2nd and shares ownership
					TS_ASSERT( Cp->ref_count() == 3 );
				}
			}
			TS_ASSERT( Cp->ref_count() == 1 );
		}
		// Atom should have been destructed
		TS_ASSERT( Atom::nAtoms == 0 );
	}

};


