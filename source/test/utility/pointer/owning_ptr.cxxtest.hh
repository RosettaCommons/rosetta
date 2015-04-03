// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/owning_ptr.cxxtest.hh
/// @brief  Unit tests for utility::pointer::owning_ptr
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Headers
#include <iostream> // Put this first to allow debug print-outs in project headers
#include <set>

// Project headers
#include <cxxtest/TestSuite.h>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>


// --- set up some helper classes/functions for these tests
// --- (these go in a special namespace to be used from this file only!)

namespace owning_ptr_test {

class Atom_ :
	public utility::pointer::ReferenceCount
{
public:

	Atom_( double const weight, double const charge ) :
		weight_( weight ),
		charge_( charge )
	{
		++nAtoms;
	}

	virtual
	~Atom_()
	{
		--nAtoms;
	}

	double weight_;
	double charge_;
	static int nAtoms; // For tracking destructions
}; // class Atom_

int Atom_::nAtoms = 0;

typedef  utility::pointer::shared_ptr< Atom_ >  AtomP_;
typedef  std::set< AtomP_ >  Atoms_;


class Molecule_ :
	public utility::pointer::ReferenceCount
{
public:

//	typedef  utility::pointer::owning_ptr< Atom_ >  AtomP_;
//	typedef  std::set< AtomP_ >  Atoms_;

	inline
	Molecule_()
	{}

	virtual
	~Molecule_()
	{}

	inline
	void
	add( Atom_ * atom_p )
	{
		atoms_.insert( atom_p );
	}

	inline
	void
	add( AtomP_ const & atom_p )
	{
		atoms_.insert( atom_p );
	}

	inline
	std::set< AtomP_ >::size_type
	n_atom() const
	{
		return atoms_.size();
	}

	Atoms_ atoms_;

}; // class Molecule_


} // namespace owning_ptr_test


// The tests live in the utility::pointer namespace, but
// grant them access to our local namespace for classes/functions.
using namespace owning_ptr_test;


class OwningPtrTests : public CxxTest::TestSuite {

	public:

	/// @brief Explicit owning pointer
	void test_shared_ptr_explicit() {
		Atom_ * Cp = new Atom_( 12.0, 0.0 ); // Raw heap pointer
		TS_ASSERT_EQUALS( Atom_::nAtoms, 1 );
		TS_ASSERT( Cp->ref_count() == 0 );
		{
			AtomP_ C( Cp ); // First owning pointer has sole ownership
			TS_ASSERT( Cp->ref_count() == 1 );
			{
				AtomP_ C2( Cp ); // Second owning pointer shares ownership
				TS_ASSERT( Cp->ref_count() == 2 );
			}
			TS_ASSERT( Cp->ref_count() == 1 );
		}
		// Atom should have been destructed
		TS_ASSERT( Atom_::nAtoms == 0 );
	}

	/// @brief Containers
	void test_shared_ptr_containers() {
		Atom_ * Cp = new Atom_( 13.0, 0.0 ); // Raw heap pointer
		Atom_ * Hp = new Atom_( 1.0, 0.0 ); // Raw heap pointer
		TS_ASSERT( Atom_::nAtoms == 2 );
		TS_ASSERT( Cp->ref_count() == 0 );
		TS_ASSERT( Hp->ref_count() == 0 );
		{
			Molecule_ m1; m1.add( Cp ); m1.add( Hp );
			TS_ASSERT( m1.n_atom() == 2 );
			TS_ASSERT( Cp->ref_count() == 1 );
			TS_ASSERT( Hp->ref_count() == 1 );
			{
				Molecule_ m2; m2.add( Cp ); m2.add( Hp );
				TS_ASSERT( m2.n_atom() == 2 );
				TS_ASSERT( Cp->ref_count() == 2 );
				TS_ASSERT( Hp->ref_count() == 2 );
				m2.add( Cp ); // Repeat: Won't be added again
				TS_ASSERT( m2.n_atom() == 2 );
				TS_ASSERT( Cp->ref_count() == 2 );
			}
			TS_ASSERT( Cp->ref_count() == 1 );
			TS_ASSERT( Hp->ref_count() == 1 );
		}
		// Atoms should have been destructed
		TS_ASSERT( Atom_::nAtoms == 0 );
	}

};


