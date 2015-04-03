// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/pointer/access_ptr.cxxtest.hh
/// @brief  Unit tests for utility::pointer::access_ptr
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Headers
#include <iostream>      // Put this first to allow debug print-outs in project headers
#include <set>
#include <cxxtest/TestSuite.h>

// Project headers
#include <utility/pointer/access_ptr.hh>


// --- set up some helper classes/functions for these tests
// --- (these go in a special namespace to be used from this file only!)
namespace access_ptr_test {

// Forward
class Bond_;

class Atom_
{

public:

	typedef  utility::pointer::weak_ptr< Bond_ >  BondP_;
	typedef  utility::pointer::weak_ptr< Bond_ const >  ConstBondP_;
	typedef  std::set< BondP_ >  Bonds_;

	Atom_( double const weight, double const charge ) :
		weight_( weight ),
		charge_( charge )
	{}

	friend
	inline
	bool
	operator ==( Atom_ const & a, Atom_ const & b )
	{
		return ( ( a.weight_ == b.weight_ ) && ( a.charge_ == b.charge_ ) );
	}

	void
	add( Bond_ & bond )
	{
//		bonds_.insert( BondP_( bond ) );
		bonds_.insert( &bond );
	}

	void
	add( BondP_ const & bp )
	{
		bonds_.insert( bp );
	}

	void
	remove( Bond_ & bond )
	{
//		bonds_.erase( BondP_( bond ) );
		bonds_.erase( &bond );
	}

	void
	remove( BondP_ const & bp )
	{
		bonds_.erase( bp );
	}

//	bool
//	has( Bond_ const & bond ) const
//	{
//		return ( bonds_.find( ConstBondP_( &bond ) ) != bonds_.end() );
//	}

	bool
	has( Bond_ & bond ) const
	{
		return ( bonds_.find( BondP_( &bond ) ) != bonds_.end() );
	}

	double weight_;
	double charge_;
	Bonds_ bonds_;
};


class Bond_
{
public:

	typedef  utility::pointer::weak_ptr< Atom_ >  AtomP_;

	Bond_( Atom_ & a, Atom_ & b ) :
		a_( a ),
		b_( b )
	{
		a.add( *this );
		b.add( *this );
	}

	friend
	inline
	bool
	operator <( Bond_ const & a, Bond_ const & b )
	{
		return ( &a < &b );
	}

	friend
	inline
	bool
	operator ==( Bond_ const & a, Bond_ const & b )
	{
		return ( ( a.a_ == b.a_ ) && ( a.b_ == b.b_ ) );
	}

	virtual
	~Bond_()
	{
		a_->remove( *this );
		b_->remove( *this );
	}

	AtomP_ a_;
	AtomP_ b_;
};


class SpecialBond_ :
	public Bond_
{
public:

	SpecialBond_( Atom_ & a, Atom_ & b ) :
		Bond_( a, b )
	{}

};


void
f_utility::pointer::weak_ptr< Bond_ > const & bp )
{
	assert( bp );
}


void
g_utility::pointer::weak_ptr< Bond_ > bp )
{
	assert( bp );
}


void
h_utility::pointer::weak_ptr< Bond_ > & bp )
{
	assert( bp );
}

} // namespace access_ptr_test


// The tests live in the utility::pointer namespace, but
// grant them access to our local namespace for classes/functions.

using namespace access_ptr_test;


class AccessPtrTests : public CxxTest::TestSuite {

	public:

	/// @brief Size test
	void test_weak_ptr_size() {
		TS_ASSERT( sizeof( utility::pointer::access_ptr< Bond_ > ) == sizeof( Bond_ * ) );
	}

	/// @brief Atom/Bond basics
	void test_weak_ptr_basics() {
		Atom_ C( 12.0, 0.0 ), H( 1.0, 0.0 );
		Bond_ bond( C, H );
		Bond_ bond2( C, H );
		TS_ASSERT( C.bonds_.size() == 2 );
		TS_ASSERT( H.bonds_.size() == 2 );
		TS_ASSERT( C.has( bond ) );
		TS_ASSERT( C.has( bond2 ) );
		TS_ASSERT( H.has( bond ) );
		TS_ASSERT( H.has( bond2 ) );
		TS_ASSERT( bond.a_ );
		TS_ASSERT( bond.b_ );
		TS_ASSERT( bond2.a_ );
		TS_ASSERT( bond2.b_ );
		TS_ASSERT( *(bond.a_) == C );
		TS_ASSERT( *(bond.b_) == H );
		TS_ASSERT( bond.a_->weight_ == C.weight_ );
		TS_ASSERT( bond.b_->weight_ == H.weight_ );
	}

	/// @brief Constructor test
	void test_weak_ptr_constructor() {
		Atom_ C( 12.0, 0.0 ), H( 1.0, 0.0 );
		Bond_ bond( C, H );
		utility::pointer::weak_ptr< Bond_ > CH( bond );
		TS_ASSERT( C.bonds_.find( CH ) != C.bonds_.end() );
		TS_ASSERT( *C.bonds_.find( CH ) == CH );
		TS_ASSERT( C.bonds_.size() == 1 );
	}

	/// @brief Derived object
	void test_weak_ptr_derived() {
		Atom_ C( 12.0, 0.0 ), H( 1.0, 0.0 );
		SpecialBond_ bond( C, H );
		utility::pointer::weak_ptr< Bond_ > CH( bond ); // Can hold pointer to derived class object
		TS_ASSERT( C.bonds_.find( CH ) != C.bonds_.end() );
		TS_ASSERT( *C.bonds_.find( CH ) == CH );
		TS_ASSERT( C.bonds_.size() == 1 );
	}

	/// @brief Const object test
	void test_weak_ptr_const_object() {
		Atom_ C( 12.0, 0.0 ), H( 1.0, 0.0 );
		Bond_ bond( C, H );
		utility::pointer::weak_ptr< Bond_ const > CH( bond );
		TS_ASSERT( CH->a_->weight_ == C.weight_ );
		//CH->a_.reset(); // Won't compile because CH holds pointer to const
	}

	/// @brief Derived object conversion
	void test_weak_ptr_derived_conversion() {
		Atom_ C( 12.0, 0.0 ), H( 1.0, 0.0 );
		Bond_ bond( C, H );
		utility::pointer::weak_ptr< Bond_ > CH( bond );
		f_( CH );
		SpecialBond_ special_bond( C, H );
		utility::pointer::weak_ptr< SpecialBond_ > sCH( special_bond );
		f_( sCH ); // OK: Creates temporary bound to const reference
		g_( sCH ); // OK: Creates temporary bound to passed value
		// h_( sCH ); // Won't compile since can't create non-const temporary
	}

	/// @brief Derived object construction and assignment
	void test_weak_ptr_derived_construction_assignment() {
		Atom_ C( 12.0, 0.0 ), H( 1.0, 0.0 );
		Bond_ bond( C, H );
		SpecialBond_ special_bond( C, H );
		utility::pointer::weak_ptr< SpecialBond_ > sCH( special_bond );
		utility::pointer::weak_ptr< Bond_ > CH( sCH );
		TS_ASSERT( CH == sCH );
		CH = sCH;
		TS_ASSERT( CH == sCH );
	}

	/// @brief Const object construction and assignment
	void test_weak_ptr_constant_construction_assignment() {
		Atom_ H( 1.0, 0.0 ), C( 12.0, 0.0 );
		utility::pointer::weak_ptr< Atom_ > Cp( C );
		utility::pointer::weak_ptr< Atom_ const > Ccp( Cp );
		TS_ASSERT( Ccp->weight_ == C.weight_ );
		// access_ptr< Atom_ > badCp( Ccp ); // Won't compile since pointer to const Atom_ isn't assignable to pointer to non-const Atom_
		utility::pointer::weak_ptr< Atom_ > Hp( H );
		Ccp = Hp;
		TS_ASSERT( Ccp->weight_ == H.weight_ );
	}

};


