// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/pointer/deep_copy.cxxtest.hhi
/// @brief  Unit tests for utility::pointer::DeepCopyOP
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Headers
#include <iostream> // Put this first to allow debug print-outs in project headers
#include <set>

// Project headers
#include <cxxtest/TestSuite.h>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/deep_copy.hh>
#include <utility/pointer/memory.hh>

//#include <cxxabi.h>

// --- set up some helper classes/functions for these tests
// --- (these go in a special namespace to be used from this file only!)

namespace deep_copy_test {

struct CloneCounts {
	int clone_counter = 0;
	int copy_counter = 0;
};

typedef utility::pointer::shared_ptr< CloneCounts > CloneCountsOP;

class ClonableBase;
typedef utility::pointer::shared_ptr< ClonableBase > ClonableBaseOP;
typedef utility::pointer::shared_ptr< ClonableBase const > ClonableBaseCOP;

class ClonableBase {
public:

	virtual ~ClonableBase() = default;

	ClonableBase( int ident = 0 ):
		ident_( ident ),
		counts_( new CloneCounts )
	{}

	ClonableBase( ClonableBase const & s ):
		ident_( s.ident_ ),
		counts_( s.counts_ ) // shared
	{
		++counts_->copy_counter; // increment for both.
	}

	virtual
	ClonableBaseOP clone() const {
		++counts_->clone_counter;
		return ClonableBaseOP( new ClonableBase( *this ) );
	}

public: // Data members

	int ident_ = 0;
	CloneCountsOP counts_;

};

ClonableBaseOP deep_copy( ClonableBase const & source ) { return source.clone(); }

class ClonableOne: public ClonableBase {
public:
	ClonableOne( int ident = 0 ):
		ClonableBase( ident + 10 )
	{}

	ClonableOne( ClonableOne const & other ):
		ClonableBase( other )
	{}

	ClonableBaseOP clone() const override {
		++counts_->clone_counter;
		return ClonableBaseOP( new ClonableOne( *this ) );
	}
};

class ClonableTwo: public ClonableBase {
public:
	ClonableTwo( int ident = 0 ):
		ClonableBase( ident + 20 )
	{}

	ClonableTwo( ClonableTwo const & other ):
		ClonableBase( other )
	{}

	ClonableBaseOP clone() const override {
		++counts_->clone_counter;
		return ClonableBaseOP( new ClonableTwo( *this ) );
	}
};

class NonClonable;
typedef utility::pointer::shared_ptr< NonClonable > NonClonableOP;
typedef utility::pointer::shared_ptr< NonClonable const > NonClonableCOP;

class NonClonable {
public:
	NonClonable( int ident = 0 ):
		ident_( ident ),
		counts_( new CloneCounts )
	{}

	NonClonable( NonClonable const & s ):
		ident_( s.ident_ ),
		counts_( s.counts_ ) // shared
	{
		++counts_->copy_counter; // increment for both.
	}

public: // Data members

	int ident_ = 0;
	CloneCountsOP counts_;

};

NonClonableOP deep_copy( NonClonable const & source ) {
	return utility::pointer::make_shared< NonClonable >( source );
}

class DeepCopied {
public:

	DeepCopied():
		shallow_c_( new ClonableOne ),
		shallow_nc_( new NonClonable ),
		deep_c_( new ClonableBase ),
		deep_cc_( new ClonableBase ),
		deep_c1_( new ClonableOne ),
		deep_nc_( new NonClonable ),
		deep_null_( nullptr )
	{}

	DeepCopied( DeepCopied const & ) = default;

public:

	ClonableBaseOP shallow_c_;
	NonClonableOP shallow_nc_;

	utility::pointer::DeepCopyOP< ClonableBase >  deep_c_;
	utility::pointer::DeepCopyOP< ClonableBase const >  deep_cc_;
	utility::pointer::DeepCopyOP< ClonableBase >  deep_c1_;
	//utility::pointer::DeepCopyOP< ClonableOne > deep_c1_; // Doesn't work - clone() doesn't return the correct type.
	utility::pointer::DeepCopyOP< NonClonable > deep_nc_;
	utility::pointer::DeepCopyOP< ClonableBase > deep_null_;
};

} // namespace deep_copy_test

using namespace deep_copy_test;

class DeepCopyOPTests : public CxxTest::TestSuite {


	public:

	/// @brief Explicit owning pointer
	void test_copying() {
		DeepCopied source;

		// Check that everything is set up.
		TS_ASSERT_EQUALS( source.shallow_c_->counts_->copy_counter, 0 );
		TS_ASSERT_EQUALS( source.shallow_nc_->counts_->copy_counter, 0 );
		TS_ASSERT_EQUALS( source.deep_c_->counts_->copy_counter, 0 );
		TS_ASSERT_EQUALS( source.deep_cc_->counts_->copy_counter, 0 );
		TS_ASSERT_EQUALS( source.deep_c1_->counts_->copy_counter, 0 );
		TS_ASSERT_EQUALS( source.deep_nc_->counts_->copy_counter, 0 );
		TS_ASSERT_EQUALS( source.shallow_c_->counts_->clone_counter, 0 );
		TS_ASSERT_EQUALS( source.shallow_nc_->counts_->clone_counter, 0 );
		TS_ASSERT_EQUALS( source.deep_c_->counts_->clone_counter, 0 );
		TS_ASSERT_EQUALS( source.deep_c1_->counts_->clone_counter, 0 );
		TS_ASSERT_EQUALS( source.deep_nc_->counts_->clone_counter, 0 );

		DeepCopied destination( source );

		// Check for shallow/deep copies with pointer identities.
		TS_ASSERT( source.shallow_c_.get() == destination.shallow_c_.get() );
		TS_ASSERT( source.shallow_nc_.get() == destination.shallow_nc_.get() );
		TS_ASSERT( source.deep_c_.get() != destination.deep_c_.get() );
		TS_ASSERT( source.deep_cc_.get() != destination.deep_cc_.get() );
		TS_ASSERT( source.deep_c1_.get() != destination.deep_c1_.get() );
		TS_ASSERT( source.deep_nc_.get() != destination.deep_nc_.get() );
		TS_ASSERT( destination.deep_null_.get() == nullptr ); // Actually, we're mainly concerned we don't get a nullptr error with usage

		// Check that we didn't accidentally smash the clonable types to the base class.
		TS_ASSERT( !dynamic_cast< ClonableOne * >( source.deep_c_.get() ) );
		TS_ASSERT( !dynamic_cast< ClonableOne * >( destination.deep_c_.get() ) );
		TS_ASSERT( !dynamic_cast< ClonableOne const * >( source.deep_cc_.get() ) );
		TS_ASSERT( !dynamic_cast< ClonableOne const * >( destination.deep_cc_.get() ) );
		TS_ASSERT( dynamic_cast< ClonableOne * >( source.deep_c1_.get() ) );
		TS_ASSERT( dynamic_cast< ClonableOne * >( destination.deep_c1_.get() ) );

		// Check the runtime types of the sub-items
		// This indirection is needed, as Clang warns about embedding the smart pointer dereferencing in the typeid()
		auto & sdc( *source.deep_c_.get() );
		auto & sdcc( *source.deep_cc_.get() );
		auto & sdc1( *source.deep_c1_.get() );
		auto & sdnc( *source.deep_nc_.get() );
		auto & ddc( *destination.deep_c_.get() );
		auto & ddcc( *destination.deep_cc_.get() );
		auto & ddc1( *destination.deep_c1_.get() );
		auto & ddnc( *destination.deep_nc_.get() );
		TS_ASSERT_EQUALS( typeid( sdc ).name(), typeid( ddc ).name() );
		TS_ASSERT_EQUALS( typeid( sdcc ).name(), typeid( ddcc ).name() );
		TS_ASSERT_EQUALS( typeid( sdc1 ).name(), typeid( ddc1 ).name() );
		TS_ASSERT_DIFFERS( typeid( ddc ).name(), typeid( ddc1 ).name() ); // First should be ClonableBase, the second ClonableOne (or their mangled equivalents)
		TS_ASSERT_EQUALS( typeid( sdnc ).name(), typeid( ddnc ).name() );

		//int status = 0;
		//std::cout << "C: " << abi::__cxa_demangle( typeid( sdc ).name(), NULL, NULL, &status ) << std::endl;
		//std::cout << "CC: " << abi::__cxa_demangle( typeid( sdcc ).name(), NULL, NULL, &status ) << std::endl;
		//std::cout << "C1: " << abi::__cxa_demangle( typeid( sdc1 ).name(), NULL, NULL, &status ) << std::endl;
		//std::cout << "NC: " << abi::__cxa_demangle( typeid( sdnc ).name(), NULL, NULL, &status ) << std::endl;

		// Check the appropriate functions were called.
		TS_ASSERT_EQUALS( destination.shallow_c_->counts_->copy_counter, 0 ); // Just the pointer was copied
		TS_ASSERT_EQUALS( destination.shallow_nc_->counts_->copy_counter, 0 ); // Just the pointer was copied
		TS_ASSERT_EQUALS( destination.deep_c_->counts_->copy_counter, 1 );
		TS_ASSERT_EQUALS( destination.deep_cc_->counts_->copy_counter, 1 );
		TS_ASSERT_EQUALS( destination.deep_c1_->counts_->copy_counter, 1 );
		TS_ASSERT_EQUALS( destination.deep_nc_->counts_->copy_counter, 1 );
		TS_ASSERT_EQUALS( destination.shallow_c_->counts_->clone_counter, 0 ); // Just the pointer was copied
		TS_ASSERT_EQUALS( destination.shallow_nc_->counts_->clone_counter, 0 ); // Just the pointer was copied
		TS_ASSERT_EQUALS( destination.deep_c_->counts_->clone_counter, 1 ); // IMPORTANT -- we need to make sure the pointer was cloned.
		TS_ASSERT_EQUALS( destination.deep_cc_->counts_->clone_counter, 1 ); // IMPORTANT -- we need to make sure the pointer was cloned.
		TS_ASSERT_EQUALS( destination.deep_c1_->counts_->clone_counter, 1 ); // IMPORTANT -- we need to make sure the pointer was cloned.
		TS_ASSERT_EQUALS( destination.deep_nc_->counts_->clone_counter, 0 );

		//std::cout << "Done running deep copy tests." << std::endl;
	}

	bool implicit_convert_ClonableBaseOP(ClonableBaseOP val) {
		return val != nullptr;
	}
	bool implicit_convert_ClonableBaseCOP(ClonableBaseCOP val) {
		return val != nullptr;
	}

	void test_conversion() {
		DeepCopied source;

		// Test we can convert from a DeepCopyOP to a sensible value.
		// (These mainly are checking that we don't get compiler errors with the usage.)
		TS_ASSERT( implicit_convert_ClonableBaseCOP( source.deep_c_ ) );
		TS_ASSERT( implicit_convert_ClonableBaseCOP( source.deep_cc_ ) );
		TS_ASSERT( implicit_convert_ClonableBaseCOP( source.deep_c1_ ) );
		TS_ASSERT( ! implicit_convert_ClonableBaseCOP( source.deep_null_ ) );

		NonClonableCOP non_clone_const( source.deep_nc_ );
		TS_ASSERT( non_clone_const != nullptr );

		TS_ASSERT( implicit_convert_ClonableBaseOP( source.deep_c_ ) );
		//TS_ASSERT( implicit_convert_ClonableBaseOP( source.deep_cc_ ) ); // This is, as should be expected, a compiler error.
		TS_ASSERT( implicit_convert_ClonableBaseOP( source.deep_c1_ ) );
		TS_ASSERT( ! implicit_convert_ClonableBaseOP( source.deep_null_ ) );

		NonClonableOP non_clone( source.deep_nc_ );
		TS_ASSERT( non_clone != nullptr );

		// Test bool conversion
		TS_ASSERT( !!( source.deep_c_ ) ); // Two NOT is getting the truth value
		TS_ASSERT( !!( source.deep_nc_ ) );
		TS_ASSERT( !( source.deep_null_ ) );

		DeepCopied destination( source );

		// Test that we can convert from OP/COP

		// equals
		ClonableBaseOP my_cb( new ClonableOne );
		ClonableBaseCOP my_const_cb( new ClonableOne );

		destination.deep_c_ = my_cb;
		destination.deep_cc_ = my_cb;

		TS_ASSERT_EQUALS( destination.deep_c_.get(), my_cb.get() );
		TS_ASSERT_EQUALS( destination.deep_cc_.get(), my_cb.get() );

		//destination.deep_c_ = my_const_cb; // Will not compile
		destination.deep_cc_ = my_const_cb;

		TS_ASSERT_EQUALS( destination.deep_cc_.get(), my_const_cb.get() );

		destination.deep_c_ = nullptr;
		destination.deep_cc_ = nullptr;

		TS_ASSERT_EQUALS( destination.deep_c_.get(), nullptr );
		TS_ASSERT_EQUALS( destination.deep_cc_.get(), nullptr );
	}
};


