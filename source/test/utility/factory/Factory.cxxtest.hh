// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/factory/Factory.cxxtest.hh
/// @brief  utility::factory::Factory test suite
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Ron Jacak (ron.jacak@gmail.com)


// Package headers
#include <cxxtest/TestSuite.h>
#include <utility/factory/Factory.hh>
#include <utility/factory/Registrant.hh>

// Project headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <string>


/// @brief Product interface class
class Product : public utility::pointer::ReferenceCount {

public: // Types

	typedef  utility::pointer::shared_ptr< Product >  ProductOP;
	typedef  utility::factory::Factory< Product >  Factory;
	typedef  std::string const  FactoryKey;
	typedef  std::string const *  FactoryKeyP;
	typedef  FactoryKey & (*FactoryKeyFxn)(); // Key function
	typedef  ProductOP  FactoryProductP;
	typedef  FactoryProductP (*FactoryCreate)(); // Creation function: 0 arg case
	typedef  int  FactoryArg; // Not used here

protected: // Creation

	/// @brief Constructor
	inline
	Product()
	{}

public: // Creation

	/// @brief Destructor
	inline
	virtual
	~Product()
	{}

public: // Methods

	/// @brief Name
	virtual
	std::string const &
	name() const = 0;

}; // Product


/// @brief Concrete Product A
class ProductA : public Product {

public: // Types

	typedef  Product  Super;

public: // Creation

	/// @brief Default constructor
	inline
	ProductA()
	{}

	/// @brief Copy constructor
	inline
	ProductA( ProductA const & a ) : Super( a )
	{}

	/// @brief Factory creation
	inline
	static
	ProductOP
	factory_create() {
		return ProductOP( new ProductA() );
	}

	/// @brief Destructor
	inline
	virtual
	~ProductA()
	{}


public: // Methods

	inline
	std::string const &
	name() const {
		static std::string const name_( "ProductA" );
		return name_;
	}

private: // Static fields

	/// @brief Factory registrant: Registers this class with Factory
	static utility::factory::Registrant< Product > const factory_registrant_;

}; // ProductA


// Static field definitions
utility::factory::Registrant< Product > const ProductA::factory_registrant_( "ProductA", ProductA::factory_create );


/// @brief Concrete Product B
class ProductB : public Product {

public: // Types

	typedef  Product  Super;

public: // Creation

	/// @brief Default constructor
	inline
	ProductB()
	{}

	/// @brief Copy constructor
	inline
	ProductB( ProductB const & b ) :
		Super( b )
	{}

	/// @brief Factory creation
	inline
	static
	ProductOP
	factory_create()
	{
		return ProductOP( new ProductB() );
	}

	/// @brief Destructor
	inline
	virtual
	~ProductB()
	{}

public: // Methods

	inline
	std::string const &
	name() const {
		static std::string const name_( "ProductB" );
		return name_;
	}

private: // Static fields

	/// @brief Factory registrant: Registers this class with Factory
	static utility::factory::Registrant< Product > const factory_registrant_;

}; // ProductB


// Static field definitions
utility::factory::Registrant< Product > const ProductB::factory_registrant_( "ProductB", ProductB::factory_create );


class FactoryTests : public CxxTest::TestSuite {

	public:

	/// @brief General tests
	void test_Factory_general() {

		TS_ASSERT( Product::Factory::create( "ProductA" )->name() == "ProductA" );
		TS_ASSERT( Product::Factory::create( "ProductB" )->name() == "ProductB" );

	}

};

