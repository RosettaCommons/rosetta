// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/resource_manager/ResourceOptionsFactory.cxxtest.hh
/// @brief test suite for basic::resource_manager::ResourceOptionsFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <basic/resource_manager/ResourceOptions.hh>
#include <basic/resource_manager/ResourceOptionsCreator.hh>
#include <basic/resource_manager/ResourceOptionsFactory.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>

//#include <basic/resource_manager/ResourceOptions

using namespace basic::resource_manager;

class DummyResourceOptions : public ResourceOptions {
public:
	DummyResourceOptions() : somevar_( 1 ) {}
	virtual
	void
	parse_my_tag(
		utility::tag::TagCOP tag
	)
	{
		somevar_ = tag->getOption< int >( "somevar", 1 );
	}

	virtual
	std::string
	type() const { return "DummyResourceOptions"; }

public:
	int somevar_;

};

class DummyResourceOptionsCreator : public ResourceOptionsCreator {
public:
	virtual std::string options_type() const { return "DummyResourceOptions"; }
	virtual ResourceOptionsOP create_options() const { return ResourceOptionsOP( new DummyResourceOptions ); }
};


class ResourceOptionsFactoryTests : public CxxTest::TestSuite {

public:

	void setUp() {
	}

	// @brief test default ctor
	void test_register_one_creator_with_ResourceOptionsFactory() {

		std::string dummy_opts( "<DummyResourceOptions somevar=5/>\n" );
		std::istringstream doptstream( dummy_opts );
		utility::tag::TagCOP tag = utility::tag::Tag::create( doptstream );

		ResourceOptionsFactory * factory = ResourceOptionsFactory::get_instance();
		factory->set_throw_on_double_registration();
		factory->factory_register( ResourceOptionsCreatorOP( new DummyResourceOptionsCreator ) );
		ResourceOptionsOP resource = factory->create_resource_options( "DummyResourceOptions", tag );
		TS_ASSERT( resource.get() ); // make sure we got back a non-null pointer
		DummyResourceOptions * dresource = dynamic_cast< DummyResourceOptions * > ( resource.get() );
		TS_ASSERT( dresource ); // make sure we got back the right resource kind
		TS_ASSERT( dresource->somevar_ == 5 );
	}
	void test_register_one_creator_twice_with_ResourceOptionsFactory() {
		try {
			ResourceOptionsFactory * factory = ResourceOptionsFactory::get_instance();
			factory->set_throw_on_double_registration();
			factory->factory_register( ResourceOptionsCreatorOP( new DummyResourceOptionsCreator ) );
			TS_ASSERT( false );
		} catch (utility::excn::Exception & e ) {
			std::string expected_error_message = "Double registration of a ResourceOptionsCreator in the ResourceOptionsFactory, named DummyResourceOptions. Are there two registrators for this options object, or have you chosen a previously assigned name to a new resource option?";
			TS_ASSERT( e.msg().find(expected_error_message) != std::string::npos );
			//std::cout << e.msg() << std::endl;
		}

	}

};
