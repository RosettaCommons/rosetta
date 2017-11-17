// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/ConstraintGeneratorFactory.cxxtest.hh
/// @brief test suite for protocols::constraint_generator::ConstraintGeneratorFactory
/// @author Tom Linsky (tlinsky at uw dot edu)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Protocol headers
#include <protocols/constraint_generator/ConstraintGenerator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorCreator.hh>
#include <protocols/constraint_generator/ConstraintGeneratorFactory.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/excn/Exceptions.hh>

// Boost headers
#include <boost/assign.hpp>

// C++ headers

using namespace protocols::constraint_generator;

class DummyConstraintGeneratorCreator;
class DummyConstraintGenerator : public ConstraintGenerator {
public:
	DummyConstraintGenerator():
		ConstraintGenerator( "DummyConstraintGenerator" ) {}

	ConstraintGeneratorOP clone() const { return ConstraintGeneratorOP( new DummyConstraintGenerator(*this) ); }

	virtual core::scoring::constraints::ConstraintCOPs
	apply( core::pose::Pose const & ) const
	{
		return boost::assign::list_of
			(core::scoring::constraints::ConstraintCOP())
			(core::scoring::constraints::ConstraintCOP())
			(core::scoring::constraints::ConstraintCOP());
	}

protected:
	virtual void
	parse_tag( utility::tag::TagCOP , basic::datacache::DataMap & )
	{}

};

class DummyConstraintGeneratorCreator : public ConstraintGeneratorCreator {
public:
	virtual ConstraintGeneratorOP
	create_constraint_generator() const {
		return ConstraintGeneratorOP( new DummyConstraintGenerator );
	}

	virtual
	std::string keyname() const { return constraint_generator_name(); }

	static
	std::string constraint_generator_name() { return "DummyConstraintGenerator"; }

	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & ) const {}
};


class ConstraintGeneratorFactoryTests : public CxxTest::TestSuite {

public:

	void setUp() {
	}

	// @brief make sure that when we register a constraint generator creator, we can get a valid constraint generator from the factory
	void test_register_one_creator_with_ConstraintGeneratorFactory() {
		std::stringstream tag_ss( "<DummyConstraintGenerator name=dummy_unit1 />" );
		utility::tag::TagCOP tag( utility::tag::Tag::create( tag_ss ) );
		ConstraintGeneratorFactory * factory = ConstraintGeneratorFactory::get_instance();
		factory->factory_register( ConstraintGeneratorCreatorOP( new DummyConstraintGeneratorCreator ) );
		basic::datacache::DataMap dm;
		ConstraintGeneratorOP generator = factory->new_constraint_generator( "DummyConstraintGenerator", tag, dm );
		TS_ASSERT( generator.get() ); // make sure we got back a non-null pointer
		DummyConstraintGenerator * dgenerator = dynamic_cast< DummyConstraintGenerator * > ( generator.get() );
		TS_ASSERT( dgenerator ); // make sure we got back the right residue selector kind

		// make sure the generator data is properly set
		TS_ASSERT_EQUALS( generator->class_name(), "DummyConstraintGenerator" );
		TS_ASSERT_EQUALS( generator->id(), "dummy_unit1" );
	}

	/// @brief make sure that the factory throws an exception if a two ConstraintGeneratorCreators with
	/// the same name are registered
	void test_register_one_creator_twice_with_ConstraintGeneratorFactory() {
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS_ANYTHING(
			ConstraintGeneratorFactory::get_instance()->factory_register( ConstraintGeneratorCreatorOP( new DummyConstraintGeneratorCreator )   )
		);
	}


	// @brief make sure that when we try to create an unregistered constraint generator, an exception is thrown
	void test_throw_on_unregistered_generator_name_ConstraintGeneratorFactory() {
		basic::datacache::DataMap dm;
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS( ConstraintGeneratorFactory::get_instance()->new_constraint_generator( "DummyConstraintGenerator2", 0, dm ),
			utility::excn::Exception );
	}

};
