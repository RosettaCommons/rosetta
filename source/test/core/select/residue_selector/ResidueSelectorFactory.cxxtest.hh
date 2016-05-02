// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/select/residue_selector/ResidueSelectorFactory.cxxtest.hh
/// @brief test suite for basic::resource_manager::ResourceOptionsFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreator.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <string>

static THREAD_LOCAL basic::Tracer TR("core.select.residue_selector.ResidueSelectorFactory.cxxtest.hh");

using namespace core::select::residue_selector;

class DummyResidueSelector : public ResidueSelector {
public:
	DummyResidueSelector() {}

	DummyResidueSelector( DummyResidueSelector const & ) {}

	ResidueSelectorOP clone() const { return ResidueSelectorOP( new DummyResidueSelector(*this) ); }

	virtual
	ResidueSubset apply( core::pose::Pose const & ) const
	{
		return ResidueSubset( 10, true );
	}

	virtual std::string get_name() const { return "DummyResidueSelector"; }

};

class DummyResidueSelectorCreator : public ResidueSelectorCreator {
public:
	virtual
	ResidueSelectorOP
	create_residue_selector() const {
		return ResidueSelectorOP( new DummyResidueSelector );
	}

	virtual
	std::string keyname() const { return "DummyResidueSelector"; }

	virtual void provide_selector_xsd( utility::tag::XMLSchemaDefinition & ) const {}
};


class ResidueSelectorFactoryTests : public CxxTest::TestSuite {

public:

	void setUp() { core_init(); }

	// @brief make sure that when we register a residue selector, we can later get it back
	void test_register_one_creator_with_ResidueSelectorFactory() {
		ResidueSelectorFactory * factory = ResidueSelectorFactory::get_instance();
		factory->set_throw_on_double_registration();
		factory->factory_register( ResidueSelectorCreatorOP( new DummyResidueSelectorCreator ) );
		basic::datacache::DataMap dm;
		ResidueSelectorOP selector = factory->new_residue_selector( "DummyResidueSelector", 0, dm );
		TS_ASSERT( selector.get() ); // make sure we got back a non-null pointer
		DummyResidueSelector * dselector = dynamic_cast< DummyResidueSelector * > ( selector.get() );
		TS_ASSERT( dselector ); // make sure we got back the right residue selector kind
	}

	/// @brief make sure that if a ResidueSelector with the name of an already-registered
	/// ResidueSelector gets registered with the factory, that the factory throws an exception
	void test_register_one_creator_twice_with_ResidueSelectorFactory() {
		try {
			ResidueSelectorFactory * factory = ResidueSelectorFactory::get_instance();
			factory->set_throw_on_double_registration();
			factory->factory_register( ResidueSelectorCreatorOP( new DummyResidueSelectorCreator ) );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_err_msg = "Factory Name Conflict: Two or more ResidueSelectorCreators registered with the name DummyResidueSelector";
			TS_ASSERT( expected_err_msg == e.msg() );

		}

	}

	// @brief make sure that when we register a residue selector, we can later get it back
	void test_throw_on_unregistered_selector_name_ResidueSelectorFactory() {
		try {
			ResidueSelectorFactory * factory = ResidueSelectorFactory::get_instance();
			basic::datacache::DataMap dm;
			ResidueSelectorOP selector = factory->new_residue_selector( "DummyResidueSelector2", 0, dm );
			TS_ASSERT( false );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::string expected_err_msg = "No ResidueSelectorCreator with the name 'DummyResidueSelector2' has been registered with the ResidueSelectorFactory";
			TS_ASSERT( expected_err_msg == e.msg() );
		}
	}


	void test_output_residue_selector_xsds() {
		utility::tag::XMLSchemaDefinition xsd;
		ResidueSelectorFactory * factory = ResidueSelectorFactory::get_instance();
		factory->define_residue_selector_xml_schema( xsd );
		TR << "XSD: " << std::endl;
		TR << xsd.full_definition() << std::endl;
	}
};
