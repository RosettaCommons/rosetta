// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file basic/resource_manager/ResourceOptionsFactory.cxxtest.hh
/// @brief test suite for basic::resource_manager::ResourceOptionsFactory
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <basic/datacache/HierarchicalDataMap.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <iostream>
#include <string>

struct Widget;
typedef utility::pointer::shared_ptr< Widget > WidgetOP;

struct Widget : public utility::pointer::ReferenceCount {
	Widget(Size id=0) : id_(id) {}
	bool operator == (WidgetOP other) { return id_ == other->id_; }
	Size id_;
};

using namespace std;
using basic::datacache::HierarchicalDataMap;
using basic::datacache::HierarchicalDataMapOP;
using utility::excn::EXCN_Msg_Exception;

class HierarchicalDataMapTests : public CxxTest::TestSuite {

public:

	WidgetOP widgets[3];
	HierarchicalDataMapOP parent, child;

	void setUp() {
		parent = HierarchicalDataMapOP( new HierarchicalDataMap );
		child = HierarchicalDataMapOP( new HierarchicalDataMap );
		child->set_parent(parent);

		widgets[0] = WidgetOP( new Widget(0) );
		widgets[1] = WidgetOP( new Widget(1) );
		widgets[2] = WidgetOP( new Widget(2) );
	}

	void teardown() {
		parent = HierarchicalDataMapOP();
		child = HierarchicalDataMapOP();
	}

	void test_successful_lookups() {
		child->set("widget", "0", widgets[0]);  // Decoy widget.
		child->set("widget", "1", widgets[1]);

		TS_ASSERT_EQUALS(child->get<WidgetOP>("widget", "1"), widgets[1]);
		TS_ASSERT_EQUALS(child->get<WidgetOP>("widget", "1", widgets[2]), widgets[1]);
		TS_ASSERT_EQUALS(child->get_or_null<WidgetOP>("widget", "1"), widgets[1]);
	}

	void test_failed_lookups() {
		set_throw_on_next_assertion_failure();
		TS_ASSERT_THROWS(child->get<WidgetOP>("widget", "1"), EXCN_Msg_Exception);
		TS_ASSERT_EQUALS(child->get<WidgetOP>("widget", "1", widgets[2]), widgets[2]);
		TS_ASSERT_EQUALS(child->get_or_null<WidgetOP>("widget", "1"), WidgetOP());
	}

	void test_optional_name_lookups() {
		child->set("widget", "", widgets[0]);
		child->set("widget", "named", widgets[1]);

		TS_ASSERT_EQUALS(child->get<WidgetOP>("widget", ""), widgets[0]);
		TS_ASSERT_EQUALS(child->get<WidgetOP>("widget", "unnamed"), widgets[0]);
		TS_ASSERT_EQUALS(child->get<WidgetOP>("widget", "named"), widgets[1]);
	}

	void test_hierarchical_lookups() {
		parent->set("widget", "", widgets[0]);
		TS_ASSERT_EQUALS(parent->get<WidgetOP>("widget", ""), widgets[0]);
		TS_ASSERT_EQUALS(child->get<WidgetOP>("widget", ""), widgets[0]);

		child->set("widget", "", widgets[1]);
		TS_ASSERT_EQUALS(parent->get<WidgetOP>("widget", ""), widgets[0]);
		TS_ASSERT_EQUALS(child->get<WidgetOP>("widget", ""), widgets[1]);
	}

};
