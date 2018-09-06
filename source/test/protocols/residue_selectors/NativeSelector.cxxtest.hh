// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/select/residue_selector/NativeSelector.cxxtest.hh
/// @brief test suite for core::select::residue_selector::NativeSelector
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Package headers
#include <protocols/residue_selectors/NativeSelector.hh>
#include <core/select/residue_selector/ResidueNameSelector.hh>
#include <protocols/simple_filters/ResidueCountFilter.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/memory.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// C++ headers
#include <string>

#ifdef SERIALIZATION
// Cereal headers
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

using namespace protocols::residue_selectors;
using namespace core::select::residue_selector;

static basic::Tracer TR( "test.core.select.NativeSelectorTests" );

class NativeSelectorTests : public CxxTest::TestSuite {

public:

	void setUp() {
		//devel_init();
		protocols_init_with_additional_options( "-native protocols/residue_selectors/two_alanines.pdb" );
	}

	void test_native_is_too_small(){
		protocols_init_with_additional_options( "-native protocols/residue_selectors/one_serine.pdb" );

		NativeSelectorOP selector = utility::pointer::make_shared< NativeSelector >();
		selector->set_residue_selector( utility::pointer::make_shared< ResidueNameSelector >( "SER", true ) );

		protocols::simple_filters::ResidueCountFilter filter;
		filter.residue_selector( selector );

		core::pose::PoseOP pose = core::import_pose::pose_from_file( "protocols/residue_selectors/two_serines.pdb" );
		auto const selection = selector->apply( * pose );
		TS_ASSERT_EQUALS( selection.size(), 2 );

		bool const val1 = selection[1];
		TS_ASSERT_EQUALS( val1, true );

		bool const val2 = selection[2];
		TS_ASSERT_EQUALS( val2, false );

		TS_ASSERT_EQUALS( core::Size( filter.report_sm( * pose ) ), 1 );
	}

	void test_native_is_too_large(){
		protocols_init_with_additional_options( "-native protocols/residue_selectors/two_serines.pdb" );

		NativeSelectorOP selector = utility::pointer::make_shared< NativeSelector >();
		selector->set_residue_selector( utility::pointer::make_shared< ResidueNameSelector >( "SER", true ) );

		protocols::simple_filters::ResidueCountFilter filter;
		filter.residue_selector( selector );

		core::pose::PoseOP pose = core::import_pose::pose_from_file( "protocols/residue_selectors/one_serine.pdb" );
		auto const selection = selector->apply( * pose );
		TS_ASSERT_EQUALS( selection.size(), 1 );
		bool const val1 = selection[1];
		TS_ASSERT_EQUALS( val1, true );
		TS_ASSERT_EQUALS( core::Size( filter.report_sm( * pose ) ), 1 );
	}

	void test_simple1(){
		NativeSelectorOP selector = utility::pointer::make_shared< NativeSelector >();
		selector->set_residue_selector( utility::pointer::make_shared< ResidueNameSelector >( "ALA", true ) );

		protocols::simple_filters::ResidueCountFilter filter;
		filter.residue_selector( selector );

		core::pose::PoseOP pose = core::import_pose::pose_from_file( "protocols/residue_selectors/two_serines.pdb" );
		TS_ASSERT_EQUALS( core::Size( filter.report_sm( * pose ) ), 2 );
	}

	void test_simple2(){
		NativeSelectorOP selector = utility::pointer::make_shared< NativeSelector >();
		selector->set_residue_selector( utility::pointer::make_shared< ResidueNameSelector >( "SER", true ) );

		protocols::simple_filters::ResidueCountFilter filter;
		filter.residue_selector( selector );

		core::pose::PoseOP pose = core::import_pose::pose_from_file( "protocols/residue_selectors/two_serines.pdb" );
		TS_ASSERT_EQUALS( core::Size( filter.report_sm( * pose ) ), 0 );
	}

	void test_serialize_simple1(){
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		NativeSelectorOP selector = utility::pointer::make_shared< NativeSelector >();
		selector->set_residue_selector( utility::pointer::make_shared< ResidueNameSelector >( "ALA", true ) );

		core::pose::PoseOP pose = core::import_pose::pose_from_file( "protocols/residue_selectors/two_serines.pdb" );

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( selector );
		}

		NativeSelectorOP instance2;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}

		TS_ASSERT( utility::pointer::dynamic_pointer_cast< NativeSelector > ( instance2 ) );

		protocols::simple_filters::ResidueCountFilter filter;
		filter.residue_selector( instance2 );
		TS_ASSERT_EQUALS( core::Size( filter.report_sm( * pose ) ), 2 );
#endif
	}

	void test_serialize_simple2(){
		TS_ASSERT( true ); // for non-serialization builds
#ifdef SERIALIZATION
		NativeSelectorOP selector = utility::pointer::make_shared< NativeSelector >();
		selector->set_residue_selector( utility::pointer::make_shared< ResidueNameSelector >( "SER", true ) );

		core::pose::PoseOP pose = core::import_pose::pose_from_file( "protocols/residue_selectors/two_serines.pdb" );

		std::ostringstream oss;
		{
			cereal::BinaryOutputArchive arc( oss );
			arc( selector );
		}

		NativeSelectorOP instance2;
		std::istringstream iss( oss.str() );
		{
			cereal::BinaryInputArchive arc( iss );
			arc( instance2 );
		}

		TS_ASSERT( utility::pointer::dynamic_pointer_cast< NativeSelector > ( instance2 ) );

		protocols::simple_filters::ResidueCountFilter filter;
		filter.residue_selector( instance2 );
		TS_ASSERT_EQUALS( core::Size( filter.report_sm( * pose ) ), 0 );
#endif
	}

};
