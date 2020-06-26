// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/moves/PyMOLMoverTest.cxxtest.hh
/// @brief  Test some expected behavior of the PyMOLMover and PyMOLObserver
/// @author Kyle Barlow (kb@kylebarlow.com)
/// @details These unit tests define expected behavior relating to how a PyMOLObserver
///          is stored inside the Pose's datacache ObserverCache

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers


// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/datacache/ObserverCache.hh>
#include <core/pose/datacache/CacheableObserverType.hh>

// Protocol Headers
#include <basic/Tracer.hh>

#include <protocols/moves/PyMOLMover.hh>

static basic::Tracer TR("PyMOLMoverTest");


class PyMOLMoverTest : public CxxTest::TestSuite {
	//Define Variables

public:

	core::pose::PoseOP original_pose;
	void setUp(){
		original_pose = utility::pointer::make_shared< core::pose::Pose >();
	}

	void tearDown(){
		original_pose = utility::pointer::make_shared< core::pose::Pose >();
	}

	void test_cmdline_01(){
		core_init_with_additional_options("-run:PyMOLMover:address 192.1.1.3 -run:PyMOLMover:port 9998");
		protocols::moves::PyMOLObserverOP original_observer = protocols::moves::AddPyMOLObserver( *original_pose, true, 0 );
		TS_ASSERT_EQUALS(original_observer->pymol().get_link().get_address(), "192.1.1.3");
		TS_ASSERT_EQUALS(original_observer->pymol().get_link().get_port(), 9998);
		TS_ASSERT_EQUALS(original_observer->pymol().get_link().get_max_packet_size(), 1500);
	}

	void test_link(){
		// WARNING!
		// 8192-512-2 is NOT a scientific number and was chosen somewhat arbitrarily!!!
		// these tests are NOT testing that the number should be 8192-512-2, but rather
		// that the max packet size is set differently when the address is "127.0.0.1"
		// vs other addresses
		core_init();
		std::string const addr1("127.0.0.1");
		unsigned int const port1(65000);
		protocols::moves::UDPSocketClient expected_init(addr1, port1, 8192-512-2);
		TR << expected_init << std::endl;
		TR << expected_init.get_address() << std::endl;
		TR << expected_init.get_port() << std::endl;
		std::string const addr2("0.0.0.0");
		unsigned int const port2(62345);
		protocols::moves::UDPSocketClient expected_init2(addr2, port2, 1555);
		TS_ASSERT_EQUALS(expected_init.get_max_packet_size(), 8192-512-2);
		TS_ASSERT_EQUALS(expected_init2.get_max_packet_size(), 1555);

		TS_ASSERT_EQUALS(expected_init.get_address(), addr1);
		TS_ASSERT_EQUALS(expected_init2.get_address(), addr2);

		TS_ASSERT_EQUALS(expected_init.get_port(), port1);
		TS_ASSERT_EQUALS(expected_init2.get_port(), port2);

		TS_ASSERT_DIFFERS(expected_init.get_address(), expected_init2.get_address());
		TS_ASSERT_DIFFERS(expected_init.get_port(), expected_init2.get_port());
		TS_ASSERT_DIFFERS(expected_init.get_max_packet_size(), expected_init2.get_max_packet_size());

		protocols::moves::PyMOLObserverOP original_observer = protocols::moves::AddPyMOLObserver( *original_pose, true, 0 );
		TS_ASSERT_EQUALS(original_observer->pymol().get_link().get_address(), expected_init.get_address());
		TS_ASSERT_EQUALS(original_observer->pymol().get_link().get_port(), expected_init.get_port());
		TS_ASSERT_EQUALS(original_observer->pymol().get_link().get_max_packet_size(), expected_init.get_max_packet_size());
		original_observer->pymol().set_link(expected_init2);
		TS_ASSERT_EQUALS(original_observer->pymol().get_link().get_address(), expected_init2.get_address());
		TS_ASSERT_EQUALS(original_observer->pymol().get_link().get_port(), expected_init2.get_port());
		TS_ASSERT_EQUALS(original_observer->pymol().get_link().get_max_packet_size(), expected_init2.get_max_packet_size());
	}

	void test_pose_clone_has_no_observer(){
		core_init();
		// Test: A clone of a pose with observer should not clone the observer attached to the original pose
		//       This is an exception to the normal deep copy process of cloning, because it prevents
		//       multiple PyMOLObservers from being created non-interactively (which would be confusing
		//       to the end user, as PyMOLObservers are often created interactively).
		protocols::moves::PyMOLObserverOP original_observer = protocols::moves::AddPyMOLObserver( *original_pose, true, 0 );
		TS_ASSERT( original_pose->observer_cache().has( core::pose::datacache::CacheableObserverType::PYMOL_OBSERVER ) );
		core::pose::Pose copy_pose = *original_pose;
		TS_ASSERT( !copy_pose.observer_cache().has( core::pose::datacache::CacheableObserverType::PYMOL_OBSERVER ) );
	}

	void test_replaced_pose_keeps_observer(){
		core_init();
		// Test: Copying a pose into original pose should preserve the observer in the original pose
		protocols::moves::PyMOLObserverOP original_observer = protocols::moves::AddPyMOLObserver( *original_pose, true, 0 );
		TS_ASSERT( original_pose->observer_cache().has( core::pose::datacache::CacheableObserverType::PYMOL_OBSERVER ) );
		core::pose::Pose new_pose = core::pose::Pose();
		*original_pose = new_pose;
		TS_ASSERT( original_pose->observer_cache().has( core::pose::datacache::CacheableObserverType::PYMOL_OBSERVER ) );
	}

	void test_pose_has_observer(){
		core_init();
		// Test: Adding a PyMOLObserver to a Pose results in the Pose having a PyMOLObserver
		protocols::moves::PyMOLObserverOP original_observer = protocols::moves::AddPyMOLObserver( *original_pose, true, 0 );
		TS_ASSERT( original_pose->observer_cache().has( core::pose::datacache::CacheableObserverType::PYMOL_OBSERVER ) );
	}

	void test_new_pose_has_no_observer(){
		core_init();
		// Test: A new, fresh, PyMOLObserver doesn't start with a PyMOLObserver
		TS_ASSERT( !original_pose->observer_cache().has( core::pose::datacache::CacheableObserverType::PYMOL_OBSERVER ) );
	}

};
