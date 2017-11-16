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
		core_init();
		original_pose = core::pose::PoseOP( new core::pose::Pose() );
	}

	void tearDown(){
		original_pose = core::pose::PoseOP( new core::pose::Pose() );
	}

	void test_pose_clone_has_no_observer(){
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
		// Test: Copying a pose into original pose should preserve the observer in the original pose
		protocols::moves::PyMOLObserverOP original_observer = protocols::moves::AddPyMOLObserver( *original_pose, true, 0 );
		TS_ASSERT( original_pose->observer_cache().has( core::pose::datacache::CacheableObserverType::PYMOL_OBSERVER ) );
		core::pose::Pose new_pose = core::pose::Pose();
		*original_pose = new_pose;
		TS_ASSERT( original_pose->observer_cache().has( core::pose::datacache::CacheableObserverType::PYMOL_OBSERVER ) );
	}

	void test_pose_has_observer(){
		// Test: Adding a PyMOLObserver to a Pose results in the Pose having a PyMOLObserver
		protocols::moves::PyMOLObserverOP original_observer = protocols::moves::AddPyMOLObserver( *original_pose, true, 0 );
		TS_ASSERT( original_pose->observer_cache().has( core::pose::datacache::CacheableObserverType::PYMOL_OBSERVER ) );
	}

	void test_new_pose_has_no_observer(){
		// Test: A new, fresh, PyMOLObserver doesn't start with a PyMOLObserver
		TS_ASSERT( !original_pose->observer_cache().has( core::pose::datacache::CacheableObserverType::PYMOL_OBSERVER ) );
	}

};
