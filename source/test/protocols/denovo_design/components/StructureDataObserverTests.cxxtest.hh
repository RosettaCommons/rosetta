// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/denovo_design/components/StructureDataObserverTests.cxxtest.hh
/// @brief  Test suite for StructureDataObserver
/// @author Tom Linsky (tlinsky@uw.edu)

// Test headers
#include <test/UMoverTest.hh>
#include <test/UTracer.hh>
#include <cxxtest/TestSuite.h>

// Project Headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/components/StructureDataObserver.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>

// Core Headers
#include <core/io/pdb/build_pose_as_is.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>

// Protocol Headers
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("StructureDataObserverTests");

using namespace protocols::denovo_design::components;

class StructureDataObserverTests : public CxxTest::TestSuite {
	//Define Variables
private:
	StructureDataFactory const * factory;

public:
	void setUp()
	{
		core_init();
		factory = StructureDataFactory::get_instance();
	}

	void tearDown() {}

	void test_attachment()
	{
		core::pose::Pose pose;
		core::io::pdb::build_pose_from_pdb_as_is( pose, "protocols/denovo_design/components/helix15.pdb" );

		std::string const sd_id = "UnitTest";
		StructureData const orig_sd = *factory->create_from_pose( pose, sd_id );

		TS_ASSERT( !pose.observer_cache().has( core::pose::datacache::CacheableObserverType::STRUCTUREDATA_OBSERVER ) );
		TS_ASSERT( !factory->observer_attached( pose ) );

		factory->attach_observer( pose, orig_sd );

		// should be attached in a slot
		TS_ASSERT( pose.observer_cache().has( core::pose::datacache::CacheableObserverType::STRUCTUREDATA_OBSERVER ) );
		TS_ASSERT( factory->observer_attached( pose ) );

		// should be correct types
		core::pose::datacache::CacheableObserverOP cache_obs =
			pose.observer_cache().get_ptr( core::pose::datacache::CacheableObserverType::STRUCTUREDATA_OBSERVER );
		TS_ASSERT( cache_obs );
		StructureDataObserverOP sd_obs = utility::pointer::static_pointer_cast< StructureDataObserver >( cache_obs );
		TS_ASSERT( sd_obs );

		sd_obs->detach_from();
		TS_ASSERT( !factory->observer_attached( pose ) );
	}

	void test_residue_deletion()
	{
		core::pose::Pose pose;
		core::io::pdb::build_pose_from_pdb_as_is( pose, "protocols/denovo_design/components/helix15.pdb" );

		std::string const sd_id = "UnitTest";
		StructureData const orig_sd = *factory->create_from_pose( pose, sd_id );
		factory->attach_observer( pose, orig_sd );
		TS_ASSERT( factory->observer_attached( pose ) );

		// delete residues 1 and 2
		// observer should recognize this and fix the StructureData
		protocols::grafting::simple_movers::DeleteRegionMover delete_residues( 1, 2 );
		delete_residues.apply( pose );

		TS_ASSERT( factory->observer_attached( pose ) );
		StructureData const new_sd = factory->get_from_const_pose( pose );
		new_sd.check_pose_consistency( pose );
		TS_ASSERT_EQUALS( orig_sd.pose_length(), new_sd.pose_length() + 2 );

	}

};

