// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/jd2/JD2ResourceManager.cxxtest.hh
/// @brief test suite for protocols::jd2::JD2ResourceManager and protocols::resource_manager::LazyResourceManager
/// @author Brian D. Weitzner brian.weitzner@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <protocols/jd2/InnerJob.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
// Package headers

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>
#include <utility/exit.hh>

// Numberic headers

// C++ headers
#include <string>

using namespace protocols::jd2;

class DummyObserver;
typedef utility::pointer::owning_ptr< DummyObserver > DummyObserverOP;
typedef utility::pointer::access_ptr< DummyObserver > DummyObserverAP;

class DummyObserver :
	public protocols::moves::Mover, // For OP
	public protocols::jd2::JobOutputterObserver {
public:
	DummyObserver() {
		call_counter_ = 0;
	}

	virtual ~DummyObserver() {};

	virtual
	void add_values_to_job( core::pose::Pose const &, protocols::jd2::Job & ) const {
		call_counter_++;
	}

	void reset() {
		call_counter_ = 0;
	}

	mutable core::Size call_counter_;
	
	// Dummy for OP
	void apply( Pose & ) {}
	virtual std::string get_name() const { return "DummyObserver"; }
	void add_ref() const {}
	void remove_ref() const {}
};

class JobOutputterTests : public CxxTest::TestSuite {

private:

public:

	void setUp() {
		protocols_init();
	}

	// @brief test default options and default locator
	void test_Observer_add_remove() {
		DummyObserverOP observer1_op = new DummyObserver;
		DummyObserverOP observer2_op = new DummyObserver;
		DummyObserverOP observer3_op = new DummyObserver;
		DummyObserver & observer1 = *observer1_op;
		DummyObserver & observer2 = *observer2_op;
		DummyObserver & observer3 = *observer3_op;
		SilentFileJobOutputter job_outputter;
		core::pose::Pose pose;
		JobOP job = new Job( InnerJobOP( new InnerJob( "job", 1 ) ), 1 );
		job_outputter.call_output_observers( pose, job  );
		TS_ASSERT( observer1.call_counter_ == 0 );
		job->add_output_observer( DummyObserverAP(observer1_op) );
		job->add_output_observer( DummyObserverAP(observer2_op) );
		job->add_output_observer( DummyObserverAP(observer3_op) );
		job_outputter.call_output_observers( pose, job  );
		TS_ASSERT( observer1.call_counter_ == 1 );
		TS_ASSERT( observer2.call_counter_ == 1 );
		TS_ASSERT( observer3.call_counter_ == 1 );
		//adding them twice should still only lead to a single call
		job->add_output_observer( DummyObserverAP(observer1_op) );
		job->add_output_observer( DummyObserverAP(observer2_op) );
		job->add_output_observer( DummyObserverAP(observer3_op) );
		job_outputter.call_output_observers( pose, job  );
		TS_ASSERT( observer1.call_counter_ == 2 );
		TS_ASSERT( observer2.call_counter_ == 2 );
		TS_ASSERT( observer3.call_counter_ == 2 );
		observer1.reset();
		observer2.reset();
		observer3.reset();
		TS_ASSERT( observer1.call_counter_ == 0 );
		TS_ASSERT( observer2.call_counter_ == 0 );
		TS_ASSERT( observer3.call_counter_ == 0 );
		//removing them should stop them from being called
		job->remove_output_observer( DummyObserverAP(observer1_op) );
		job->remove_output_observer( DummyObserverAP(observer2_op) );
		job->remove_output_observer( DummyObserverAP(observer3_op) );
		job_outputter.call_output_observers( pose, job  );
		TS_ASSERT( observer1.call_counter_ == 0 );
		TS_ASSERT( observer2.call_counter_ == 0 );
		TS_ASSERT( observer3.call_counter_ == 0 );
	}
};
