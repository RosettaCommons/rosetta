// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jobdist/PlainSilentFileJobDistributor.cxxtest.hh
/// @brief  testing PlainSilentFileJobDistributor
/// @author David Kim

// Test headers
#include <cxxtest/TestSuite.h>

#include <test/core/init_util.hh>

#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>

#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>

//Auto Headers
#include <core/io/silent/EnergyNames.fwd.hh>
#include <utility/vector1.hh>



using basic::T;
using basic::Error;
using basic::Warning;

// Package Headers

static basic::Tracer TR("protocol.jobdist.PlainSilentFileJobDistributor.cxxtest");

class PlainSilentFileJobDistributorTest : public CxxTest::TestSuite
{

public:
	PlainSilentFileJobDistributorTest() {};

	// Shared initialization goes here.
	void setUp() {
		core_init();
	}

	// Shared finalization goes here.
	void tearDown() {
	}


	void test_next_job() {
		using protocols::jobdist::BasicJob;
		using protocols::jobdist::BasicJobOP;
		using protocols::jobdist::PlainSilentFileJobDistributor;
		using namespace basic::options;
		using namespace OptionKeys;

		// There should be 7 structures in the test silent file below.
		option[out::file::silent].value("protocols/jobdist/PlainSilentFileJobDistributor_test.out");
		// The curr_nstruct (next structure) should be 8
		int actual_curr_nstruct = 8;
		// Lets loop through 9 models
		int const nstruct = 9;

		// setup JobDistributor stuff
		utility::vector1< BasicJobOP > input_jobs;
		BasicJobOP job( new BasicJob("" /*no input tag*/, "test_job", nstruct) );
		input_jobs.push_back( job );
		PlainSilentFileJobDistributor jobdist( input_jobs );
		BasicJobOP curr_job;
		int curr_nstruct;

		jobdist.startup();
		while ( jobdist.next_job(curr_job, curr_nstruct) ) {
			//std::cout << "jobdist curr_nstruct: " << curr_nstruct << " actual curr_nstruct: " << actual_curr_nstruct << std::endl;
			TS_ASSERT( curr_nstruct == actual_curr_nstruct );
			actual_curr_nstruct++;
		}
		jobdist.shutdown();

	}

};

