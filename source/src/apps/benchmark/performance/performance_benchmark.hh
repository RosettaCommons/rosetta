// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/apps/benchmark/performance/performance_benchmark.hh
///
/// @brief  Base class for mini rosetta benchmark system
/// @author Sergey Lyskov

#ifndef INCLUDED_apps_benchmark_performance_benchmark_hh
#define INCLUDED_apps_benchmark_performance_benchmark_hh

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <vector>
#include <string>

static THREAD_LOCAL basic::Tracer TR( "benchmark" );


class PerformanceBenchmark
{
public:
	PerformanceBenchmark(std::string name) : result_(0), time_(0.0), name_(name) {
		allBenchmarks().push_back(this);
		for ( unsigned int i=0; i<name_.size(); i++ ) {
			if ( name_[i]=='.' || name_[i]==' ' ) name_[i]='_';
		}
	}
	virtual ~PerformanceBenchmark() {}

	virtual void setUp() {}
	virtual void run(core::Real /*scaleFactor*/) {}
	virtual void tearDown() {}

	/// Execute benchmark cycle i.e.: setUp() - run() - tearDown()
	/// return number of seconds that was used to performe 'run' step.
	double execute(core::Real scaleFactor);
	std::string name() { return name_; }

public:
	static void executeOneBenchmark(
		std::string const & name,
		core::Real scaleFactor=1);
	static void executeAllBenchmarks(core::Real scaleFactor=1);
	//static void perform_until_set_found ( PerformanceBenchmark * B, core::Real scaleFactor );
	static std::string getReport();
	static std::string getOneReport(std::string const & name);

private:
	int result_;
	double time_;
	std::string name_; ///< name of the benchmark, must corelate to namespace ie: core.pose

	/// function for keepig record of all created benchmark classes.
	static std::vector< PerformanceBenchmark * > & allBenchmarks();
};


#endif
