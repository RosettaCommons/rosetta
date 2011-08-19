// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/benchmark.hh
///
/// @brief  Base class for mini rosetta benchmark system
/// @author Sergey Lyskov

#ifndef INCLUDED_apps_benchmark_benchmark_hh
#define INCLUDED_apps_benchmark_benchmark_hh

#include <vector>
#include <string>

class Benchmark
{
public:
	Benchmark(std::string name) : result_(0.), name_(name) {
		allBenchmarks().push_back(this);
		for(unsigned int i=0; i<name_.size(); i++) {
			if( name_[i]=='.' || name_[i]==' ' ) name_[i]='_';
		}
	};
	virtual ~Benchmark() {};

	virtual void setUp() {};
	virtual void run(int /*scaleFactor*/) {};
	virtual void tearDown() {};

	/// Execute benchmark cycle i.e.: setUp() - run() - tearDown()
	/// return number of seconds that was used to performe 'run' step.
	double execute(int scaleFactor);
	std::string name() { return name_; };

public:
	static void executeAllBenchmarks(int scaleFactor=1);
	static std::string getReport();

private:
	double result_;
	std::string name_; ///< name of the benchmark, must corelate to namespace ie: core.pose

	/// function for keepig record of all created benchmark classes.
	static std::vector<Benchmark *> &allBenchmarks();
};


#endif // INCLUDED_demo_rosetta_benchmark_HH
