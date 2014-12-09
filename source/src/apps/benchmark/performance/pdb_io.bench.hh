// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/pdb_io.bench.cc
///
/// @brief  Performance benchmark for PDB input and output
/// @author Matthew O'Meara

#include <apps/benchmark/performance/performance_benchmark.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <iostream>
#include <fstream>

#include <utility/vector1.hh>


using namespace core;

class PDB_IOBenchmark : public PerformanceBenchmark
{
public:
	pose::Pose pose;

	PDB_IOBenchmark(std::string name) : PerformanceBenchmark(name) {};

	virtual void setUp() {

		std::ifstream pdb("test_in.pdb");

		pdb.seekg(0, std::ios::end);
		Size length = pdb.tellg();
		pdb.seekg(0, std::ios::beg);

		pdb_string_.resize(length);
		pdb.read(&pdb_string_[0], length);
	}

	virtual void run(core::Real scaleFactor) {
		core::Size reps( (core::Size)(100*scaleFactor) );
		if( reps == 0 ) { reps = 1; } // do at least one rep, regardless of scale factor
		for(core::Size i=0; i<reps; i++) {
			core::import_pose::pose_from_pdbstring(pose, pdb_string_);
		}
	};

	virtual void tearDown() {};

	std::string pdb_string_;
};
