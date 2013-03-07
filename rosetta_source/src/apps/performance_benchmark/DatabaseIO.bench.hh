// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/pdb_io.bench.cc
///
/// @brief  Performance benchmark for database input output
/// @author Tim Jacobs

#include <apps/performance_benchmark/performance_benchmark.hh>

#include <core/pose/Pose.hh>

#include <fstream>

#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/sql_database/types.hh>

#include <protocols/features/helixAssembly/ConcurrencyTest.hh>

// Boost Headers
#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp>
#include <boost/uuid/uuid_io.hpp>

class DatabaseIOBenchmark : public PerformanceBenchmark
{
public:
	utility::sql_database::sessionOP db_session_;
	protocols::features::helixAssembly::ConcurrencyTest test_feature_;

	DatabaseIOBenchmark(std::string name) : PerformanceBenchmark(name) {};

	virtual void setUp() {
		db_session_ =
			get_db_session(
				utility::sql_database::DatabaseMode::sqlite3, "db_io_benchmark.db3");

		test_feature_.write_schema_to_db(db_session);
	}

	virtual void run(core::Real scaleFactor) {

		core::pose::Pose pose;

		boost::uuids::basic_random_generator<numeric::random::RandomGenerator>
		uuids_rng(numeric::random::RG);
		boost::uuids::uuid struct_id = uuids_rng();

		test_feature_.report_features(pose, struct_id, db_session_);
	};

	virtual void tearDown() {};

};
