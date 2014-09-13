// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Headers {{{1
#include <protocols/canonical_sampling/DbTrajectoryReader.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentFileData.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector.hh>
#include <utility/tools/make_vector1.hh>
#include <basic/database/sql_utils.hh>
#include <basic/Tracer.hh>

// External Headers
#include <boost/foreach.hpp>
#include <cppdb/frontend.h>

// C++ Headers
#include <sstream>
#include <algorithm>

// }}}1

namespace protocols {
namespace canonical_sampling {

// Global Names {{{1
using namespace std;

using core::Real;
using core::Size;
using core::pose::Pose;
using utility::vector1;
using utility::tools::make_vector1;
using basic::database::table_exists;
using basic::database::safely_prepare_statement;
using basic::database::safely_read_from_database;
using cppdb::statement;
using cppdb::result;
// }}}1

DbTrajectoryReader::DbTrajectoryReader() { // {{{1
	db_session_ = basic::database::get_db_session();
	job_id_ = core::Size(-1); // bazzoli: used cast to avoid warning
}

DbTrajectoryReader::DbTrajectoryReader(Size job_id) { // {{{1
	db_session_ = basic::database::get_db_session();
	job_id_ = job_id;
}

void DbTrajectoryReader::set_job_id(Size job_id) { // {{{1
	job_id_ = job_id;
}

Size DbTrajectoryReader::get_num_iterations() const { // {{{1
	return get_iterations().back();
}

vector1<Size> DbTrajectoryReader::get_iterations() const { // {{{1
	if (! table_exists(db_session_, "trajectories")) {
		utility_exit_with_message("no 'trajectories' table in the database.");
	}

	vector1<Size> iterations;
	string select_string =
		"SELECT iteration "
		"FROM trajectories "
		"WHERE job_id = ?;";

	statement select_statement =
		safely_prepare_statement(select_string, db_session_);
	select_statement.bind(1, job_id_);
	result select_result = safely_read_from_database(select_statement);

	if (! select_result.next()) {
		stringstream error_message;
		error_message << "Unable to locate job with id '" << job_id_ << "'";
		utility_exit_with_message(error_message.str());
	}

	while (! select_result.empty()) {
		Size iteration;
		select_result >> iteration;
		select_result.next();
		iterations.push_back(iteration);
	}

	// I imagine this list will already be sorted by the database, but there's no
	// reason not to be safe.  This class is meant to be used in after-the-fact
	// analysis scripts, so it's not performance-critical.

	sort(iterations.begin(), iterations.end());
	return iterations;
}

Pose DbTrajectoryReader::get_pose(Size iteration) const { // {{{1
	if (! table_exists(db_session_, "trajectories")) {
		utility_exit_with_message("no 'trajectories' table in this database.");
	}

	string select_string =
		"SELECT silent_pose "
		"FROM trajectories "
		"WHERE job_id = ? AND iteration = ?;";

	statement select_statement =
		safely_prepare_statement(select_string, db_session_);
	select_statement.bind(1, job_id_);
	select_statement.bind(2, iteration);
	result select_result = safely_read_from_database(select_statement);

	if (! select_result.next()) {
		stringstream error_message;
		error_message << "Unable to locate iteration '" << iteration << "' ";
		error_message << "in job with id '" << job_id_ << "'.";
		utility_exit_with_message(error_message.str());
	}

	string silent_pose;
	select_result >> silent_pose;
	select_result.next();

	Pose pose;
	core::io::silent::SilentFileData silent_file;
	stringstream string_stream(silent_pose);
	vector1<string> tags = make_vector1("db");
	silent_file.read_stream(string_stream, tags, true);
	silent_file["db"]->fill_pose(pose);
	silent_file["db"]->energies_into_pose(pose);

	// I can't get the silent file machinery to set the "total score" of the
	// pose, so I'm setting it manually.  I borrowed these two lines from the
	// ScoreFunction class, but this is fragile.  I also wouldn't be surprised if
	// the silent file machinery left the pose incomplete in other ways, too.

	Real score = silent_file["db"]->get_energy("score");
	pose.energies().total_energy() = score;
	pose.energies().total_energies()[core::scoring::total_score] = score;

	if (! select_result.empty()) {
		stringstream error_message;
		error_message << "More than one pose found with job id '" << job_id_;
		error_message << "' and iteration '" << iteration << "'.";
		utility_exit_with_message(error_message.str());
	}

	return pose;
}

vector1<Pose> DbTrajectoryReader::get_poses() const { // {{{1
	vector1<Pose> poses;
	vector1<Size> iterations = get_iterations();

	BOOST_FOREACH(Size i, iterations) {
		Pose pose = get_pose(i);
		poses.push_back(pose);
	}
	return poses;
}


// }}}1

} // namespace canonical_sampling
} // namespace protocols
