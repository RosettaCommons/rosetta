// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Headers {{{1
#include <devel/init.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentFileData.hh>

// Protocol headers
#include <protocols/trajectory/DbTrajectoryReader.hh>
#include <protocols/canonical_sampling/PDBTrajectoryRecorder.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <basic/database/sql_utils.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// External headers
#include <boost/foreach.hpp>
#include <boost/noncopyable.hpp>
#include <cppdb/frontend.h>

// C++ headers
#include <sstream>
#include <algorithm>

// Global Names {{{1
#define foreach BOOST_FOREACH

using namespace std;
using namespace basic::options;

using core::Size;
using core::Real;
using core::pose::Pose;
using protocols::trajectory::DbTrajectoryReader;
using protocols::canonical_sampling::PDBTrajectoryRecorder;
using utility::vector1;
using utility::tools::make_vector1;
using utility::sql_database::sessionOP;
using basic::database::table_exists;
using basic::database::safely_prepare_statement;
using basic::database::safely_read_from_database;
using cppdb::statement;
using cppdb::result;

// Options {{{1

OPT_1GRP_KEY(Integer, movie, job)
OPT_1GRP_KEY(File, movie, out)
OPT_1GRP_KEY(Integer, movie, frames)

// }}}1

int main(int argc, char **argv) {
	option.add(OptionKeys::movie::job, "Job id");
	option.add(OptionKeys::movie::out, "PDB movie path").def("movie.pdb");
	option.add(OptionKeys::movie::frames, "Number of frames").def(0);

	devel::init(argc, argv);

	if (not option[OptionKeys::movie::job].active()) {
		utility_exit_with_message("Use -movie:job to specify a job id.");
	}
	if (not option[OptionKeys::inout::dbms::database_name].user()) {
		option[OptionKeys::inout::dbms::database_name].value("sandbox.db");
	}

	DbTrajectoryReader reader;
	PDBTrajectoryRecorder writer;
	
	reader.set_job_id(option[OptionKeys::movie::job]());
	writer.file_name(option[OptionKeys::movie::out]().name());
	
	vector1<Size> iterations = reader.get_iterations();
	Size const frames = option[OptionKeys::movie::frames]();
	Size const frequency = (frames == 0) ?
		1 : max<Size>(1, iterations.size() / frames);
	Size const expected = iterations.size() / frequency;
	Size visited = 0, recorded = 0;
	
	foreach (Size i, iterations) {
	cout << 3 << endl;
		if (visited++ % frequency != 0) continue;
		Pose pose = reader.get_pose(i); writer.apply(pose);
		cout << "\r[" << ++recorded << "/" << expected << "]" << flush;
	}
	cout << endl;
}
