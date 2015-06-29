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
#include <core/conformation/Residue.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreType.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/trajectory/DbTrajectoryReader.hh>

// Utility headers
#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/inout.OptionKeys.gen.hh>

// External headers
#include <boost/foreach.hpp>

// C++ headers
#include <sstream>
#include <algorithm>

// Global Names {{{1

using namespace std;
using namespace basic::options;

using core::Size;
using core::Real;
using core::pose::Pose;
using protocols::loops::Loop;
using protocols::loops::Loops;
using protocols::loops::loop_rmsd;
using protocols::trajectory::DbTrajectoryReader;
using utility::vector1;

// Options {{{1

OPT_1GRP_KEY(Integer, query, job)

// }}}1

int main(int argc, char **argv) { // {{{1
	option.add(OptionKeys::query::job, "Job id");

	devel::init(argc, argv);

	if (not option[OptionKeys::query::job].active()) {
		string error = "Use '-query:job' to specify a job id.";
		utility_exit_with_message(error);
	}
	if (not option[OptionKeys::inout::dbms::database_name].user()) {
		option[OptionKeys::inout::dbms::database_name].value("sandbox.db");
	}

	DbTrajectoryReader reader(option[OptionKeys::query::job]());
	vector1<Size> iterations = reader.get_iterations();

	Loop loop = Loop(153, 164);
	Loops loops; loops.add_loop(loop);
	Pose reference = reader.get_pose(0);

	cout << "header ";
	cout << iterations.size() << " ";
	cout << loop.start() << " ";
	cout << loop.stop() << endl;
	
	BOOST_FOREACH(Size i, iterations) {
		cout << "  " << i << endl;
		Pose pose = reader.get_pose(i);
		Real score = pose.energies().total_energy();
		Real rmsd = loop_rmsd(pose, reference, loops);

		// Print out the current iteration.
		
		cout << "iter " << i << endl;

		// Print out score and RMSD values for the pose.

		cout << "pose " << score << " " << rmsd << endl;

		// Print out torsions for each sampled residue.

		for (Size j = loop.start(); j <= loop.stop(); j++) {
			cout << "residue " << j << " ";
			cout << pose.phi(j) << " ";
			cout << pose.psi(j) << " ";
			cout << pose.omega(j) << " ";
			for (Size k = 1; k <= pose.residue(j).nchi(); k++) {
				cout << pose.chi(k, j) << " ";
			}
			cout << endl;
		}
	}
}

