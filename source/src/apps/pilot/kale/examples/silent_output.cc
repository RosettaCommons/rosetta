// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <devel/init.hh>
#include <boost/foreach.hpp>

#include <sstream>

using namespace std;
using namespace core::io::silent;
using core::Size;
using core::Real;
using core::pose::Pose;
using core::import_pose::pose_from_pdb;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using utility::vector1;
using utility::tools::make_vector1;

int main(int argc, char** argv) {
	devel::init(argc, argv);

	// Load pose from pdb
	Pose pose; pose_from_pdb(pose, "structures/linear/5.1ubq.pdb");
	ScoreFunctionOP score_function = core::scoring::get_score_function();
	score_function->score(pose);
	cout << "Original Score: " << pose.energies().total_energy() << endl;
	pose.dump_pdb("pdb_from_pose.pdb");

	// Write pose to silent string
	string silent_pose;
	{
		stringstream string_stream;
		SilentFileData silent_file;
		SilentStructOP silent_data = new BinarySilentStruct(pose, "db");
		silent_file._write_silent_struct(*silent_data, string_stream);
		silent_pose = string_stream.str();
	}

	cout << silent_pose << endl;

	// Read pose from silent
	{	
		Pose temp_pose;
		SilentFileData silent_file;
		stringstream string_stream(silent_pose);
		vector1<string> tags = make_vector1("db");
		silent_file.read_stream(string_stream, tags, true);

		silent_file["db"]->fill_pose(temp_pose);
		silent_file["db"]->energies_into_pose(temp_pose);
		Real total_score = silent_file["db"]->get_energy("score");

		temp_pose.energies().total_energy() = total_score;
		temp_pose.energies().total_energies()[core::scoring::total_score] = total_score;

		temp_pose.dump_pdb("pdb_from_silent.pdb");
	}

}

