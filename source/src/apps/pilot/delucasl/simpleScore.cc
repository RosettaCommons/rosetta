// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/delucasl/simpleScore.cc
///
/// @brief This scores a pose
/// @detail Run this script with the following arguments:
/// 1) in::file::s <list of one PDB to score>
/// 2) in::path::database <list of one database root directory>
/// 3) score::weights <name of weight file (without extension)>
/// @author Sam DeLuca (samuel.l.deluca@vanderbilt.edu)

#include <utility/vector0.hh>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/database/open.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


int main(int argc, char* argv[])
{
    try {
	devel::init(argc, argv);

	utility::vector0<std::string> pdbs;
	{
		using namespace basic::options::OptionKeys;
		using basic::options::option;
		pdbs = option[in::file::l]();
	}

	core::scoring::ScoreFunction scoreFunction;
	{
		using namespace basic::options::OptionKeys;
		using basic::options::option;
		const std::string weightstag(option[score::weights]());
		scoreFunction.initialize_from_file(basic::database::full_name("scoring/weights/"+weightstag+".wts"));
	}

	core::pose::Pose pose;
	for(int structIndex = 0; structIndex < pdbs.size(); ++structIndex)
	{
		core::import_pose::pose_from_pdb(pose,pdbs[structIndex]);
		pose.dump_scored_pdb(pdbs[structIndex],scoreFunction);
	}
	/*
	 {
	 std::string pdb = pdbs[0];
	 core::import_pose::pose_from_pdb(pose,pdb);
	 }

	 {
	 const std::string output("output.pdb");
	 pose.dump_scored_pdb(output,scoreFunction);
	 }
	 */
    } catch ( utility::excn::EXCN_Base const & e ) {
                             std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                }
	return(0);
}

