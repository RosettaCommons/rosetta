// -*-
// mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t
// -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/** @page simpleScore
	This simple script scores a PDB
	It reads in a PDB, scores it and prints out the scored PDB
	Run it by typing:
		scorePDB.cc -in::file::s <pdb file> -in::path::database <DB root dir>
*/

/// @file   apps/pilot/lemmon/simpleScore.cc
///
/// @brief This is to illustrate scoring a PDB
/// @detail Run this script with the following arguments:
/// 1) in::file::s <list of one PDB to score>
/// 2) in::path::database <list of one database root directory>
/// 3) in::file::extra_res_fa <list of extra params files, one per residue type>
/// @author Gordon Lemmon (glemmon@gmail.com)

#include <utility/vector0.hh>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <basic/database/open.hh>
#include <utility/excn/Exceptions.hh>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>



//////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
    try {
  devel::init(argc, argv);

  utility::vector0<std::string> pdbs;
  {// process the options
    using namespace basic::options::OptionKeys;
    using basic::options::option;
    pdbs= option[in::file::s]();
  }
	core::scoring::ScoreFunction scoreFunction;
	{
		const std::string weightsPath("/scoring/weights/standard.wts");
		scoreFunction.initialize_from_file(basic::database::full_name(weightsPath));
	}
  core::pose::Pose pose; // starts NULL, coords *never* modified!
	{
		std::string pdb=pdbs[0];
		core::import_pose::pose_from_pdb( pose, pdb);
	}
	{
		const std::string output("output.pdb");
		pose.dump_scored_pdb(output, scoreFunction);
	}
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
    }
    return 0;
}
