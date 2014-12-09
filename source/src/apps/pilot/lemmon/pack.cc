// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/** @page pack
	This code reads in, packs, and prints a PDB
	Run it by typing:
	packPDB -in::file::s <pdb file> -in::path::database <DB root dir>
*/

/// @file   apps/pilot/lemmon/pack.cc
///
/// @brief This is to illustrate packing residues in a PDB
/// @detail Run this script with the following arguments:
/// 1) in::file::s <list of one PDB with at least 2 residues to pack>
/// 2) in::path::database <list of one database root directory>
/// 3) in::file::extra_res_fa <list of extra params files, one per residue type>
/// @author Gordon Lemmon (glemmon@gmail.com)

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>

#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>


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
  core::pose::Pose pose; // starts NULL, coords *never* modified!
	{
		std::string pdb=pdbs[0];
		core::import_pose::pose_from_pdb(pose, pdb);
	}
	protocols::simple_moves::PackRotamersMover mover;
	mover.apply(pose);
	{
		const std::string output("output.pdb");
		core::scoring::ScoreFunctionCOP score_fxn= mover.score_function();
		pose.dump_scored_pdb(output, *score_fxn);
	}

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
    }
    return 0;
}
