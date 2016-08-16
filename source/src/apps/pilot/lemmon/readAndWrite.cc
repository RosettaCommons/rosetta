// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/** @page readAndWrite
	This simple file describes reading and writing a PDB with Rosetta
	Run it like this:
	"readAndWritePDB.cc -in::file::s <PDB file name> -in::path::database <database root dir>"

	If you provide a PDB with a ligand that has multiple residues, this code will "glue" them together
	Of course you would need to provide information about these additional residues:
		in::file::extra_res_fa <list of extra params files, one per residue type>
*/

/// @file   apps/pilot/lemmon/readAndWrite.cc
///
/// @brief This is to illustrate reading a PDB
/// @detail Run this script with the following arguments:
/// 1) in::file::s <list of one PDB to score>
/// 2) in::path::database <list of one database root directory>
/// 3) in::file::extra_res_fa <list of extra params files, one per residue type>
/// @author Gordon Lemmon (glemmon@gmail.com)

#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>

#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>


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
		core::import_pose::pose_from_file( pose, pdb, core::import_pose::PDB_file);
	}
	{
		const std::string output("output.pdb");
		pose.dump_pdb(output);
	}
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
    }
    return 0;
}
