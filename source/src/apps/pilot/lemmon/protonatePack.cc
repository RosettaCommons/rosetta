// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @page protonatePack

/// @file   apps/pilot/lemmon/protonatePack.cc
///
/// @brief This is to illustrate reading in a ligand with multiple residues and
///  then reassembling the ligand (using Rosetta).
/// @detail Run this script with the following arguments:
/// 1) in::file::s <list of PDBs with ligand residues to glue>
/// 2) in::path::database <root director of minirosetta database>
/// 3) in::file::extra_res_fa <list of extra params files, one per residue type>
/// @author Gordon Lemmon (glemmon@gmail.com)


#include <devel/init.hh>
#include <protocols/ligand_docking/LigandBaseProtocol.hh>
#include <Objecxx/FArray1D.hh>

#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>


// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>



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
  std::string pdb=pdbs[0];
  core::import_pose::pose_from_pdb( pose, pdb);

	protocols::ligand_docking::LigandBaseProtocol protocol();

	// the next lines uses Eon's code to deal with protonation states
	//ObjexxFCL::FArray1D_bool allow_repack(num_residues, false);
	//core::pack::task::PackerTaskOP task = protocol.make_packer_task(pose, allow_repack, ligand_protonation);
//OR
	// Let the ligand_docking protocol do the work
	bool ligand_protonation=true;
	bool include_all_rsds=true;
	core::Real sc_padding= 0;
	int jump_id=1;

	core::pack::task::PackerTaskOP task=
		protocol.make_packer_task(pose, jump_id, sc_padding, allow_repack);
//OR
	//core::pack::task::PackerTaskOP task=
	//protocol.make_packer_task_ligand_only(pose, jump_id, ligand_protonation);

  const std::string output("output.pdb");
  pose.dump_pdb(output);
    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
    }
    return 0;
}
