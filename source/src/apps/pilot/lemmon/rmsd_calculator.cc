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

#include <devel/init.hh>
//#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/option.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/options/util.hh>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

void
chain_rmsd(core::pose::Pose & after, core::pose::Pose & before, core::Size const jump_id, bool superimpose){


	core::pose::Pose before_ligand;
	core::pose::Pose after_ligand;
	{ ///TODO make this section a function
		core::Size const before_first_residue = before.fold_tree().downstream_jump_residue(jump_id);
		core::Size const after_first_residue = after.fold_tree().downstream_jump_residue(jump_id);
		core::Size const before_chain= before.chain(before_first_residue);
		core::Size const after_chain= before.chain(after_first_residue);
		core::pose::Pose before_copy(before);
		core::pose::PoseOP after_copy= new core::pose::Pose(after);
		before_copy.remove_constraints();/// TODO fix split_by_chain to avoid this
		after_copy->remove_constraints();/// TODO fix split_by_chain to avoid this
		before_ligand= before_copy.split_by_chain(before_chain);
		after_ligand= after_copy->split_by_chain(after_chain);
	}

	if(superimpose){
		core::Real ligand_rms_with_super= core::scoring::rmsd_with_super(
				before_ligand,
				after_ligand,
				core::scoring::is_ligand_heavyatom
		);
		std::cout<< "chain rmsd super: "<< ligand_rms_with_super<< std::endl;
	}
	else{
		core::Real ligand_rms_no_super= core::scoring::rmsd_no_super(
				before_ligand,
				after_ligand,
				core::scoring::is_ligand_heavyatom
		);
		std::cout<< "chain rmsd no super: "<< ligand_rms_no_super<< std::endl;
	}
}

void
automorph_rmsd(core::pose::Pose & after, core::pose::Pose & before, core::Size res_id, bool superimpose){
	assert(before.size() == after.size());
	core::Real rms = core::scoring::automorphic_rmsd(
			before.residue(before.size()),
			after.residue(after.size()),
			superimpose
	);
	std::cout << "RMS: "<< rms<< std::endl;
}

/////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
	devel::init(argc, argv);

	std::string native_name= basic::options::option[basic::options::OptionKeys::in::file::native]();
	core::pose::PoseOP native = core::import_pose::pose_from_file(native_name, core::import_pose::PDB_file);

	utility::vector1< std::string > pdb_names= basic::options::start_files();

	std::string chain = basic::options::option[basic::options::OptionKeys::run::chain]();
	core::Size chain_id= core::pose::get_chain_id_from_chain(chain, *native);
	core::Size res_id= native->conformation().chain_begin(chain_id);
	core::Size jump_id= core::pose::get_jump_id_from_chain_id(chain_id, *native);

	utility::vector0<std::string>::iterator begin= pdb_names.begin();
	for(; begin != pdb_names.end(); ++begin){
		std::cout<< *begin << std::endl;
		core::pose::PoseOP pose = core::import_pose::pose_from_file(*begin, core::import_pose::PDB_file); // no super
		chain_rmsd(*pose, *native, jump_id, false); // no super
		chain_rmsd(*pose, *native, jump_id, true); // super
		automorph_rmsd(*pose, *native, res_id, false);// no super
		automorph_rmsd(*pose, *native, res_id, true);// super
	}

	return(0);
}
