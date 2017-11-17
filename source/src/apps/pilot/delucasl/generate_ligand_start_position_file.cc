// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/delucasl/generate_ligand_start_position_file.cc
/// @author Sam DeLuca

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Residue.hh>

#include <devel/init.hh>

#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <utility/json_spirit/json_spirit_writer.h>
#include <utility/json_spirit/json_spirit_writer_options.h>
#include <utility/tools/make_vector.hh>
#include <utility/io/ozstream.hh>
#include <string>
#include <iostream>

OPT_KEY(String,ligand_chain)
OPT_KEY(String,start_position_file)

int main(int argc, char*argv[])
{
	try {
		NEW_OPT(ligand_chain,"the name of the ligand chain","");
		NEW_OPT(start_position_file,"the name of the file to output to","");

		devel::init(argc,argv);


		char ligand_id = basic::options::option[basic::options::OptionKeys::ligand_chain]()[0];
		std::string outfile_path =  basic::options::option[basic::options::OptionKeys::start_position_file]();

		utility::vector1<std::string> file_names = basic::options::start_files();

		std::vector<utility::json_spirit::Value> json_data;

		for ( utility::vector1<std::string>::iterator file_it = file_names.begin(); file_it != file_names.end(); ++file_it ) {
			core::pose::PoseOP current_pose(core::import_pose::pose_from_file(*file_it,false, core::import_pose::PDB_file));

			std::string protein_hash = core::pose::get_sha1_hash_excluding_chain(ligand_id,*current_pose);
			if ( !core::pose::has_chain(ligand_id,*current_pose) ) {
				utility_exit_with_message(*file_it + " does not contain chain "+ ligand_id);
			}
			core::Size ligand_chain_id = core::pose::get_chain_id_from_chain(ligand_id,*current_pose);
			core::conformation::ResidueCOPs ligand_residues(core::pose::get_chain_residues(*current_pose,ligand_chain_id));
			std::cout << *file_it <<std::endl;
			if ( ligand_residues.size() != 1 ) {
				std::cout <<ligand_residues.size() <<std::endl;
				utility_exit_with_message("you can only have one residue in a ligand right now sorry");
			}

			core::Vector nbr_atom_coords = ligand_residues[1]->nbr_atom_xyz();

			utility::json_spirit::Pair in_tag("input_tag",*file_it);
			utility::json_spirit::Pair x_coord("x",nbr_atom_coords.x());
			utility::json_spirit::Pair y_coord("y",nbr_atom_coords.y());
			utility::json_spirit::Pair z_coord("z",nbr_atom_coords.z());
			utility::json_spirit::Pair hash_record("hash",protein_hash);

			utility::json_spirit::Value data_point(utility::tools::make_vector(in_tag,x_coord,y_coord,z_coord,hash_record));
			json_data.push_back(data_point);

		}

		utility::io::ozstream outfile;
		outfile.open(outfile_path.c_str());
		outfile <<utility::json_spirit::write(json_data,utility::json_spirit::pretty_print);
		outfile.close();
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
