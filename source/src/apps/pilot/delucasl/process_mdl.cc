// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/apps/pilot/delucasl/process_mdl.cc
/// @author Sam DeLuca
/// @brief  AtomType and assign neighbor atoms for a list of mdl files, output fixed mdl files and pdbs

#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/sdf/mol_parser.hh>
#include <core/chemical/sdf/mol_writer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <utility/string_util.hh>
#include <list>
#include <string>

int main(int argc, char*argv[])
{
    try {
  devel::init(argc, argv);
  std::list<std::string> file_list;
  if(basic::options::option[basic::options::OptionKeys::in::file::l].user())
  {
	  std::string list_file = basic::options::option[basic::options::OptionKeys::in::file::l]()[1];


	  utility::io::izstream list;
	  list.open(list_file,std::_S_in);
	  while(!list.eof())
	  {
		  std::string line;
		  getline(list,line);
		  if(line.size() > 0)
		  {
			  file_list.push_back(line);
		  }
	  }
	  list.close();
  }else if(basic::options::option[basic::options::OptionKeys::in::file::s].user())
  {
	  std::string mol_file = basic::options::option[basic::options::OptionKeys::in::file::s]()[1];
	  file_list.push_back(mol_file);
  }


  core::chemical::AtomTypeSetCAP atom_types =
		  core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
  core::chemical::MMAtomTypeSetCAP mm_atom_types =
		  core::chemical::ChemicalManager::get_instance()->mm_atom_type_set("fa_standard");
  core::chemical::orbitals::OrbitalTypeSetCAP orbital_types =
		  core::chemical::ChemicalManager::get_instance()->orbital_type_set("fa_standard");
  core::chemical::ElementSetCAP element_set =
		  core::chemical::ChemicalManager::get_instance()->element_set("default");

  std::list<std::string>::iterator file_list_it;
  for(file_list_it = file_list.begin(); file_list_it != file_list.end();++file_list_it)
  {
	  std::string pathname = *file_list_it;
	  std::string filename = *(utility::string_split(pathname, '/').rbegin());
	  std::string base_name = *(utility::string_split(filename,'.').begin());

	  std::cout <<pathname <<std::endl;
	  core::chemical::sdf::MolFileParser parser(pathname);
	  parser.parse_mol_file(atom_types,element_set,mm_atom_types,orbital_types);
	  core::chemical::ResidueTypeOP residue = parser.GetResidueTypeOP();
	  std::cout <<"parsed molfile "<<std::endl;

	  core::chemical::sdf::MolWriter writer;
	  writer.output_residue(base_name+"_typed.mol",residue);

	  core::conformation::Residue new_residue(*residue,true);
	  core::pose::Pose pose;

	  pose.append_residue_by_jump(new_residue,1);
	  pose.dump_pdb(base_name+"_typed.pdb");
	  std::cout <<"wrote molfile" <<std::endl;

  }
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
}

