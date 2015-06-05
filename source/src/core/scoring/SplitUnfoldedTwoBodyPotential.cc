// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/SplitUnfoldedTwoBodyPotential.cc
/// @brief  Reads in and stores the two body energies for each residue type for the split unfolded energy.
/// @author Riley Simmons-Edler (rse231@nyu.edu)


#include <core/scoring/SplitUnfoldedTwoBodyPotential.hh>

#include <core/pose/Pose.hh>

#include <core/conformation/Residue.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreTypeManager.hh>
#include <core/scoring/EnergyMap.hh>

#include <core/chemical/ResidueType.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>

#include <ObjexxFCL/string.functions.hh>

#include <iostream>
#include <string>
#include <map>


namespace core {
namespace scoring {


  SplitUnfoldedTwoBodyPotential::SplitUnfoldedTwoBodyPotential(std::string filename)
  {
    read_database_file(filename);
  }

  SplitUnfoldedTwoBodyPotential::SplitUnfoldedTwoBodyPotential(std::string filename,std::string atom_type_label_name)
  {
		atom_type_label_set_used_=atom_type_label_name;
    read_database_file(filename);
  }

  void SplitUnfoldedTwoBodyPotential::get_restype_emap(const chemical::ResidueType & restype,EnergyMap & emap) const
  {
    std::string resname=restype.name();
    std::map<std::string,EnergyMap>::const_iterator iter=residue_two_body_energies_.find(resname);
    if(iter != residue_two_body_energies_.end())
      {
				EnergyMap retmap=iter->second;
				emap=retmap;
      }
    else
      {
				//	residue_two_body_energies_[resname]=calculate_residue_emap(restype);
				//GRR CONST METHODS NOT LETTING ME SAVE THIS!!!!!
				//				residue_two_body_energies_.insert(std::pair<std::string,core::scoring::EnergyMap>(resname,calculate_residue_emap(restype)));
				//				residue_two_body_energies_[resname]=calculate_residue_emap(restype);
				//emap=residue_two_body_energies_.find(resname)->second;
				emap=calculate_residue_emap(restype);
      }
  }

  EnergyMap SplitUnfoldedTwoBodyPotential::get_weights() const
  {
    return residue_score_term_weights_;
  }

  EnergyMap SplitUnfoldedTwoBodyPotential::calculate_residue_emap(const chemical::ResidueType & restype) const
  {
    Size natoms=restype.natoms();
    EnergyMap retmap;
    for(Size i=1;i<natoms;i++)
      {
				//Determine what atom type set we are using and get the appropriate atom name for this atom.
				std::string curratomname;
				if(atom_type_label_set_used_=="rosetta")
					{
						curratomname=restype.atom_type(i).name();
					}
				else if(atom_type_label_set_used_=="mm")
					{
						curratomname=restype.atom(i).mm_name();
					}
				else if(atom_type_label_set_used_=="elemental")
					{
						curratomname=restype.atom_type(i).element();
					}
				else if(atom_type_label_set_used_=="pdb")
					{
						curratomname=restype.atom_name(i);
						std::string::iterator end_pos = std::remove(curratomname.begin(), curratomname.end(), ' ');
						curratomname.erase(end_pos, curratomname.end());
					}
				else if(atom_type_label_set_used_=="unique" && restype.atom_type(i).name()!="VIRT") //skip virtual atoms, since we won't find distributions for those anyway and they keep throwing error messages. :(
					{
						//using unpatched names since patched residues don't work with unique atom types.
						curratomname=restype.name3();
						curratomname.append("/");
						std::ostringstream tostring;
						tostring<<i;
						curratomname+=tostring.str();
					}
				else
					{
						curratomname=restype.atom(i).mm_name();
					}
				//				std::cout<<curratomname<<std::endl;

	//now need to add the energies for this atom to the running total
	//not sure this is the best way to do this...
	//	EnergyMap curratommap=atom_two_body_energies_[curratomname];
	//hacky workaround to avoid using the [] operator, which const functions disallow for some reason.
				std::map<std::string,core::scoring::EnergyMap>::const_iterator curratomiter=atom_two_body_energies_.find(curratomname);
				if(curratomiter==atom_two_body_energies_.end())
					{
						std::cout<<"ENERGY FOR ATOM TYPE "<<curratomname<<" NOT FOUND!"<<std::endl;
					}
				EnergyMap curratommap=curratomiter->second;
				/*for(EnergyMap::iterator iter=curratommap.begin();iter!=curratommap.end();iter++)
					{
						std::cout<<iter<<std::endl;
						//	    retmap.set(iter->first,retmap.get(iter->first)+iter->second);
						}*/
				//set both ref and regular values(regular ones are used to generate combined two body term in SplitUnfoldedTwoBodyEnergy).
				retmap.set(fa_atr_ref,retmap.get(fa_atr_ref)+curratommap.get(fa_atr));
				retmap.set(fa_rep_ref,retmap.get(fa_rep_ref)+curratommap.get(fa_rep));
				retmap.set(fa_sol_ref,retmap.get(fa_sol_ref)+curratommap.get(fa_sol));
				retmap.set(fa_elec_ref,retmap.get(fa_elec_ref)+curratommap.get(fa_elec));
				retmap.set(hbond_ref,retmap.get(hbond_ref)+curratommap.get(hbond_sc));
				retmap.set(dslf_fa13_ref,retmap.get(dslf_fa13_ref)+curratommap.get(dslf_fa13));

				retmap.set(fa_atr,retmap.get(fa_atr)+curratommap.get(fa_atr));
				retmap.set(fa_rep,retmap.get(fa_rep)+curratommap.get(fa_rep));
				retmap.set(fa_sol,retmap.get(fa_sol)+curratommap.get(fa_sol));
				retmap.set(fa_elec,retmap.get(fa_elec)+curratommap.get(fa_elec));
				retmap.set(hbond_sc,retmap.get(hbond_sc)+curratommap.get(hbond_sc));
				retmap.set(dslf_fa13,retmap.get(dslf_fa13)+curratommap.get(dslf_fa13));

      }
    return retmap;
  }

  //this function heavily adapted from UnfoldedStatePotential.cc
  void SplitUnfoldedTwoBodyPotential::read_database_file(std::string filename)
  {
		std::cout<<filename<<std::endl;
    //check for file existence
    if(!utility::file::file_exists(filename)) 
      {
	utility_exit_with_message("Cannot find file '"+filename+"'");
      }

    //check for data "goodness" :)
    utility::io::izstream data(filename);
    if (!data.good()) // data is no good! :(
      {
	utility_exit_with_message("Cannot open file '"+filename+"'");
      }

    // read in all lines in file
    utility::vector1<std::string> lines;
    std::string line;
    while(getline(data,line)) 
      {
	std::istringstream l(line);
	if(line.size() < 1 || line[0] == '#') continue; // skip comment lines
	lines.push_back(line);
      }
    data.close();

    // parse the first line that contains the score types
    std::string const & first_line(lines[1]);
    std::istringstream h(first_line);
    std::string tag;
    utility::vector1<ScoreType> atom_score_types;

    h >> tag;
    if(tag!="ATOM") 
      {
	utility_exit_with_message("Error parsing first line of '"+filename+"'");
      }
    while (!h.fail()) 
      {
	h >> tag;
	if (h.fail()) break;
	if (ScoreTypeManager::is_score_type(tag))  
	  {
	    atom_score_types.push_back(ScoreTypeManager::score_type_from_name(tag));
	  } 
	else 
	  {
	    utility_exit_with_message("Error, score type '"+tag+"' specified in '"+filename+"' doesn't exist.");
	  }
      }

    // parse the second line that contains the weights
    std::string const & second_line(lines[2]);
    std::istringstream k(second_line);
    Size const ntypes(atom_score_types.size());
    Real weight;
    
    k >> tag;
    if(tag!="WEIGHT") 
      {
	utility_exit_with_message("Error parsing second line of '"+filename+"'");
      }
    for(Size i=1;i<=ntypes;++i) 
      {
	k >> weight;
	if(k.fail()) 
	  {
	    utility_exit_with_message("Error, number of energies doesn't match number of score types in '"+filename+"'");
	  }
	residue_score_term_weights_[atom_score_types[i]] = weight; // add energy to appropriate score type in temp emap
      }
    
    
    // parse the rest of the file
    Size const nlines(lines.size());

    for(Size i=3;i<=nlines;++i) 
      {
	std::string const & temp_line(lines[i]);
	std::istringstream l(temp_line);
	std::string tlc;
	Real val;
	EnergyMap emap;
	l >> tlc;
	ObjexxFCL::lpad(tlc, 1); // make sure we have *something*(add a blank space if there is no atom type name, should never happen under normal conditions though)
	for(Size i=1;i<=ntypes;++i) 
	  {
	    l >> val;
	    if(l.fail()) 
	      {
		utility_exit_with_message("Error, number of energies doesn't match number of score types in '"+filename+"'");
	      }
	    emap[atom_score_types[i]] = val; // add energy to appropriate score type in temp emap
	  }
      
	// add tlc/emap pair to map
	atom_two_body_energies_[tlc] = emap;
      }	
		
  }


} //scoring
} //core
