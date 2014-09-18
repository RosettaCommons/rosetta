// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Liz Kellogg ekellogg@u.washington.edu

// libRosetta headers

#include <core/types.hh>

#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueMatcher.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <protocols/jd2/ScoreMap.hh>

#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

// AUTO-REMOVED #include <basic/options/util.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>


#include <numeric/random/random.hh>

#include <protocols/ddg/ddGData.hh>

#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

// C++ headers
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>






namespace protocols {
namespace ddG{

ddGData::ddGData(){
	file_to_read_ = "";
	curr_mut_filename = "";
	curr_wt_filename = "";
	experimental = 0.0;
	curr_mut_read = true;
	curr_wt_read = true;
}

ddGData::ddGData(std::string filename){
	file_to_read_ = filename;
	curr_mut_filename = "";
	curr_wt_filename = "";
	experimental = 0.0;
	curr_mut_read = true;
	curr_wt_read = true;
}

//returns false if file hasn't been opened
bool ddGData::end(){
	if(!inputstream.is_open() && (file_to_read_.compare("") != 0)){
		inputstream.open(file_to_read_.c_str());
	}
	return (inputstream.eof());
}

void ddGData::get_next_filenames(){
	if(!inputstream.is_open()){
		inputstream.open(file_to_read_.c_str());
	}
	if(inputstream.is_open() && !inputstream.eof() && curr_mut_read && curr_wt_read){
		inputstream >> curr_wt_filename >> curr_mut_filename >> experimental;
	}
	curr_mut_read = false;
	curr_wt_read = false;
}

utility::vector1<pose::Pose> ddGData::read_mut_data()
{
	utility::vector1<pose::Pose> curr_mut;
	core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	if(!inputstream.eof() && !curr_mut_read){
		curr_mut_read = true;
	}
	utility::file::FileName mut(curr_mut_filename);
	if(mut.ext().compare("list") == 0){
		std::ifstream istream;
		istream.open(curr_mut_filename.c_str());
		if(istream.is_open()){
			while(!istream.eof()){
				pose::Pose newpose;
				std::string curr_filename;
				istream >> curr_filename;
				if(curr_filename.compare("") != 0){
					std::cout << "opening file: " << curr_filename << std::endl;
					core::import_pose::pose_from_pdb(newpose,(*rsd_set),curr_filename,false);
					curr_mut.push_back(newpose);
				}
			}
		}
	}else if(mut.ext().compare("out") == 0){
		core::io::silent::SilentFileData sfd;
		sfd.set_filename(curr_mut_filename);
		utility::vector1<std::string> tags = sfd.tags();
		for(unsigned int i=1;i <= tags.size(); i++){
			core::io::silent::SilentStructOP ss = sfd[tags[i]];
			pose::Pose newpose;
			ss->fill_pose(newpose,(*rsd_set));
			curr_mut.push_back(newpose);
		}
	}else{
		//extension not recognized
	}
	return curr_mut;
}

utility::vector1<pose::Pose> ddGData::read_wt_data(){
	core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	utility::vector1<pose::Pose> curr_wt;
	if(!inputstream.eof() && !curr_wt_read){
		curr_wt_read = true;
	}
	utility::file::FileName wt(curr_wt_filename);
	if(wt.ext().compare("list") == 0){
		std::ifstream istream;
		istream.open(curr_wt_filename.c_str());
		if(istream.is_open()){
			while(!istream.eof()){
				pose::Pose newpose;
				std::string curr_filename;
				istream >> curr_filename;
				if(curr_filename.compare("") != 0){
					core::import_pose::pose_from_pdb(newpose,(*rsd_set),curr_filename,false);
					curr_wt.push_back(newpose);
				}
			}
		}
	}else if(wt.ext().compare("out") == 0){
		core::io::silent::SilentFileData sfd;
		sfd.set_filename(curr_wt_filename);
		utility::vector1<std::string> tags = sfd.tags();
		for(unsigned int i=1;i <= tags.size(); i++){
			core::io::silent::SilentStructOP ss = sfd[tags[i]];
			pose::Pose newpose;
			ss->fill_pose(newpose,(*rsd_set));
			curr_wt.push_back(newpose);
		}
	}else{
		//extension not recognized
	}
	return curr_wt;
}

core::Real ddGData::read_exp_data(){
	return experimental;
}


} //namespace ddG
} //namespace protocols
