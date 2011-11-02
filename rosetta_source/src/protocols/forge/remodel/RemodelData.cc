// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//  CVS information:
//  $Revision: 13011 $
//  $Date: 2007-02-21 17:17:13 -0800 (Wed, 21 Feb 2007) $
//  $Author: possu $

// Rosetta Headers
#include <basic/Tracer.hh>
#include <core/io/pdb/pose_io.hh>

#include <protocols/forge/remodel/RemodelData.hh>
//#include <devel/remodel/helpMenu.hh>

//for DSSP
#include <protocols/jumping/Dssp.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/remodel.OptionKeys.gen.hh>

//fragset
// AUTO-REMOVED #include <core/fragment/OrderedFragSet.hh>

// AUTO-REMOVED #include <protocols/viewer/viewers.hh>

// Auto-header: duplicate removed #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/chemical/ResidueType.hh>
 // for switch typeset


// for resfile command map
#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/ResfileReader.fwd.hh>

// AUTO-REMOVED #include <protocols/forge/components/VarLengthBuild.hh>

// AUTO-REMOVED #include <protocols/loops/LoopMover_QuickCCD_Moves.hh>
// AUTO-REMOVED #include <protocols/forge/build/BuildInstruction.hh> // REQUIRED FOR WINDOWS

/*
//yab headers
#include "AtomPoint.hh"
#include "BoundingBox.hh"
#include "epigraft_functions.hh"
#include "Octree.hh"
#include "rootstock_types.hh"
#include "ccd_functions.hh"
*/
// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>

// C++ Headers
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <map>
#include <set>


// Utility Headers
// AUTO-REMOVED #include <utility/basic_sys_util.hh>
#include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/io/ocstream.hh>
// AUTO-REMOVED #include <utility/file/file_sys_util.hh>
#include <utility/vector1.hh>

#include <core/import_pose/import_pose.hh>



//////////////// REMODEL

namespace protocols{
namespace forge{
namespace remodel{

static basic::Tracer TR_REMODEL("REMODELd");

void
protocols::forge::remodel::RemodelData::splitString(std::string str, std::string delim, std::vector<std::string> & results) {
	int cutAt;
	while((cutAt = (int)str.find_first_of(delim)) != int(str.npos) ) {
		if(cutAt > 0) {
			results.push_back(str.substr(0,cutAt));
		}
		str = str.substr(cutAt+1);
	}
	if(str.length() > 0) {
		results.push_back(str);
	}
}

protocols::forge::remodel::RemodelData::RemodelData(){
	has_design_info_ = false;
	design_neighbor = false;
	auto_design = false;
	natro_movemap_.set_chi(true);
}

void
protocols::forge::remodel::RemodelData::getLoopsToBuildFromFile()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

// read blueprint file and load everything into the maps
	std::string filename(option[basic::options::OptionKeys::remodel::blueprint]());

	if (filename == ""){
		TR_REMODEL << "can't find blueprint file for remodel!" << std::endl;
	}
	utility::io::izstream data(filename.c_str());
	if (!data) {
		TR_REMODEL << "Can't open blueprint file " << filename << std::endl;
		utility::exit(EXIT_FAILURE, __FILE__, __LINE__);
	}
	std::string line;

	//extension management
	std::string ext_ss_buffer;

	std::ostringstream oss; // for resfile parsing
	//std::ostringstream oss_switch; // for detecting design info
	oss << "NATRO" << std::endl; // preseve non designed to starting rotamer
	oss << "start" << std::endl; // mark start for resfile parser

	int index = 1;
	bool mark_start = false;
	while (getline( data, line)){
		std::istringstream line_stream(line);
		std::vector<std::string> split_info;
	std::ostringstream oss_switch; // for detecting design info
		this->splitString(line_stream.str(), " ", split_info);

		// skip comment lines
		if (split_info[0].at(0) == '#'){
			continue;
		}

		protocols::forge::remodel::LineObject line;
		line.isDesignable = false; // initialize design default to false
	//	TR << "index: " << index << std::endl;
		line.index = index;

		// debug
		//		line_stream >> line.original_index >> line.resname >> line.sstype >> skip;
		//		this->blueprint.push_back(line);
		//
		index++;
		// could have initialized blueprint after the split, oh well...
		std::istringstream(split_info[0]) >> line.original_index;
		if (line.original_index != 0 && mark_start == false){
			pdb_start = line.original_index;
			mark_start = true;
		}
		if (line.original_index != 0){
			pdb_stop = line.original_index;
		}

		line.resname = split_info[1];
		line.sstype = split_info[2];

		// error checking, disallow '#' in columns 1 and 2
		if (split_info[1].at(0) == '#' || split_info[2].at(0) == '#') {
			std::ostringstream err_message;
			err_message << "ERROR: comment marker '#' cannot be in residue or ss column at line:\n";
			err_message << line_stream.str();
			utility::exit(__FILE__, __LINE__, err_message.str());
		}

		if (split_info.size() > 3){ // has design info
			// skip comments at end of line
			if (split_info[3].at(0) == '#'){
				continue;
			}

			for (std::vector<std::string>::iterator it=split_info.begin(), end=split_info.end(); it != end ; it++){
				if ((*it).substr(0,3) == "CST"){
					TR_REMODEL << "constraint found " << *it <<  std::endl;
					line.has_constraints=true;
					line.constraint_definition.push_back(*it);
				}
				if ((*it).substr(0,3) == "DM_"){
					disulfMobileRange.push_back(line.index);
					if (disulfMobileRange.size() > 2){
						std::ostringstream err_message;
						err_message << "ERROR: Disulfide mobile range assigment contains " <<  disulfMobileRange.size() << " elements." << std::endl;
					  utility::exit(__FILE__, __LINE__, err_message.str());
					}
				}
				if ((*it).substr(0,3) == "DS_"){
					disulfLandingRange.push_back(line.index);
					if (disulfLandingRange.size() > 2){
						std::ostringstream err_message;
						err_message << "ERROR: Disulfide landing range assigment contains " <<  disulfLandingRange.size() << " elements." << std::endl;
					  utility::exit(__FILE__, __LINE__, err_message.str());
					}
				}
			}

			std::map< std::string, core::pack::task::ResfileCommandOP > resfile_command_map = core::pack::task::create_command_map();

			oss << line.index << " _ " ;
			bool pickaa = false;
			for (int i = 3; i< (int)split_info.size();  i++){
				if (split_info[i].substr(0,3) != "CST"){
					oss << split_info[i] << " " ;
					oss_switch << split_info[i];
				}
				if (split_info[i].substr(0,5) == "PIKAA"){
					//toggle on manual residue selection switch
					pickaa = true;
					continue;
				}
				if (pickaa){ // the column following PIKAA
					for (int j = 0; j < (int)split_info[i].size();++j){ // only find string element size

						core::chemical::AA aa;
						//char one_letter_name(core::chemical::aa_from_oneletter_code(aa));
						char one_letter_name = split_info[i].substr(j,1).c_str()[0];
						aa = core::chemical::aa_from_oneletter_code(one_letter_name);
											TR_REMODEL << "  design position to " << one_letter_name << " " << aa << std::endl;
						line.aminoAcidList.push_back(aa);
					}
					pickaa=false; // turns it right off so doesn't get other columns
				}

				if ( split_info[i] == "NATRO"){
					TR_REMODEL << "NATRO movemap setup: turning off chi move for refinement stage: " << line.index << std::endl;
					natro_movemap_.set_chi(line.index, false);
				}
			}
			oss << std::endl;

			//find out that there's info other than CST and turn on manual modes
			if (oss_switch.str() != ""){
				TR_REMODEL << "oss_switch: " << oss_switch.str() << std::endl;
					has_design_info_ = true;
			}


			std::cout << "DEBUG parsed STRING " << oss.str() << std::endl;
			this->parsed_string_for_resfile = oss.str();
			if (option[ OptionKeys::remodel::repeat_structuer].user()){
				for (int i = 1; i<= option[ OptionKeys::remodel::repeat_structuer ]; i++){
					this->parsed_string_for_resfile.append(this->parsed_string_for_resfile);
				}
			}

			//TR_REMODEL << "manual design overwrite position: " << line.index << std::endl;
			//this->design_mode = 3; //default manual mode
			/*if (option[Remodel::Design::design_neighbors]()){
				// fully manual design mode automatically switched on when you assign residues by hand
				this->design_mode = 4;
			}
			if (option[Remodel::Design::neighbor_repack]()){
				// bc repack neigbors
				this->design_mode = 5;
			}
			*/
			line.isDesignable = true;
			line.design_type = split_info[3];
			if (!resfile_command_map[line.design_type]){
				TR_REMODEL << "WARNING: unknown packer token: " << line.design_type << std::endl;
			}

			//debug
			//TR << resfile_command_map[split_info[3]] << " resfile command map to " << split_info[3] << std::endl;

/* BUGGY
			if (split_info.size() > 4 && resfile_command_map[split_info[3]]) { // has manual amino acid assignment
				for ( int i = 4 ; i < (int)split_info.size(); i++) {
					// skip comments at end of line
					if (split_info[i].at(0) == '#'){
						break;
					}
					core::chemical::AA aa;
					char one_letter_name(core::chemical::aa_from_oneletter_code(aa));
					one_letter_name = split_info[i].c_str()[0];
					aa = core::chemical::aa_from_oneletter_code(one_letter_name);
										TR_REMODEL << "  design position to " << one_letter_name << " " << aa << std::endl;
					line.aminoAcidList.push_back(aa);
				}
			}
			*/


		}
		this->blueprint.push_back(line);
	}

	if ( has_design_info_ ){
		design_mode = 3;
	}

	//process blueprint to initialize all the needed strings/vectors
	std::vector<protocols::forge::remodel::LineObject>::iterator iter;
	for ( iter = this->blueprint.begin(); iter != this->blueprint.end(); iter++) {
		this->sequence.append(iter->resname);
		this->ss.append(iter->sstype);
	}

	TR_REMODEL << "sequence: " << std::endl << this->sequence << std::endl;
	TR_REMODEL << "sstype  : " << std::endl << this->ss << std::endl;

}

void
protocols::forge::remodel::RemodelData::updateWithDsspAssignment(ObjexxFCL::FArray1D_char & dsspSS){
	for (int i = 0; i < (int)ss.size(); i++){
		int idx = this->blueprint[i].original_index;
		char const * ss_chars = ss.c_str();
		if (ss_chars[i] != '.'){
			dssp_updated_ss.append(1, ss_chars[i]);
		}
		else {
			dssp_updated_ss.append(1, dsspSS(idx));
		}
	}
	//turn upper case if not already so
	std::transform(dssp_updated_ss.begin(), dssp_updated_ss.end(), dssp_updated_ss.begin(), ::toupper);
	TR_REMODEL << "dssp_updated_ss:" << std::endl << dssp_updated_ss << std::endl;
}

void
protocols::forge::remodel::RemodelData::collectInsertionPose(){

	core::import_pose::pose_from_pdb( insertPose, basic::options::option[basic::options::OptionKeys::remodel::domainFusion::insert_segment_from_pdb]());
	insertionSize = (int)insertPose.total_residue();
	protocols::jumping::Dssp dssp(insertPose);
	ObjexxFCL::FArray1D_char dsspSS((int)insertPose.total_residue());
	dssp.dssp_reduced(dsspSS);
	for (int i = 1; i <= (int)dsspSS.size(); i++){
		insertionSS.push_back(dsspSS(i));
	}
	std::cout << "insertion SS: " << insertionSS << std::endl;
}

} //namespace remodel
} //namespace forge
} //namespace protocols
