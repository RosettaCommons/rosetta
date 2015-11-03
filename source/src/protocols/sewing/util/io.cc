// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file util.cc
///
/// @brief
/// @author Tim Jacobs

//Unit
#include <protocols/sewing/util/io.hh>

//Core headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>

//Utility headers
#include <basic/Tracer.hh>

#include <numeric/xyzVector.hh>

#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

#include <boost/foreach.hpp>

namespace protocols {
namespace sewing  {

static basic::Tracer TR( "protocols.sewing.io" );

void
write_native_residue_file(
	NativeRotamersMap native_residue_map,
	std::string filename
){
	utility::io::ozstream native_residue_file;
	native_residue_file.open(filename);

	for ( NativeRotamersMap::const_iterator pos_it=native_residue_map.begin();
			pos_it != native_residue_map.end(); ++pos_it ) {
		native_residue_file << "RESNUM " << pos_it->first << std::endl;
		utility::vector1< std::pair<bool, core::conformation::ResidueOP> >::const_iterator res_it = pos_it->second.begin();
		utility::vector1< std::pair<bool, core::conformation::ResidueOP> >::const_iterator res_it_end = pos_it->second.end();
		for ( ; res_it != res_it_end; ++res_it ) {
			native_residue_file << "RESIDUE " << res_it->second->type().name() << " " << res_it->first << std::endl;
			if ( res_it->first ) {
				for ( core::Size i=1; i<=res_it->second->natoms(); ++i ) {
					native_residue_file << "ATOM " << i << " " <<
						res_it->second->xyz(i).x() << " " <<
						res_it->second->xyz(i).y() << " " <<
						res_it->second->xyz(i).z() << " " <<
						std::endl;
				}
			}
		}
	}
	native_residue_file.close();
}//write_native_residue_file

NativeRotamersMap
read_native_residue_file(
	std::string filename
) {
	utility::io::izstream file(filename);
	if ( !file.good() ) {
		utility_exit_with_message("Could not find NativeRotamersFile file with name: " + filename);
	}

	core::chemical::ResidueTypeSetCOP res_type_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

	NativeRotamersMap rotamers;
	std::string line;

	core::Size cur_resnum = 0;
	core::conformation::ResidueOP cur_residue = 0;
	bool save_chis = false;
	while ( getline( file, line ) ) {
		utility::vector1<std::string> tokens = utility::string_split(line);
		assert(tokens.size() > 0);

		if ( tokens[1]=="RESNUM" ) {
			if ( cur_resnum != 0 ) {
				rotamers[cur_resnum].push_back(std::make_pair(save_chis, cur_residue));
			}
			std::istringstream(tokens[2]) >> cur_resnum;
		} else if ( tokens[1]=="RESIDUE" ) {
			if ( cur_residue != 0 ) {
				rotamers[cur_resnum].push_back(std::make_pair(save_chis, cur_residue));
			}
			core::chemical::ResidueType const & res_type = res_type_set->name_map(tokens[2]);
			cur_residue = core::conformation::ResidueFactory::create_residue(res_type);
			save_chis = utility::string2int(tokens[3]);
		} else if ( tokens[1]=="ATOM" ) {
			core::Size atom_index;
			std::istringstream(tokens[2]) >> atom_index;

			core::Real x,y,z;
			std::istringstream(tokens[3]) >> x;
			std::istringstream(tokens[4]) >> y;
			std::istringstream(tokens[5]) >> z;
			numeric::xyzVector<core::Real> xyz(x, y, z);
			cur_residue->set_xyz(atom_index, xyz);
		}
	}
	rotamers[cur_resnum].push_back(std::make_pair(save_chis, cur_residue));
	TR << "Read in Rotamers for " << rotamers.size() << " positions" << std::endl;
	return rotamers;
}


void
write_hashing_scores_to_file(
	ScoreResults const & scores,
	std::string filename
){
	utility::io::ozstream file;
	file.open_append(filename);

	for ( ScoreResults::const_iterator it = scores.begin(); it != scores.end(); ++it ) {
		BasisPair const & basis_residues = it->first;

		std::map< SegmentPair, core::Size > const & segment_matches = it->second.segment_match_counts;
		std::map< SegmentPair, core::Size >::const_iterator seg_it = segment_matches.begin();
		std::map< SegmentPair, core::Size >::const_iterator seg_it_end = segment_matches.end();
		core::Size sum=0;
		std::stringstream model_1_segments;
		std::stringstream model_2_segments;
		for ( ; seg_it != seg_it_end; ++seg_it ) {
			sum += seg_it->second;
			model_1_segments << seg_it->first.first << " ";
			model_2_segments << seg_it->first.second << " ";
		}
		core::Real average_segment_score = 0;
		if ( segment_matches.size() > 0 ) {
			average_segment_score = sum/segment_matches.size();
		}
		file
			<< basis_residues.first.model_id << " "
			<< basis_residues.first.resnum << " "
			<< model_1_segments.str()
			<< basis_residues.second.model_id << " "
			<< basis_residues.second.resnum << " "
			<< model_2_segments.str()
			<< average_segment_score << std::endl;
	}
	file.close();
}

/*
std::string
see_whether_model_is_H_bonded_by_terminal_strands(
Model model){
if ( !utility::file::file_exists( model.pdb_code_ ) ){
std::stringstream err;
err << "You must provide pdb files to make pose as defined in your input model file";
utility_exit_with_message(err.str());
}

std::string model_is_H_bonded_by_terminal_strands = "no"; // just initial assignment

// we need to deal with models whose both terminal segments are strands
if (model.segments_[1].dssp_ != 'E' || model.segments_[model.segments_.size()].dssp_ != 'E'){
if(TR.Debug.visible()){
TR << "(model.segments_[1].dssp_ != 'E' || model.segments_[model.segments_.size()].dssp_ != 'E')" << std::endl;
}
return model_is_H_bonded_by_terminal_strands;
}
// TR << "model.pdb_code_: " << model.pdb_code_ << std::endl;


core::pose::Pose pose;
core::import_pose::pose_from_pdb( pose, model.pdb_code_ );
if(TR.Debug.visible()){
TR << "pose.pdb_info()->name(): " << pose.pdb_info()->name() << std::endl;
TR << "model.segments_[1].residues_.front().resnum_: " << model.segments_[1].residues_.front().resnum_ << std::endl;
TR << "model.segments_[1].residues_.back().resnum_: " << model.segments_[1].residues_.back().resnum_ << std::endl;
TR << "model.segments_[model.segments_.size()].residues_.front().resnum_: " << model.segments_[model.segments_.size()].residues_.front().resnum_ << std::endl;
TR << "model.segments_[model.segments_.size()].residues_.back().resnum_: " << model.segments_[model.segments_.size()].residues_.back().resnum_ << std::endl;
}

protocols::features::strand_assembly::SandwichFragment strand_i(model.segments_[1].residues_.front().resnum_, model.segments_[1].residues_.back().resnum_);
protocols::features::strand_assembly::SandwichFragment strand_j(model.segments_[model.segments_.size()].residues_.front().resnum_, model.segments_[model.segments_.size()].residues_.back().resnum_);


core::Size return_of_find_sheet_anti = find_sheet(
pose,
strand_i,
strand_j,
1, //  find anti-parallel sheet
3.5, //  min_CA_CA_dis_,
6.2, //  max_CA_CA_dis_,
120 //  min_C_O_N_angle_
);

// TR << "return_of_find_sheet_anti: " << return_of_find_sheet_anti << std::endl;

if (return_of_find_sheet_anti == 1){
model_is_H_bonded_by_terminal_strands = "antiparallel";
return model_is_H_bonded_by_terminal_strands;
}

else {
core::Size return_of_find_sheet_para = find_sheet(
pose,
strand_i,
strand_j,
0, // find parallel sheet
3.5, //  min_CA_CA_dis_,
6.2, //  max_CA_CA_dis_,
120 //  min_C_O_N_angle_
);
//   TR << "return_of_find_sheet_para: " << return_of_find_sheet_para << std::endl;
if (return_of_find_sheet_para == 1){
model_is_H_bonded_by_terminal_strands = "parallel";
}
return model_is_H_bonded_by_terminal_strands;
}
}//see_whether_model_is_H_bonded_by_terminal_strands
*/

std::string
see_whether_model_is_H_bonded_by_terminal_strands(
	Model model,
	std::string P_PA){
	if ( !utility::file::file_exists( model.pdb_code_ ) ) {
		std::stringstream err;
		err << "You must provide pdb files to make pose as defined in your input model file";
		utility_exit_with_message(err.str());
	}

	std::string model_is_H_bonded_by_terminal_strands = "no"; // just initial assignment

	// we need to deal with models whose both terminal segments are strands
	if ( model.segments_[1].dssp_ != 'E' || model.segments_[model.segments_.size()].dssp_ != 'E' ) {
		if ( TR.Debug.visible() ) {
			TR << "(model.segments_[1].dssp_ != 'E' || model.segments_[model.segments_.size()].dssp_ != 'E')" << std::endl;
		}
		return model_is_H_bonded_by_terminal_strands;
	}

	core::pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, model.pdb_code_ );
	if ( TR.Debug.visible() ) {
		TR << "pose.pdb_info()->name(): " << pose.pdb_info()->name() << std::endl;
		TR << "model.segments_[1].residues_.front().resnum_: " << model.segments_[1].residues_.front().resnum_ << std::endl;
		TR << "model.segments_[1].residues_.back().resnum_: " << model.segments_[1].residues_.back().resnum_ << std::endl;
		TR << "model.segments_[model.segments_.size()].residues_.front().resnum_: " << model.segments_[model.segments_.size()].residues_.front().resnum_ << std::endl;
		TR << "model.segments_[model.segments_.size()].residues_.back().resnum_: " << model.segments_[model.segments_.size()].residues_.back().resnum_ << std::endl;
	}

	protocols::features::strand_assembly::SandwichFragment strand_i(model.segments_[1].residues_.front().resnum_, model.segments_[1].residues_.back().resnum_);
	protocols::features::strand_assembly::SandwichFragment strand_j(model.segments_[model.segments_.size()].residues_.front().resnum_, model.segments_[model.segments_.size()].residues_.back().resnum_);

	if ( P_PA == "antiparallel" ) {
		core::Size return_of_find_sheet_anti = find_sheet(
			pose,
			strand_i,
			strand_j,
			1, //  find anti-parallel sheet
			3.5, //  min_CA_CA_dis_,
			6.2, //  max_CA_CA_dis_,
			120, //  min_C_O_N_angle_
			false); //care_smaller_sheet
		if ( return_of_find_sheet_anti == 1 ) {
			model_is_H_bonded_by_terminal_strands = "antiparallel";
		}
		return model_is_H_bonded_by_terminal_strands;
	} else if ( P_PA == "parallel" ) {
		core::Size return_of_find_sheet_para = find_sheet(
			pose,
			strand_i,
			strand_j,
			0, // find parallel sheet
			3.5, //  min_CA_CA_dis_,
			6.2, //  max_CA_CA_dis_,
			120, //  min_C_O_N_angle_
			false); //care_smaller_sheet
		if ( return_of_find_sheet_para == 1 ) {
			model_is_H_bonded_by_terminal_strands = "parallel";
		}
		return model_is_H_bonded_by_terminal_strands;
	}
	return model_is_H_bonded_by_terminal_strands;
}//see_whether_model_is_H_bonded_by_terminal_strands

}//sewing
}//protocols
