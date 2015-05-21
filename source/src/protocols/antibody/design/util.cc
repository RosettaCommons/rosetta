// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/design/util.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/database/CDRSetOptions.hh>
#include <protocols/antibody/database/CDRSetOptionsParser.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>
#include <protocols/antibody/design/CDRGraftDesignOptions.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <protocols/grafting/util.hh>
#include <protocols/loops/util.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

#include <boost/algorithm/string.hpp>
#include <string>
#include <utility/exit.hh>
#include <iostream>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>

static basic::Tracer TR("protocols.antibody.design.util");

namespace protocols {
namespace antibody {
namespace design {
	using namespace core::pack::task::operation;
	using namespace protocols::antibody;
	using namespace utility;
	

void
insert_cdr_into_antibody(AntibodyInfoCOP ab_info, CDRNameEnum const cdr, core::pose::Pose & pose, core::pose::Pose & cdr_piece, core::Size overhang) {
	core::Size cdr_start = ab_info->get_CDR_start(cdr, pose);
	core::Size cdr_end = ab_info->get_CDR_end(cdr, pose);
	
	protocols::grafting::delete_region(cdr_piece, 1, overhang);
	protocols::grafting::delete_region(cdr_piece, cdr_piece.total_residue() - overhang + 1, cdr_piece.total_residue());
	
	//core::Size insert_length = cdr_piece.total_residue();
	protocols::grafting::delete_region(pose, cdr_start+1, cdr_end - 1);
	pose = protocols::grafting::insert_pose_into_pose(pose, cdr_piece, cdr_start, cdr_start+1);

	pose.pdb_info()->copy(*(cdr_piece.pdb_info()), 1, cdr_piece.total_residue(), cdr_start+1);
	pose.pdb_info()->obsolete(false);
	
}


void
get_all_graft_permutations(
		vector1<core::Size > & cdr_set_totals,
		vector1< vector1< core::Size> > & all_permutations,
		vector1<core::Size>current_index,
		core::Size const cdr_num
)
{
	//Current index is what is being worked on. 
	
	if (cdr_num > cdr_set_totals.size()){return;}
	
	//No CDR in CDR set.  Set index 0, move on to next CDR
	if (cdr_set_totals[cdr_num]==0){
		current_index[cdr_num]=0;
		//Inner most loop, add current_index and return;
		if (cdr_num==cdr_set_totals.size()){
			all_permutations.push_back(current_index);
			return;
		}
		get_all_graft_permutations(cdr_set_totals, all_permutations, current_index, cdr_num+1);
	}
	else{
		for(core::Size i = 1; i<=cdr_set_totals[cdr_num]; ++i){
			current_index[cdr_num]=i;
			//Inner most loop, add current_index and return;
			if (cdr_num==cdr_set_totals.size()){
				all_permutations.push_back(current_index);
				return;
			}
			get_all_graft_permutations(cdr_set_totals, all_permutations, current_index, cdr_num+1);
		
		}
	}
}

DesignTypeEnum
design_type_from_string(std::string const design_type){
	std::string type = design_type;
	boost::to_upper(type);
	
	if (type == "FLXBB"){
		return flxbb;
	}
	else if(type == "FIXBB" || type == "FIXEDBB"){
		return fixbb;
	}
	else if(type == "RELAXED_DESIGN" || type == "RELAX_DESIGN" || type == "RELAX_TF"){
		return relaxed_design;
	}
	else{
		utility_exit_with_message("DesignType unrecognized.  Please check AntibodyDesign settings.");
	}
}

std::string
design_type_to_string(DesignTypeEnum const design_type){
	if (design_type == flxbb){
		return "FLXBB";
	}
	else if(design_type==fixbb){
		return "FIXBB";
	}
	else if(design_type==relaxed_design){
		return "RELAXED_DESIGN";
	}
	else{
		utility_exit_with_message("DesignType unrecognized.  Please check AntibodyDesign settings.");
	}
}

std::string
get_dock_chains_from_ab_dock_chains(AntibodyInfoCOP ab_info, std::string ab_dock_chains){
	vector1<char> chains;
	for (core::Size i = 0; i <= ab_dock_chains.length()-1; ++i){
		char ab_chain = ab_dock_chains[i];
		if (ab_chain == 'A'){
			vector1<char>antigen = ab_info->get_antigen_chains();
			if (antigen.size() == 0){
				TR <<" Antigen not present to dock. skipping addition of Antigen. - setting L_H dock instead" << std::endl;
				return "L_H";
				
			}
			
			for (core::Size x = 1; x <= antigen.size(); ++x){
				chains.push_back(antigen[x]);
			}
		}
		else if (ab_chain == 'L' || ab_chain == 'H' || ab_chain == '_'){
			chains.push_back(ab_chain);
		}
		else{
			utility_exit_with_message("ab_dock_chains must be L H or A: " + ab_dock_chains);
		}
	}
	std::string dock_chains(chains.begin(), chains.end());
	return dock_chains;
}

vector1<PDBNumbering>
get_pdb_numbering_from_string(vector1<std::string> const & pdb_residues) {
	vector1<PDBNumbering> numbering;
	
	for (core::Size i = 1; i <= pdb_residues.size(); ++i){
		vector1<std::string> res_str = utility::string_split(pdb_residues[i], ':');
		PDBNumbering number;
		if (res_str.size()== 1){
			assert(res_str[1].length() >= 2);
			
			number.icode = ' ';
			number.chain = res_str[1][res_str[1].length() - 1];
			std::stringstream(res_str[1].substr(0, res_str[1].length() - 1)) >> number.resnum;
			
		}
		else if (res_str.size() == 2){
			assert(res_str[1].length() >= 2);
			assert(res_str[2].length() == 1);
			
			number.icode = res_str[2][0];
			number.chain = res_str[1][res_str[1].length() - 1];
			std::stringstream(res_str[1].substr(0, res_str[1].length() - 1)) >> number.resnum;
		}
		else {
			utility_exit_with_message("Cannot convert string to pdb string: "+pdb_residues[i]);
		}
		numbering.push_back(number);
	}
	return numbering;
}

vector1<bool>
get_resnum_from_pdb_numbering(core::pose::Pose const & pose, vector1<PDBNumbering> const & pdb_residues){
	
	vector1<bool> residues(pose.total_residue(), false);
	for (core::Size i = 1; i <= pdb_residues.size(); ++i){
		PDBNumbering numbering = pdb_residues[i];
		core::Size resnum = pose.pdb_info()->pdb2pose(numbering.chain, numbering.resnum, numbering.icode);
		residues[resnum] = true;
	}
	return residues;
}

std::pair<bool, core::Size>
check_cb(core::pose::Pose const & pose, utility::vector1<bool> const & residues){
	assert(residues.size() == pose.total_residue());
	
	std::pair<bool, core::Size> cb = std::make_pair(false, 0);
	for (core::Size i = 1; i <= residues.size(); ++i){
		if (! residues[i]) continue;
		cb = protocols::loops::has_severe_pep_bond_geom_issues(pose, i, true, true, 1.5, 15, 15);
		if (cb.first) {
			return cb;
		}
	}
	
	return cb;
}

std::pair<bool, core::Size>
check_cb(core::pose::Pose const & pose, protocols::loops::Loops const & loops){
	std::pair<bool, core::Size> cb = std::make_pair(false, 0);
	for (protocols::loops::Loops::const_iterator it = loops.begin(); it != loops.end(); ++it){
		cb = protocols::loops::has_severe_pep_bond_geom_issues(pose, it->cut(), true, true, 1.5, 15, 15);
		if (cb.first == true){
			return cb;
		}
	}
	return cb;
}


core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_region(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose,
	AntibodyRegionEnum region)
{
	RestrictResidueToRepackingOP restrict( new RestrictResidueToRepacking() );
	for (core::Size i = 1; i <= pose.total_residue(); ++i){
		if (ab_info->get_region_of_residue(pose, i) == region){
			restrict->include_residue(i);
		}
	}
	return restrict;
}

core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_antigen(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose)
{
	return disable_design_region(ab_info, pose, antigen_region);
}

core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_framework(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose) 
{
	return disable_design_region(ab_info, pose, framework_region);
}

core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_cdrs(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose)
{
	return disable_design_region(ab_info, pose, cdr_region);
}

core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_cdr(
	AntibodyInfoCOP ab_info,
	CDRNameEnum cdr,
	const core::pose::Pose & pose) {

	//One restrict op per CDR.  That way we can pop them off the TF  individually if we need to.
	RestrictResidueToRepackingOP restrict( new RestrictResidueToRepacking() );
	core::Size start = ab_info->get_CDR_start(cdr, pose);
	core::Size end = ab_info->get_CDR_end(cdr, pose);
	for (core::Size i = start; i <= end; ++i){
		restrict->include_residue(i);
	}
	return restrict;
}

core::pack::task::operation::RestrictResidueToRepackingOP
disable_conserved_framework_positions(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose)
{
	RestrictResidueToRepackingOP restrict( new RestrictResidueToRepacking() );
	utility::vector1< core::Size > conserved_positions;
	
	if (core::pose::has_chain( "L", pose )) {
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'L', 23, ' ', false ) );
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'L', 43, ' ', false ) );
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'L', 106, ' ', false ) );
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'L', 139, ' ', false ) );
	}
	
	if (core::pose::has_chain( "H", pose )) {
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'H', 23, ' ', false ) );
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'H', 43, ' ', false ) );
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'H', 106, ' ', false ) );
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'H', 139, ' ', false ) );
	}
	
	for ( core::Size i = 1; i <= conserved_positions.size(); ++i ){
		if ( conserved_positions[ i ] != 0 ) {
			restrict->include_residue( conserved_positions[ i ] );
		}
	}
	return restrict;
}

AntibodyCDRSetOptions
get_cdr_set_options(){
	CDRSetOptionsParser parser = CDRSetOptionsParser();
	AntibodyCDRSetOptions options_settings;
	
	if (basic::options::option [basic::options::OptionKeys::antibody::design::instructions].user()){
		std::string filename = basic::options::option [basic::options::OptionKeys::antibody::design::instructions]();
		for (core::Size i = 1; i <= 6; ++i){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			options_settings.push_back(parser.parse_default_and_user_options(cdr, filename));
		}
	}
	else{
		std::string filename = basic::options::option [basic::options::OptionKeys::antibody::design::base_instructions]();
		for (core::Size i = 1; i <= 6; ++i){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			options_settings.push_back(parser.parse_options(cdr, filename));
		}
	}
	return options_settings;
	
}

AntibodyCDRGraftDesignOptions
get_graft_design_options(){
	
	CDRGraftDesignOptionsParser parser = CDRGraftDesignOptionsParser();
	AntibodyCDRGraftDesignOptions options_settings;
	
	if (basic::options::option [basic::options::OptionKeys::antibody::design::instructions].user()){
		std::string filename = basic::options::option [basic::options::OptionKeys::antibody::design::instructions]();
		for (core::Size i = 1; i <= 6; ++i){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			options_settings.push_back(parser.parse_default_and_user_options(cdr, filename));
		}
	}
	else {
		std::string filename = basic::options::option [basic::options::OptionKeys::antibody::design::base_instructions]();
		for (core::Size i = 1; i <= 6; ++i){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			options_settings.push_back(parser.parse_options(cdr, filename));
		}
	}
	

	return options_settings;
}

AntibodyCDRSeqDesignOptions
get_seq_design_options(){
	CDRSeqDesignOptionsParser parser = CDRSeqDesignOptionsParser();
	AntibodyCDRSeqDesignOptions options_settings;
	
	if ( basic::options::option [basic::options::OptionKeys::antibody::design::instructions].user() ){
		std::string filename = basic::options::option [basic::options::OptionKeys::antibody::design::instructions]();
		for (core::Size i = 1; i <=6; ++i){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			options_settings.push_back(parser.parse_default_and_user_options(cdr, filename));
		}
	}
	else {
		std::string filename = basic::options::option [basic::options::OptionKeys::antibody::design::base_instructions]();
		for (core::Size i = 1; i <=6; ++i){
			CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
			options_settings.push_back(parser.parse_options(cdr, filename));
		}
	}
	return options_settings;
}

} //design
} //antibody
} //protocols
