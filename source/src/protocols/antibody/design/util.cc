// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/util.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/util.hh>
#include <boost/algorithm/string.hpp>
#include <string>
#include <utility/exit.hh>
#include <iostream>

namespace protocols {
namespace antibody {
namespace design {
	using namespace protocols::antibody;
	using namespace utility;
	
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


} //design
} //antibody
} //protocols
