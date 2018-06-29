// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/select/residue_selector/SSElementSelector.cc
/// @brief  The SSElementSelector selects based on the secondary element using DSSP.
/// @author TJ Brunette (tjbrunette@gmail.com)

// Unit headers
#include <core/select/residue_selector/SSElementSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/select/residue_selector/util.hh>

// Project headers
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

static basic::Tracer TR( "core.select.residue_selector.SSElementSelector" );

namespace core {
namespace select {
namespace residue_selector {


SSElementSelector::SSElementSelector() = default;
SSElementSelector::~SSElementSelector() = default;


SSElementSelector::SSElementSelector( SSElementSelector const &) = default;

ResidueSelectorOP SSElementSelector::clone() const { return ResidueSelectorOP( new SSElementSelector(*this) ); }

utility::vector1<SSElementSelector::SSElement> SSElementSelector::parse_ss(core::pose::Pose const & pose) const{
	utility::vector1<SSElementSelector::SSElement> ss_elements;
	core::scoring::dssp::Dssp dssp( pose );
	dssp.dssp_reduced();
	//initilize-------------------------------
	//default case --------
	Size start_pose_res = 1;
	Size end_pose_res = pose.size();
	//symmetric case  ------
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		end_pose_res = core::pose::symmetry::symmetry_info(pose)->num_independent_residues();
	}
	//chain entered case  ------
	bool skip_first_SS=false;
	if ( chain_!="" ) {
		//check chain existence
		if ( !has_chain(chain_[0],pose) ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "chain does not exist");
		}
		utility::vector1<core::Size> resnums = get_resnums_for_chain(pose, chain_[0]);
		start_pose_res=resnums[1];
		end_pose_res=resnums[resnums.size()];
		TR.Debug << "Chain " << chain_ << " selected; resi " << start_pose_res << "-" << end_pose_res << std::endl;
		if ( chain_ != "A" ) { //need to skip the first SS in the list later because of chain offset if not chain A
			skip_first_SS=true;
		}
	}
	//pass 1----------------------------------------
	char lastSecStruct = dssp.get_dssp_secstruct( start_pose_res );
	Size startSS = 1;
	for ( core::Size ii = start_pose_res+1; ii <= end_pose_res; ++ii ) {
		if ( dssp.get_dssp_secstruct(ii)!=lastSecStruct ) {
			if ( !skip_first_SS ) { //skip SS assignment for the first loop to compensate for chain offset
				Size endSS = ii-1;
				TR.Debug << "type: " << lastSecStruct << " startSS: " << startSS << " endSS: " << endSS << std::endl;
				SSElementSelector::SSElement ss_element_tmp(startSS,endSS,lastSecStruct);
				ss_elements.push_back(ss_element_tmp);
				startSS=ii;
			} else {
				TR.Debug << "Chain " << chain_ << " selected, first SS skipped (chain offset)." << std::endl;
				startSS=ii;
				skip_first_SS=false;
			}
		}
		lastSecStruct = dssp.get_dssp_secstruct(ii);
	}
	SSElementSelector::SSElement ss_element_tmp(startSS,end_pose_res,lastSecStruct);
	ss_elements.push_back(ss_element_tmp);
	//pass 2----------------------------------------Remove first and last loop. Could be done faster in above loop, but for future modifications I'm doing it below.
	if ( !reassign_short_terminal_loop_ ) {
		return(ss_elements);
	}

	bool delete_n_term_loop=false;
	bool delete_c_term_loop=false;
	Size ss_length = ss_elements[1].end_res-ss_elements[1].start_res;
	if ( (ss_length<reassign_short_terminal_loop_)&&(ss_elements[1].type=="L") ) {
		ss_elements[2].start_res =ss_elements[1].start_res;
		delete_n_term_loop=true;
	}
	ss_length = ss_elements[ss_elements.size()].end_res-ss_elements[ss_elements.size()].start_res;
	if ( (ss_length<reassign_short_terminal_loop_)&&(ss_elements[ss_elements.size()].type=="L") ) {
		ss_elements[ss_elements.size()-1].end_res =ss_elements[ss_elements.size()].end_res;
		delete_c_term_loop=true;
	}
	Size start_pos = 1;
	if ( delete_n_term_loop ) {
		start_pos=reassign_short_terminal_loop_;
	}
	Size end_pos = ss_elements.size();
	if ( delete_c_term_loop ) {
		end_pos = ss_elements.size()-1;
	}
	utility::vector1<SSElementSelector::SSElement> ss_elements_v2;
	for ( Size ii=start_pos; ii<=end_pos; ++ii ) {
		ss_elements_v2.push_back(ss_elements[ii]);
	}
	return(ss_elements_v2);
}

SSElementSelector::SSElement SSElementSelector::get_SSElement(utility::vector1<SSElementSelector::SSElement> ss_elements, int goal_position, std::string type,std::string description) const {
	int ss_ct=0;
	Size location=0;
	if ( goal_position>0 ) {
		for ( Size ii=1; ii<=ss_elements.size(); ++ii ) {
			if ( ss_elements[ii].type == type ) {
				ss_ct++;
				if ( goal_position==ss_ct ) {
					location=ii;
				}
			}
		}
	} else {
		for ( Size ii=ss_elements.size(); ii>=1; --ii ) {
			if ( ss_elements[ii].type == type ) {
				ss_ct++;
				if ( std::abs(goal_position)==ss_ct ) {
					location=ii;
				}
			}
		}
	}
	if ( location==0 ) {
		std::string error = "There are no positions corresponding to " + description + " check your SS counts";
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, error);
	}
	return(ss_elements[location]);
}

utility::vector1<Size> SSElementSelector::convert_string_to_residues(utility::vector1<SSElement> ss_elements, std::string description) const {
	//"2,L,M" middle of second loop
	//"2,L,S" start of second loop
	//"2,L"
	//"-2,L"
	//"n_term"
	//"c_term"
	//"middle"
	//n-term case------------
	utility::vector1<Size> selected_residues;
	if ( description =="n_term" ) {
		selected_residues.push_back(ss_elements[1].start_res);
		return(selected_residues);
	}
	//c-term case------------
	if ( description == "c_term" ) {
		selected_residues.push_back(ss_elements[ss_elements.size()].end_res);
		return(selected_residues);
	}
	//middle case------------
	if ( description == "middle" ) {
		Size n_term =ss_elements[1].start_res;
		Size c_term =ss_elements[ss_elements.size()].end_res;
		Size middle =(c_term-n_term)/2+n_term;
		selected_residues.push_back(middle);
	}
	//"2,L,M" and "2,L" case------------
	utility::vector1 <std::string> split_string = utility::string_split(description,',',std::string());
	int goal_position =utility::string2int(split_string[1]);
	std::string ss_type = split_string[2];
	SSElement ss_element_selected = get_SSElement(ss_elements,goal_position,ss_type,description);
	//"2,L" case
	if ( split_string.size()==2 ) {
		for ( Size ii=ss_element_selected.start_res; ii<=ss_element_selected.end_res; ++ii ) {
			selected_residues.push_back(ii);
		}
		return(selected_residues);
	}
	//"2,L,S", "2,L,E", "2,L,M" cases
	if ( split_string.size()==3 ) {
		std::string descriptor = split_string[3];
		Size middle = (ss_element_selected.end_res-ss_element_selected.start_res)/2+ss_element_selected.start_res;
		if ( descriptor =="S" ) {
			selected_residues.push_back(ss_element_selected.start_res);
		}
		if ( descriptor == "E" ) {
			selected_residues.push_back(ss_element_selected.end_res);
		}
		if ( descriptor == "M" ) {
			selected_residues.push_back(middle);
		}
		return(selected_residues);
	}
	throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Something is wrong with" + description);
	return(selected_residues);
}
ResidueSubset SSElementSelector::combine_residue_selections(utility::vector1<Size> start_selection_residues, utility::vector1<Size>  end_selection_residues,core::pose::Pose const & pose) const{
	ResidueSubset subset( pose.size(), false );
	Size nTerm_start=0;
	Size cTerm_start=0;
	Size nTerm_end=0;
	Size cTerm_end=0;
	nTerm_start = start_selection_residues[1];
	for ( Size ii=1; ii<=start_selection_residues.size(); ++ii ) {
		subset[start_selection_residues[ii]]=true;
		cTerm_start = start_selection_residues[ii]; //can assume the vector is ordered
	}
	if ( end_selection_residues.size()>0 ) {
		nTerm_end= end_selection_residues[1]; // also assume ordered vector
		for ( Size ii=1; ii<=end_selection_residues.size(); ++ii ) {
			subset[end_selection_residues[ii]]=true;
			cTerm_end = end_selection_residues[ii]; //can assume the vector is ordered
		}
		//positive indeces
		for ( Size jj=cTerm_start+1; jj<nTerm_end; ++jj ) {
			subset[jj]=true;
		}
		//neg indeces
		for ( Size jj=cTerm_end+1; jj<nTerm_start; ++jj ) {
			subset[jj]=true;
		}
	}
	return(subset);
}



ResidueSubset
SSElementSelector::apply( core::pose::Pose const & pose ) const
{
	utility::vector1<SSElement> ss_elements = parse_ss(pose);
	utility::vector1<Size> start_selection_residues = convert_string_to_residues(ss_elements,start_);
	utility::vector1<Size> end_selection_residues;
	if ( end_ != "" ) {
		end_selection_residues = convert_string_to_residues(ss_elements,end_);
	}
	ResidueSubset subset = combine_residue_selections(start_selection_residues,end_selection_residues,pose);
	return subset;
}

void SSElementSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &
)
{
	start_ = tag->getOption< std::string >("selection");
	end_ = tag->getOption< std::string >("to_selection","");//switched naming because it can get confusing when using negative indices
	chain_ = tag->getOption<std::string >("chain","");
	reassign_short_terminal_loop_ = tag->getOption< core::Size >("reassign_short_terminal_loop",2);
}


std::string SSElementSelector::get_name() const {
	return SSElementSelector::class_name();
}

std::string SSElementSelector::class_name() {
	return "SSElement";
}

void
SSElementSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes + XMLSchemaAttribute( "selection", xs_string , "H=helix,L=Loop,S=Sheet,N=n_terminal,C=terminal can define only start" )
		+ XMLSchemaAttribute::attribute_w_default( "to_selection", xs_string , "H=helix,L=Loop,S=Sheet,N=n_terminal,C=terminal","" )
		+ XMLSchemaAttribute::attribute_w_default( "chain", xs_string , "chain letter","" )
		+ XMLSchemaAttribute::attribute_w_default( "reassign_short_terminal_loop", xsct_non_negative_integer , "if terminal less than X residues loop is reassigned to neighboring SS element (default=2; 0=no reassignment).","2");
	xsd_type_definition_w_attributes( xsd, class_name(), "a selector for choosing parts of secondary structure", attributes );
}


ResidueSelectorOP
SSElementSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new SSElementSelector );
}

std::string
SSElementSelectorCreator::keyname() const {
	return SSElementSelector::class_name();
}

void
SSElementSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	SSElementSelector::provide_xml_schema( xsd );
}


} //namespace residue_selector
} //namespace select
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::select::residue_selector::SSElementSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( start_ ) ); // std::string
	arc( CEREAL_NVP( end_ ) ); // std::string
	arc( CEREAL_NVP( chain_) ); //std::string
	arc( CEREAL_NVP(reassign_short_terminal_loop_) ); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::select::residue_selector::SSElementSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( start_ ); // std::string
	arc( end_ ); // std::string
	arc( chain_ ); //std::string
	arc( reassign_short_terminal_loop_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::SSElementSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::SSElementSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_SSElementSelector )
#endif // SERIALIZATION
