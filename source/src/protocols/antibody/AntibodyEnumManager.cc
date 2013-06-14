// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/AntibodyEnumManager.cc
/// @brief Functions for AntibodyEnumerators
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <map>
#include <string>

#include <utility/vector1.hh>
#include <utility/PyAssert.hh>

namespace protocols {
namespace antibody {


AntibodyEnumManager::AntibodyEnumManager() {
	setup();
}

AntibodyEnumManager::~AntibodyEnumManager() {}

void
AntibodyEnumManager::setup() {

	///Initialize the length of the vectors
	cdr_name_to_string_.resize(num_cdr_loops);
	numbering_scheme_to_string_.resize(AntibodyNumberingSchemeEnum_total);
	h3_base_type_to_string_.resize(H3BaseTypeEnum_total);
	packing_angle_to_string_.resize(PackingAngleEnum_total);

	///Manually construct the variables

	cdr_name_to_string_[h1] = "H1";
	cdr_name_to_string_[h2] = "H2";
	cdr_name_to_string_[h3] = "H3";
	cdr_name_to_string_[l1] = "L1";
	cdr_name_to_string_[l2] = "L2";
	cdr_name_to_string_[l3] = "L3";

	cdr_name_to_enum_["H1"] = h1;
	cdr_name_to_enum_["H2"] = h2;
	cdr_name_to_enum_["H3"] = h3;
	cdr_name_to_enum_["L1"] = l1;
	cdr_name_to_enum_["L2"] = l2;
	cdr_name_to_enum_["L3"] = l3;

	///Just in case these are used in files.
	cdr_name_to_enum_["h1"] = h1;
	cdr_name_to_enum_["h2"] = h2;
	cdr_name_to_enum_["h3"] = h3;
	cdr_name_to_enum_["l1"] = l1;
	cdr_name_to_enum_["l2"] = l2;
	cdr_name_to_enum_["l3"] =l3;

	numbering_scheme_to_string_[Aroop] = "Aroop";
	numbering_scheme_to_string_[Chothia] = "Chothia";
	numbering_scheme_to_string_[Kabat] ="Kabat";
	numbering_scheme_to_string_[Enhanced_Chothia] = "Enhanced_Chothia";
	numbering_scheme_to_string_[AHO] = "AHO";
	numbering_scheme_to_string_[Modified_AHO] = "Modified_AHO";
	numbering_scheme_to_string_[IMGT] = "IMGT";

	numbering_scheme_to_enum_["Aroop"] = Aroop;
	numbering_scheme_to_enum_["Chothia"] = Chothia;
	numbering_scheme_to_enum_["Kabat"] = Kabat;
	numbering_scheme_to_enum_["Enhanced_Chothia"] = Enhanced_Chothia;
	numbering_scheme_to_enum_["AHO"] = AHO;
	numbering_scheme_to_enum_["Modified_AHO"] = Modified_AHO;
	numbering_scheme_to_enum_["IMGT"] = IMGT;

	h3_base_type_to_string_[Kinked] = "KINKED";
	h3_base_type_to_string_[Extended] = "EXTENDED";
	h3_base_type_to_string_[Neutral] = "NEUTRAL";
	h3_base_type_to_string_[Unknown] = "UNKNOWN";

	h3_base_type_to_enum_["KINKED"] = Kinked;
	h3_base_type_to_enum_["EXTENDED"] = Extended;
	h3_base_type_to_enum_["NEUTRAL"] = Neutral;
	h3_base_type_to_enum_["UNKNOWN"] = Unknown;

	packing_angle_to_string_[VL_sheet_1] = "VL_sheet_1";
	packing_angle_to_string_[VL_sheet_2] = "VL_sheet_2";
	packing_angle_to_string_[VH_sheet_1] = "VH_sheet_1";
	packing_angle_to_string_[VH_sheet_2] = "VH_sheet_2";

	packing_angle_to_enum_["VL_sheet_1"] = VL_sheet_1;
	packing_angle_to_enum_["VL_sheet_2"] = VL_sheet_2;
	packing_angle_to_enum_["VH_sheet_1"] = VH_sheet_1;
	packing_angle_to_enum_["VH_sheet_2"] = VH_sheet_2;

}

CDRNameEnum
AntibodyEnumManager::cdr_name_string_to_enum(std::string const & cdr_name) const {

	//This is here due to const correctness issues with [] operator
	std::map< std::string, CDRNameEnum >::const_iterator iter( cdr_name_to_enum_.find( cdr_name ) );
	//utility::PyAssert((iter != cdr_name_to_enum_.end()), "CDR Name not found");
	return iter->second;
}

std::string
AntibodyEnumManager::cdr_name_enum_to_string(CDRNameEnum const cdr_name) const {
	return cdr_name_to_string_[cdr_name];
}

bool
AntibodyEnumManager::cdr_name_is_present(std::string const & cdr_name) const {

	std::map< std::string, CDRNameEnum >::const_iterator iter( cdr_name_to_enum_.find( cdr_name ) );
	return iter != cdr_name_to_enum_.end();
}

AntibodyNumberingSchemeEnum
AntibodyEnumManager::numbering_scheme_string_to_enum(std::string const & numbering_scheme) const {

	//This is here due to const correctness issues with [] operator
	std::map< std::string, AntibodyNumberingSchemeEnum >::const_iterator iter( numbering_scheme_to_enum_.find( numbering_scheme ) );
	//utility::PyAssert((iter != numbering_scheme_to_enum_.end()), "Numbering not found");
	return iter->second;
}

std::string
AntibodyEnumManager::numbering_scheme_enum_to_string(AntibodyNumberingSchemeEnum const numbering_scheme) const {
	return numbering_scheme_to_string_[numbering_scheme];
}

H3BaseTypeEnum
AntibodyEnumManager::h3_base_type_string_to_enum(std::string const & base_type) const {

	//This is here due to const correctness issues with [] operator
	std::map< std::string, H3BaseTypeEnum >::const_iterator iter( h3_base_type_to_enum_.find( base_type ) );
	//utility::PyAssert((iter != h3_base_type_to_enum_.end()), "H3 base type not found");
	return iter->second;
}

std::string
AntibodyEnumManager::h3_base_type_enum_to_string(H3BaseTypeEnum const base_type) const {
	return h3_base_type_to_string_[base_type];
}

PackingAngleEnum
AntibodyEnumManager::packing_angle_string_to_enum(std::string const & angle_type) const {

	//This is here due to const correctness issues with [] operator
	std::map< std::string, PackingAngleEnum >::const_iterator iter( packing_angle_to_enum_.find( angle_type ) );
	//utility::PyAssert((iter != packing_angle_to_enum_.end()), "Packing Angle not found");
	return iter->second;
}

std::string
AntibodyEnumManager::packing_angle_enum_to_string(PackingAngleEnum const angle_type) const {
	return packing_angle_to_string_[angle_type];
}

}
}

