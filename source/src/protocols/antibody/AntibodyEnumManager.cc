// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/AntibodyEnumManager.cc
/// @brief Functions for AntibodyEnumerators
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <map>
#include <string>

#include <utility/vector1.hh>
#include <utility/py/PyAssert.hh>

namespace protocols {
namespace antibody {


AntibodyEnumManager::AntibodyEnumManager() {
	setup();
}

AntibodyEnumManager::~AntibodyEnumManager() = default;

void
AntibodyEnumManager::setup() {

	///Initialize the length of the vectors
	cdr_name_to_string_.resize(CDRNameEnum_proto_total);
	numbering_scheme_to_string_.resize(AntibodyNumberingSchemeEnum_total);
	cdr_definition_to_string_.resize(CDRDefinitionEnum_total);
	light_chain_type_to_string_.resize(LightChainTypeEnum_total);
	h3_base_type_to_string_.resize(H3BaseTypeEnum_total);
	packing_angle_to_string_.resize(PackingAngleEnum_total);
	cdr_landmark_to_string_.resize(CDRLandmarkEnum_total);
	antibody_region_to_string_.resize(AntibodyRegionEnum_total);

	///Manually construct the variables

	////////////////// CDR names ///////////////////////////////////////////////
	cdr_name_to_string_[h1] = "H1";
	cdr_name_to_string_[h2] = "H2";
	cdr_name_to_string_[h3] = "H3";
	cdr_name_to_string_[l1] = "L1";
	cdr_name_to_string_[l2] = "L2";
	cdr_name_to_string_[l3] = "L3";

	cdr_name_to_string_[proto_l4] = "Proto_L4";
	cdr_name_to_string_[proto_h4] = "Proto_H4";

	cdr_name_to_enum_["H1"] = h1;
	cdr_name_to_enum_["H2"] = h2;
	cdr_name_to_enum_["H3"] = h3;
	cdr_name_to_enum_["L1"] = l1;
	cdr_name_to_enum_["L2"] = l2;
	cdr_name_to_enum_["L3"] = l3;

	cdr_name_to_enum_["L4"] = proto_l4;
	cdr_name_to_enum_["H4"] = proto_h4;
	cdr_name_to_enum_["PROTO_L4"] = proto_l4;
	cdr_name_to_enum_["PROTO_H4"] = proto_h4;

	///Just in case these are used in files.
	cdr_name_to_enum_["h1"] = h1;
	cdr_name_to_enum_["h2"] = h2;
	cdr_name_to_enum_["h3"] = h3;
	cdr_name_to_enum_["l1"] = l1;
	cdr_name_to_enum_["l2"] = l2;
	cdr_name_to_enum_["l3"] = l3;

	cdr_name_to_enum_["l4"] = l4;
	cdr_name_to_enum_["h4"] = h4;
	cdr_name_to_enum_["proto_l4"] = proto_l4;
	cdr_name_to_enum_["proto_h4"] = proto_h4;
	cdr_name_to_enum_["Proto_L4"] = proto_l4;
	cdr_name_to_enum_["Proto_H4"] = proto_h4;

	////////////////// Numbering Schemes ///////////////////////////////////////
	numbering_scheme_to_string_[Chothia_Scheme] = "Chothia_Scheme";
	numbering_scheme_to_string_[Kabat_Scheme] ="Kabat_Scheme";
	numbering_scheme_to_string_[Enhanced_Chothia_Scheme] = "Enhanced_Chothia_Scheme";
	numbering_scheme_to_string_[AHO_Scheme] = "AHO_Scheme";
	numbering_scheme_to_string_[IMGT_Scheme] = "IMGT_Scheme";

	numbering_scheme_to_enum_["Chothia_Scheme"] = Chothia_Scheme;
	numbering_scheme_to_enum_["Kabat_Scheme"] = Kabat_Scheme;
	numbering_scheme_to_enum_["Enhanced_Chothia_Scheme"] = Enhanced_Chothia_Scheme;
	numbering_scheme_to_enum_["AHO_Scheme"] = AHO_Scheme;
	numbering_scheme_to_enum_["IMGT_Scheme"] = IMGT_Scheme;

	////////////////// CDR Definitions  ///////////////////////////////////////
	cdr_definition_to_string_[Chothia] = "Chothia";
	cdr_definition_to_string_[Aroop] = "Aroop";
	cdr_definition_to_string_[North] = "North";
	cdr_definition_to_string_[Kabat] = "Kabat";
	cdr_definition_to_string_[Martin] = "Martin";

	cdr_definition_to_enum_["Chothia"] = Chothia;
	cdr_definition_to_enum_["Aroop"] = Aroop;
	cdr_definition_to_enum_["North"] = North;
	cdr_definition_to_enum_["Kabat"] = Kabat;
	cdr_definition_to_enum_["Martin"] = Martin;

	///////////////// LightChain Types ////////////////////////////////////////////
	light_chain_type_to_string_[lambda] = "lambda";
	light_chain_type_to_string_[kappa] = "kappa";
	light_chain_type_to_string_[lambda6] = "lambda6";
	light_chain_type_to_string_[unknown] = "unknown";

	light_chain_type_to_enum_["lambda"] = lambda;
	light_chain_type_to_enum_["kappa"] = kappa;
	light_chain_type_to_enum_["lambda6"] = lambda6;
	light_chain_type_to_enum_["unknown"] = unknown;

	////////////////// H3 Base Type ////////////////////////////////////////////
	h3_base_type_to_string_[Kinked] = "KINKED";
	h3_base_type_to_string_[Extended] = "EXTENDED";
	h3_base_type_to_string_[Neutral] = "NEUTRAL";
	h3_base_type_to_string_[Unknown] = "UNKNOWN";

	h3_base_type_to_enum_["KINKED"] = Kinked;
	h3_base_type_to_enum_["EXTENDED"] = Extended;
	h3_base_type_to_enum_["NEUTRAL"] = Neutral;
	h3_base_type_to_enum_["UNKNOWN"] = Unknown;

	///////////////// Packing Angle ////////////////////////////////////////////
	packing_angle_to_string_[VL_sheet_1] = "VL_sheet_1";
	packing_angle_to_string_[VL_sheet_2] = "VL_sheet_2";
	packing_angle_to_string_[VH_sheet_1] = "VH_sheet_1";
	packing_angle_to_string_[VH_sheet_2] = "VH_sheet_2";

	packing_angle_to_enum_["VL_sheet_1"] = VL_sheet_1;
	packing_angle_to_enum_["VL_sheet_2"] = VL_sheet_2;
	packing_angle_to_enum_["VH_sheet_1"] = VH_sheet_1;
	packing_angle_to_enum_["VH_sheet_2"] = VH_sheet_2;

	///////////////// Antibody Landmarks ////////////////////////////////////////////
	cdr_landmark_to_enum_["cdr_start"] = cdr_start;
	cdr_landmark_to_enum_["CDR_START"] =  cdr_start;
	cdr_landmark_to_enum_["cdr_stop"] = cdr_end;
	cdr_landmark_to_enum_["cdr_end"] = cdr_end;

	cdr_landmark_to_string_[cdr_start] = "cdr_start";
	cdr_landmark_to_string_[cdr_end] = "cdr_end";

	///////////////// Antibody Regions ////////////////////////////////////////////
	antibody_region_to_enum_["antigen_region"] = antigen_region;
	antibody_region_to_enum_["cdr_region"] = cdr_region;
	antibody_region_to_enum_["framework_region"] = framework_region;

	antibody_region_to_string_[antigen_region] = "antigen_region";
	antibody_region_to_string_[cdr_region] = "cdr_region";
	antibody_region_to_string_[framework_region] = "framework_region";

}

CDRNameEnum
AntibodyEnumManager::cdr_name_string_to_enum(std::string const & cdr_name) const {

	//This is here due to const correctness issues with [] operator
	auto iter( cdr_name_to_enum_.find( cdr_name ) );
	//utility::PyAssert((iter != cdr_name_to_enum_.end()), "CDR Name not found");
	return iter->second;
}

std::string
AntibodyEnumManager::cdr_name_enum_to_string(CDRNameEnum const cdr_name) const {
	return cdr_name_to_string_[cdr_name];
}

bool
AntibodyEnumManager::cdr_name_is_present(std::string const & cdr_name) const {

	auto iter( cdr_name_to_enum_.find( cdr_name ) );
	return iter != cdr_name_to_enum_.end();
}


////////////////// Numbering Schemes ///////////////////////////////////////
AntibodyNumberingSchemeEnum
AntibodyEnumManager::numbering_scheme_string_to_enum(std::string const & numbering_scheme) const {

	//This is here due to const correctness issues with [] operator
	auto iter( numbering_scheme_to_enum_.find( numbering_scheme ) );
	//utility::PyAssert((iter != numbering_scheme_to_enum_.end()), "Numbering not found");
	return iter->second;
}

std::string
AntibodyEnumManager::numbering_scheme_enum_to_string(AntibodyNumberingSchemeEnum const numbering_scheme) const {
	return numbering_scheme_to_string_[numbering_scheme];
}

bool
AntibodyEnumManager::numbering_scheme_is_present( std::string const & numbering_scheme) const {
	auto iter( numbering_scheme_to_enum_.find( numbering_scheme ) );
	return iter != numbering_scheme_to_enum_.end();
}


////////////////// CDR Definitions  ///////////////////////////////////////
CDRDefinitionEnum
AntibodyEnumManager::cdr_definition_string_to_enum(std::string const & cdr_definition) const {

	//This is here due to const correctness issues with [] operator
	auto iter( cdr_definition_to_enum_.find( cdr_definition ) );
	//utility::PyAssert((iter != numbering_scheme_to_enum_.end()), "Numbering not found");
	return iter->second;
}

std::string
AntibodyEnumManager::cdr_definition_enum_to_string(CDRDefinitionEnum const cdr_definition ) const {
	return cdr_definition_to_string_[cdr_definition];
}

bool
AntibodyEnumManager::cdr_definition_is_present(std::string const & cdr_definition) const {
	auto iter( cdr_definition_to_enum_.find( cdr_definition ) );
	return iter != cdr_definition_to_enum_.end();
}


///////////////// LightChain Types ////////////////////////////////////////////
LightChainTypeEnum
AntibodyEnumManager::light_chain_type_string_to_enum(const std::string& light_chain) const {
	auto iter( light_chain_type_to_enum_.find( light_chain ) );
	return iter->second;
}

std::string
AntibodyEnumManager::light_chain_type_enum_to_string(const LightChainTypeEnum light_chain) const {
	return light_chain_type_to_string_[light_chain];
}

bool
AntibodyEnumManager::light_chain_type_is_present(const std::string& light_chain) const {
	auto iter( light_chain_type_to_enum_.find( light_chain) );
	return iter != light_chain_type_to_enum_.end();
}


////////////////// H3 Base Type ////////////////////////////////////////////
H3BaseTypeEnum
AntibodyEnumManager::h3_base_type_string_to_enum(std::string const & base_type) const {

	//This is here due to const correctness issues with [] operator
	auto iter( h3_base_type_to_enum_.find( base_type ) );
	//utility::PyAssert((iter != h3_base_type_to_enum_.end()), "H3 base type not found");
	return iter->second;
}

std::string
AntibodyEnumManager::h3_base_type_enum_to_string(H3BaseTypeEnum const base_type) const {
	return h3_base_type_to_string_[base_type];
}


///////////////// Packing Angle ////////////////////////////////////////////
PackingAngleEnum
AntibodyEnumManager::packing_angle_string_to_enum(std::string const & angle_type) const {

	//This is here due to const correctness issues with [] operator
	auto iter( packing_angle_to_enum_.find( angle_type ) );
	//utility::PyAssert((iter != packing_angle_to_enum_.end()), "Packing Angle not found");
	return iter->second;
}

std::string
AntibodyEnumManager::packing_angle_enum_to_string(PackingAngleEnum const angle_type) const {
	return packing_angle_to_string_[angle_type];
}


///////////////// Antibody Landmarks ////////////////////////////////////////////
std::string
AntibodyEnumManager::cdr_landmark_enum_to_string(CDRLandmarkEnum const landmark) const {
	return cdr_landmark_to_string_[landmark];
}

CDRLandmarkEnum
AntibodyEnumManager::cdr_landmark_string_to_enum(std::string const &landmark) const {
	auto iter( cdr_landmark_to_enum_.find( landmark ) );
	return iter->second;
}

bool
AntibodyEnumManager::cdr_landmark_is_present(std::string const &landmark) const {
	auto iter( cdr_landmark_to_enum_.find( landmark ) );
	return iter != cdr_landmark_to_enum_.end();
}

///////////////// Antibody Regions ////////////////////////////////////////////
std::string
AntibodyEnumManager::antibody_region_enum_to_string(AntibodyRegionEnum const antibody_region) const {
	return antibody_region_to_string_[ antibody_region ];
}

AntibodyRegionEnum
AntibodyEnumManager::antibody_region_string_to_enum(std::string const & antibody_region) const {
	auto iter( antibody_region_to_enum_.find( antibody_region ) );
	return iter->second;
}

}
}

