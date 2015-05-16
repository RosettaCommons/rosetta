// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/AntibodyEnum.hh
/// @brief Enumerators + TypeDefs for Antibody namespace
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Jianqing Xu (xubest@gmail.com)

#ifndef  INCLUDED_protocols_antibody_AntibodyEnum_hh
#define INCLUDED_protocols_antibody_AntibodyEnum_hh

#include <utility/vector1.hh>
#include <core/types.hh>
#include <iostream>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A collection of Enumerators + TypeDefs to the Antibody Namespace.
// These keep AntibodyInfo light, allow the use of vectors instead of maps, and keep us from passing around strings.
// All Enums should have a total for resizing vectors and creating iterations.
// NOTE: Please see AntibodyEnumManager for conversions if you add anything here.
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace antibody {


enum CDRNameEnum {
	h1 = 1,
	h2,
	h3,
	l1,
	l2,
	l3,
	
	/// Proto, not quite CDR loops for design and analysis ///
	/// Also known as the DE loop ///
	
	proto_h4,
	proto_l4,
	
	///////// Convenience //////////////////
	start_cdr_loop = h1,
	H_chain_last_loop = h3,
	L_chain_last_loop = l3,
	camelid_last_loop = h3,
	num_cdr_loops = l3,
	CDRNameEnum_start = h1,
	CDRNameEnum_total = l3,
	
	//Proto CDR Loops
	CDRNameEnum_proto_start = proto_h4,
	CDRNameEnum_proto_total = proto_l4,
	l4 = proto_l4,
	h4 = proto_h4
	
};

enum AntibodyNumberingSchemeEnum {
	Chothia_Scheme = 1,
	Kabat_Scheme,
	Enhanced_Chothia_Scheme,
	IMGT_Scheme,
	AHO_Scheme,

	AntibodyNumberingSchemeEnum_start = Chothia_Scheme,
	AntibodyNumberingSchemeEnum_total = AHO_Scheme,
	NONE
};

enum CDRDefinitionEnum {
	Chothia = 1,
	Aroop,
	Kabat,
	Martin,
	North,
	
	CDRDefinitionEnum_start = Chothia,
	CDRDefinitionEnum_total = North
	
};

//Light chain types
enum LightChainTypeEnum {
	lambda =1,
	kappa,
	lambda6,
	unknown,
	LightChainTypeEnum_start = lambda,
	LightChainTypeEnum_total = unknown
};

//Overall region of the antibody.
// Note: 
//  _region designation helps for unique variable naming (aka cdr variable is passed around all over)
enum AntibodyRegionEnum {
	antigen_region = 1,
	cdr_region,
	framework_region,
	AntibodyRegionEnum_start = antigen_region,
	AntibodyRegionEnum_total = framework_region
};

///Main enumerator for AntibodyNumbering.
enum CDRLandmarkEnum {
	cdr_start = 1,
	cdr_end,
	CDRLandmarkEnum_total = cdr_end
};

enum H3BaseTypeEnum {
	Kinked = 1,
	Extended,
	Neutral,
	Unknown,
	H3BaseTypeEnum_start = Kinked,
	H3BaseTypeEnum_total = Unknown
};

///These are used to determine the VL_VH packing angle.
enum PackingAngleEnum {
	VL_sheet_1 = 1,
	VL_sheet_2,
	VH_sheet_1,
	VH_sheet_2,

	PackingAngleEnum_start = VL_sheet_1,
	PackingAngleEnum_total = VH_sheet_2
};

	
}//antibody
}//protocols
#endif	//#ifndef INCLUDED_protocols/antibody/AntibodyEnum_HH

