// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/AntibodyEnum.hh
/// @brief Enumerators + TypeDefs for Antibody namespace
/// @author Jianqing Xu (xubest@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

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
using utility::vector1;
using core::Size;

enum CDRNameEnum {
    h1 = 1,
    h2,
    h3,
    l1,
    l2,
    l3,
    ///////// Convenience //////////////////
    start_cdr_loop = h1,
    H_chain_last_loop = h3,
    L_chain_last_loop = l3,
    camelid_last_loop = h3,
    num_cdr_loops = l3,
    CDRNameEnum_start = h1,
    CDRNameEnum_total = l3
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

