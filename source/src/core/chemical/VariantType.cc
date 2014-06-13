// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @file core/chemical/VariantType.cc
///
/// @brief
/// VariantTypes are used by Patches.
///
/// @details
/// VariantTypes are utilized by Patches. All the type does is add a name
/// that can be used later on in different protocols. It also helps the patch
/// system keep track of what residues are patched with what type. The type has
/// no magical meaning. It is a name handler.
///
/// @author
/// Phil Bradley
/// Steven Combs - comments only
///
/// @last_modified October 22 2010
/////////////////////////////////////////////////////////////////////////
/// @file
/// @author Phil Bradley

// // Unit headers
#include <core/chemical/VariantType.hh>

// // Package headers

// Project headers

// Utility headers
//#include <utility/vector1.hh>
//#include <utility/pointer/owning_ptr.hh>
//#include <utility/pointer/ReferenceCount.hh>

// C++ headers

namespace core {
namespace chemical {

/// @brief C-terminus cap
VariantType const UPPER_TERMINUS( "UPPER_TERMINUS" );
/// @brief N-terminus cap
VariantType const LOWER_TERMINUS( "LOWER_TERMINUS" );

/// @brief C-terminus truncation
VariantType const UPPERTERM_TRUNC( "UPPERTERM_TRUNC" );
/// @brief N-terminus truncation
VariantType const LOWERTERM_TRUNC( "LOWERTERM_TRUNC" );

/// @brief for use during loop modeling, at positions before a cutpoint
VariantType const CUTPOINT_LOWER( "CUTPOINT_LOWER" );
/// @brief for use during loop modeling, at positions after a cutpoint
VariantType const CUTPOINT_UPPER( "CUTPOINT_UPPER" );
///
VariantType const DISULFIDE( "DISULFIDE" );
/// @brief Variant type used for branched polymers and glycosylations.
VariantType const BRANCH_POINT("BRANCH_POINT");
/// @brief Variant type used for branched polymers and glycosylations.
VariantType const BRANCH_LOWER_TERMINUS("BRANCH_LOWER_TERMINUS");
///
VariantType const ADDUCT( "ADDUCT" );

VariantType const METHYLATION( "METHYLATION" );

VariantType const CENTROID_HA( "CENTROID_WITH_HA" );

VariantType const PROTONATED( "PROTONATED");

VariantType const DEPROTONATED( "DEPROTONATED" );

/// @brief Generic variant type that allows for differential scoring of a set of residues/rotamers
VariantType const SPECIAL_ROT( "SPECIAL_ROT" );

VariantType const VIRTUAL_PHOSPHATE( "VIRTUAL_PHOSPHATE" );
VariantType const VIRTUAL_RNA_RESIDUE( "VIRTUAL_RNA_RESIDUE" );
VariantType const VIRTUAL_O2PRIME_HYDROGEN( "VIRTUAL_O2PRIME_HYDROGEN" );

// The following are added by Andy M. Chen in July 2009 to be used for PTM patches/variants
VariantType const PHOSPHORYLATION( "PHOSPHORYLATION" );
VariantType const ACETYLATION( "ACETYLATION" );
VariantType const SULFATION( "SULFATION" );
VariantType const CARBOXYLATION( "CARBOXYLATION" );
VariantType const HYDROXYLATION( "HYDROXYLATION" );
VariantType const DIMETHYLATION( "DIMETHYLATION" );
VariantType const TRIMETHYLATION( "TRIMETHYLATION" );
VariantType const DIIODINATION( "DIIODINATION" );

/// @brief Acetylated N-terminus cap, written for creating amino acid dipeptides for NCAA rotamer libraries
VariantType const ACETYLATED_NTERMINUS( "ACETYLATED_NTERMINUS" );
/// @brief Methylated C-terminus cap, written for creating amino acid dipeptides for NCAA rotamer libraries
VariantType const METHYLATED_CTERMINUS( "METHYLATED_CTERMINUS" );
VariantType const SC_ORBITALS("SC_ORBITALS");

///@brief Cap extensions at termini to include peptide bonds, written for stepwise assembly (SWA) code.
VariantType const N_ACETYLATION( "N_ACETYLATION" );
VariantType const C_METHYLAMIDATION( "C_METHYLAMIDATION" );

///@ brief only the repulsive energy will be considered during structure calculations
VariantType const REPLONLY("REPLONLY");

// @breif N-terminal connect and C-terminal connect
VariantType const NTERM_CONNECT( "NTERM_CONNECT" );
VariantType const CTERM_CONNECT( "CTERM_CONNECT" );

///@ brief oop_pre patch, used for oligooxopiperazines (OOPs)
VariantType const OOP_PRE("OOP_PRE");
///@ brief oop_post patch, used for oligooxopiperazines (OOPs)
VariantType const OOP_POST("OOP_POST");

///@ brief hbs_pre patch, used for hydrogen bond surrogates
VariantType const HBS_PRE("HBS_PRE");
///@ brief hbs_post patch, used for hydrogen bond surrogates
VariantType const HBS_POST("HBS_POST");


// carbohydrate-specific variants
///@ brief variant for any saccharide residue modified at the 1 position
VariantType const C1_MODIFIED_SUGAR("C1_MODIFIED_SUGAR");

///@ brief variant for any saccharide residue modified at the 2 position
VariantType const C2_MODIFIED_SUGAR("C2_MODIFIED_SUGAR");

///@ brief variant for any saccharide residue modified at the 3 position
VariantType const C3_MODIFIED_SUGAR("C3_MODIFIED_SUGAR");

///@ brief variant for any saccharide residue modified at the 4 position
VariantType const C4_MODIFIED_SUGAR("C4_MODIFIED_SUGAR");

///@ brief variant for any saccharide residue modified at the 5 position
VariantType const C5_MODIFIED_SUGAR("C5_MODIFIED_SUGAR");

///@ brief variant for any saccharide residue modified at the 6 position
VariantType const C6_MODIFIED_SUGAR("C6_MODIFIED_SUGAR");

///@ brief variant for any saccharide residue modified at the 7 position
VariantType const C7_MODIFIED_SUGAR("C7_MODIFIED_SUGAR");

///@ brief variant for any saccharide residue modified at the 8 position
VariantType const C8_MODIFIED_SUGAR("C8_MODIFIED_SUGAR");

///@ brief variant for any saccharide residue modified at the 9 position
VariantType const C9_MODIFIED_SUGAR("C9_MODIFIED_SUGAR");


///@brief This is used for chemically conjugable residues (LYX, CYX)
///used for sidechain conjugation (like ubiquitination)
VariantType const SIDECHAIN_CONJUGATION("SIDECHAIN_CONJUGATION");

} // chemical
} // core
