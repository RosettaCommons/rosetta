// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin VariantType.cc
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

/// @brief for use during loop modeling, at positions before a cutpoint
VariantType const CUTPOINT_LOWER( "CUTPOINT_LOWER" );
/// @brief for use during loop modeling, at positions after a cutpoint
VariantType const CUTPOINT_UPPER( "CUTPOINT_UPPER" );
///
VariantType const DISULFIDE( "DISULFIDE" );
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
VariantType const VIRTUAL_O2STAR_HYDROGEN( "VIRTUAL_O2STAR_HYDROGEN" );

// The following are added by Andy M. Chen in July 2009 to be used for PTM patches/variants
VariantType const PHOSPHORYLATION( "PHOSPHORYLATION" );
VariantType const ACETYLATION( "ACETYLATION" );
VariantType const SULFATION( "SULFATION" );
VariantType const CARBOXYLATION( "CARBOXYLATION" );
VariantType const HYDROXYLATION( "HYDROXYLATION" );
VariantType const DIMETHYLATION( "DIMETHYLATION" );
VariantType const TRIMETHYLATION( "TRIMETHYLATION" );
VariantType const DIIODINATION( "DIIODINATION" );

/// @brief Actyleated N-terminus cap, written for creating amino acid dipeptides for NCAA rotamer libraries
VariantType const ACTYLATED_NTERMINUS( "ACTYLATED_NTERMINUS" );
/// @brief Methylated C-terminus cap, written for creating amino acid dipeptides for NCAA rotamer libraries
VariantType const METHYLATED_CTERMINUS( "METHYLATED_CTERMINUS" );
VariantType const SC_ORBITALS("SC_ORBITALS");

///@brief Cap extensions at termini to include peptide bonds, written for stepwise assembly (SWA) code.
VariantType const N_ACETYLATION( "N_ACETYLATION" );
VariantType const C_METHYLAMIDATION( "C_METHYLAMIDATION" );

///@ brief only the repulsive energy will be considered during structure calculations
VariantType const REPLONLY("REPLONLY");

///@ brief oop_pre patch, used for oligooxopiperazines (OOPs)
VariantType const OOP_PRE("OOP_PRE");
///@ brief oop_post patch, used for oligooxopiperazines (OOPs)
VariantType const OOP_POST("OOP_POST");

///@brief This is used for chemically conjugable residues (LYX, CYX)
///used for sidechain conjugation (like ubiquitination)
VariantType const SIDECHAIN_CONJUGATION("SIDECHAIN_CONJUGATION");

} // chemical
} // core
