// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/VariantType.hh
/// @author Phil Bradley
// see VariantType.cc for explanation

#ifndef INCLUDED_core_chemical_VariantType_hh
#define INCLUDED_core_chemical_VariantType_hh

// Unit headers
#include <core/chemical/VariantType.fwd.hh>


namespace core {
namespace chemical {

extern VariantType const UPPER_TERMINUS_VARIANT;
extern VariantType const LOWER_TERMINUS_VARIANT;
extern VariantType const UPPERTERM_TRUNC_VARIANT;
extern VariantType const LOWERTERM_TRUNC_VARIANT;
extern VariantType const CUTPOINT_UPPER;
extern VariantType const CUTPOINT_LOWER;
extern VariantType const DISULFIDE;
extern VariantType const BRANCH_POINT_VARIANT;
extern VariantType const BRANCH_LOWER_TERMINUS_VARIANT;
extern VariantType const METHYLATION;
extern VariantType const ADDUCT_VARIANT;
extern VariantType const CENTROID_HA;
extern VariantType const PROTONATED;
extern VariantType const DEPROTONATED;
extern VariantType const SPECIAL_ROT;
extern VariantType const VIRTUAL_PHOSPHATE;
extern VariantType const VIRTUAL_RNA_RESIDUE;
extern VariantType const VIRTUAL_O2PRIME_HYDROGEN;

// The following are added by Andy M. Chen in July 2009 to be used for PTM patches/variants
extern VariantType const PHOSPHORYLATION;
extern VariantType const ACETYLATION;
extern VariantType const SULFATION;
extern VariantType const CARBOXYLATION;
extern VariantType const HYDROXYLATION;
extern VariantType const DIMETHYLATION;
extern VariantType const TRIMETHYLATION;
extern VariantType const DIIODINATION;

extern VariantType const ACETYLATED_NTERMINUS_VARIANT;
extern VariantType const METHYLATED_CTERMINUS_VARIANT;

// this has different geometry/atoms then ACETYLATED_NTERMINUS above:
extern VariantType const N_ACETYLATION;
// this is distinct from METHYLATED_CTERMINUS above:
extern VariantType const C_METHYLAMIDATION;

extern VariantType const SC_ORBITALS_VARIANT;

extern VariantType const REPLONLY;

extern VariantType const NTERM_CONNECT;
extern VariantType const CTERM_CONNECT;

extern VariantType const OOP_PRE;
extern VariantType const OOP_POST;

extern VariantType const HBS_PRE;
extern VariantType const HBS_POST;

// carbohydrate-specific variants
extern VariantType const C1_MODIFIED_SUGAR;
extern VariantType const C2_MODIFIED_SUGAR;
extern VariantType const C3_MODIFIED_SUGAR;
extern VariantType const C4_MODIFIED_SUGAR;
extern VariantType const C5_MODIFIED_SUGAR;
extern VariantType const C6_MODIFIED_SUGAR;
extern VariantType const C7_MODIFIED_SUGAR;
extern VariantType const C8_MODIFIED_SUGAR;
extern VariantType const C9_MODIFIED_SUGAR;

//This is used for chemically conjugable residues (LYX, CYX) used for sidechain conjugation (like ubiquitination)
extern VariantType const SIDECHAIN_CONJUGATION;

} // chemical
} // core

#endif
