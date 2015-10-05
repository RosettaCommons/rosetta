// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/legacy_types.cc
/// @brief  types that were in use for patching by mm_params in 2015.
/// @author Rhiju Das

#include <core/chemical/legacy_types.hh>

namespace core {
namespace chemical {

  utility::vector1< VariantType > const &
  variant_types_list_LEGACY() {
    static  utility::vector1< VariantType > variant_types_list_LEGACY_;
    static bool init( false );
    if ( !init ) {
      variant_types_list_LEGACY_.push_back( UPPER_TERMINUS_VARIANT );
      variant_types_list_LEGACY_.push_back( LOWER_TERMINUS_VARIANT );
      variant_types_list_LEGACY_.push_back( NTERM_CONNECT );
      variant_types_list_LEGACY_.push_back( CTERM_CONNECT );
      variant_types_list_LEGACY_.push_back( CUTPOINT_LOWER );
      variant_types_list_LEGACY_.push_back( CUTPOINT_UPPER );
      variant_types_list_LEGACY_.push_back( METHYLATION );
      variant_types_list_LEGACY_.push_back( METHYLATED_NTERM_VARIANT );
      variant_types_list_LEGACY_.push_back( PHOSPHORYLATION );
      variant_types_list_LEGACY_.push_back( ACETYLATION );
      variant_types_list_LEGACY_.push_back( SULFATION );
      variant_types_list_LEGACY_.push_back( CARBOXYLATION );
      variant_types_list_LEGACY_.push_back( DIMETHYLATION );
      variant_types_list_LEGACY_.push_back( TRIMETHYLATION );
      variant_types_list_LEGACY_.push_back( DIIODINATION );
      variant_types_list_LEGACY_.push_back( ACETYLATED_NTERMINUS_VARIANT );
      variant_types_list_LEGACY_.push_back( METHYLATED_CTERMINUS_VARIANT );
      variant_types_list_LEGACY_.push_back( VIRTUAL_DNA_PHOSPHATE );
      variant_types_list_LEGACY_.push_back( VIRTUAL_PHOSPHATE );
      variant_types_list_LEGACY_.push_back( VIRTUAL_BACKBONE_EXCEPT_C1PRIME );
      variant_types_list_LEGACY_.push_back( BULGE );
      variant_types_list_LEGACY_.push_back( VIRTUAL_O2PRIME_HYDROGEN );
      variant_types_list_LEGACY_.push_back( THREE_PRIME_END_OH );
      variant_types_list_LEGACY_.push_back( FIVE_PRIME_END_OH );
      variant_types_list_LEGACY_.push_back( FIVE_PRIME_END_PHOSPHATE );
      variant_types_list_LEGACY_.push_back( OOP_PRE );
      variant_types_list_LEGACY_.push_back( OOP_POST );
      variant_types_list_LEGACY_.push_back( HBS_PRE );
      variant_types_list_LEGACY_.push_back( HBS_POST );
      variant_types_list_LEGACY_.push_back( A3B_HBS_PRE );
      variant_types_list_LEGACY_.push_back( A3B_HBS_POST );
      variant_types_list_LEGACY_.push_back( TRIAZOLAMERN );
      variant_types_list_LEGACY_.push_back( TRIAZOLAMERC );
      variant_types_list_LEGACY_.push_back( DISULFIDE );
      variant_types_list_LEGACY_.push_back( HYDROXYLATION1 );
      variant_types_list_LEGACY_.push_back( HYDROXYLATION2 );
      variant_types_list_LEGACY_.push_back( REPLS_BB );
      variant_types_list_LEGACY_.push_back( VIRTUAL_RESIDUE_VARIANT );
      variant_types_list_LEGACY_.push_back( VIRTUAL_NTERM );
      variant_types_list_LEGACY_.push_back( VIRTUAL_BB );
      variant_types_list_LEGACY_.push_back( SHOVE_BB );
      init = true;
    }
    return variant_types_list_LEGACY_;
  }

} // namespace chemical
} // namespace core
