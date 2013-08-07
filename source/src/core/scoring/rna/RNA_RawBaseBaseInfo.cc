// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_RawBaseBaseInfo.cc
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay
/// @author Rhiju Das


// Unit headers
#include <core/scoring/rna/RNA_RawBaseBaseInfo.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_RawBaseBaseInfo.fwd.hh>

// Package headers
#include <core/chemical/rna/RNA_Util.hh>

// Project headers
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/pose/Pose.hh>

// Utility headers
// AUTO-REMOVED #include <numeric/xyzMatrix.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>

#include <utility/vector1.hh>


// C++

using namespace core::chemical::rna;


///////////////////////////////////////////////////////
// Keep track of some base geometry that is
// useful for RNA scoring.
///////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace rna {


/// @details Copy constructors must copy all data, not just some...
RNA_RawBaseBaseInfo::RNA_RawBaseBaseInfo( RNA_RawBaseBaseInfo const & src ) :
	CacheableData()
{
  base_pair_array_ = src.base_pair_array_;
  base_axis_array_ = src.base_axis_array_;
  base_stagger_array_ = src.base_stagger_array_;
  base_stack_array_ = src.base_stack_array_;
  base_stack_axis_array_ = src.base_stack_axis_array_;
  base_geometry_orientation_array_ = src.base_geometry_orientation_array_;
  base_geometry_height_array_ = src.base_geometry_height_array_;
  calculated_ = src.calculated_;
}

void
RNA_RawBaseBaseInfo::resize( Size const & total_residue )
{

	if ( base_pair_array_.size1() == total_residue ) return;

  base_pair_array_.dimension( total_residue, total_residue, 3 );
  base_stagger_array_.dimension( total_residue, total_residue, 3 );
  base_axis_array_.dimension( total_residue, total_residue, 3 );
  base_stack_array_.dimension( total_residue, total_residue );
  base_stack_axis_array_.dimension( total_residue, total_residue );
  base_geometry_orientation_array_.dimension( total_residue, total_residue );
  base_geometry_height_array_.dimension( total_residue, total_residue );
}

	////////////////////////////////////////////////////////
void
RNA_RawBaseBaseInfo::zero()
{
  base_pair_array_ = 0.0;
  base_stagger_array_ = 0.0;
  base_axis_array_ = 0.0;
  base_stack_array_ = 0.0;
  base_stack_axis_array_ = 0.0;
  base_geometry_orientation_array_ = 0.0;
  base_geometry_height_array_ = 0.0;
}

	////////////////////////////////////////////////////////
void
RNA_RawBaseBaseInfo::copy_values( RNA_RawBaseBaseInfo const & src, Size const & i, Size const & j )
{
	for (Size k = 1; k <= NUM_EDGES; k++ ){
		base_pair_array_(i,j,k) = src.base_pair_array_(i,j,k);
		base_stagger_array_(i,j,k) = src.base_stagger_array_(i,j,k);
		base_axis_array_(i,j,k) = src.base_axis_array_(i,j,k);
	}
	base_stack_array_(i,j) = src.base_stack_array_(i,j);
  base_stack_axis_array_(i,j) = src.base_stack_axis_array_(i,j);
  base_geometry_orientation_array_(i,j) = src.base_geometry_orientation_array_(i,j);
  base_geometry_height_array_(i,j) = src.base_geometry_height_array_(i,j);
}

}
}
}
