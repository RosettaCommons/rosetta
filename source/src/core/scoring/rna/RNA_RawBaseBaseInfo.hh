// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_RawBaseBasePotential.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_rna_RNA_RawBaseBaseInfo_hh
#define INCLUDED_core_scoring_rna_RNA_RawBaseBaseInfo_hh

#include <core/types.hh>

// Package headers
//#include <core/conformation/Residue.fwd.hh>

// Project headers
#include <basic/datacache/CacheableData.hh>

// Utility headers

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>

// C++

namespace core {
namespace scoring {
namespace rna {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA centroid information inside the pose.
//// Rhiju move this to its own namespace!
//// Also, should probably use EnergyGraph instead of FArrays -- much smaller memory footprint (!)
////
class RNA_RawBaseBaseInfo : public basic::datacache::CacheableData {

public:

	RNA_RawBaseBaseInfo(): calculated_(false) {};

  RNA_RawBaseBaseInfo( RNA_RawBaseBaseInfo const & src );

  basic::datacache::CacheableDataOP
  clone() const
  {
    return new RNA_RawBaseBaseInfo( *this );
  }

  Size
  size() const {
    return base_pair_array_.size1();
  }

  void
  resize( Size const & total_residue );

  void
  zero();

	void
	copy_values( scoring::rna::RNA_RawBaseBaseInfo const & src, Size const & i, Size const & j );

  // Undefinded, comented out to make python bindings complile
  //void
  //initialize( pose::Pose const & pose );

  bool
  calculated() const
  {
    return calculated_;
  }

  bool &
  calculated()
  {
    return calculated_;
  }

  void
  set_calculated( bool const & setting)
  {
    calculated_ = setting;
  }

  ObjexxFCL::FArray3D< Real > & base_pair_array()  { return  base_pair_array_; }
  ObjexxFCL::FArray3D< Real > & base_axis_array()  { return  base_axis_array_; }
  ObjexxFCL::FArray3D< Real > & base_stagger_array()  { return  base_stagger_array_; }
  ObjexxFCL::FArray2D< Real > & base_stack_array()  { return  base_stack_array_; }
  ObjexxFCL::FArray2D< Real > & base_stack_axis_array()  { return  base_stack_axis_array_; }
  ObjexxFCL::FArray2D< Real > & base_geometry_orientation_array()  { return  base_geometry_orientation_array_; }
  ObjexxFCL::FArray2D< Real > & base_geometry_height_array()  { return  base_geometry_height_array_; }

  ObjexxFCL::FArray3D< Real > const & base_pair_array() const { return  base_pair_array_; }
  ObjexxFCL::FArray3D< Real > const & base_axis_array() const { return  base_axis_array_; }
  ObjexxFCL::FArray3D< Real > const & base_stagger_array() const { return  base_stagger_array_; }
  ObjexxFCL::FArray2D< Real > const & base_stack_array() const { return  base_stack_array_; }
  ObjexxFCL::FArray2D< Real > const & base_stack_axis_array() const { return  base_stack_axis_array_; }
  ObjexxFCL::FArray2D< Real > const & base_geometry_orientation_array() const { return  base_geometry_orientation_array_; }
  ObjexxFCL::FArray2D< Real > const & base_geometry_height_array() const { return  base_geometry_height_array_; }


private:

  // For now, direct copy of what was in rosetta++.
  //The third dimension here refers to the edge of the base that is pairing.
  // Note that these are "scratch" arrays, not filtered to avoid,
  // e.g. one base edge forming multiple base pairs.
  ObjexxFCL::FArray3D< Real > base_pair_array_;
  ObjexxFCL::FArray3D< Real > base_axis_array_;
  ObjexxFCL::FArray3D< Real > base_stagger_array_;
  ObjexxFCL::FArray2D< Real > base_stack_array_;
  ObjexxFCL::FArray2D< Real > base_stack_axis_array_;

  ObjexxFCL::FArray2D< Real > base_geometry_orientation_array_;
  ObjexxFCL::FArray2D< Real > base_geometry_height_array_;

  bool calculated_;

};

}
}
}

#endif
