// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/RNA_FilteredBaseBasePotential.hh
/// @brief  Statistically derived RNA potential
/// @author Rhiju Das

#ifndef INCLUDED_core_scoring_rna_RNA_FilteredBaseBaseInfo_hh
#define INCLUDED_core_scoring_rna_RNA_FilteredBaseBaseInfo_hh

#include <core/types.hh>

// Package headers
//#include <core/conformation/Residue.fwd.hh>

// Project headers
#include <core/scoring/rna/RNA_RawBaseBaseInfo.fwd.hh>
#include <core/scoring/rna/RNA_BaseDoubletClasses.hh>
#include <core/scoring/rna/RNA_DataInfo.hh>

#include <basic/datacache/CacheableData.hh>

// Utility headers

// ObjexxFCL
#include <ObjexxFCL/FArray2D.hh>

#include <utility/vector1.hh>


// C++

namespace core {
namespace scoring {
namespace rna {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA centroid information inside the pose.
//// Rhiju move this to its own namespace!
class RNA_FilteredBaseBaseInfo : public basic::datacache::CacheableData {

public:

	RNA_FilteredBaseBaseInfo();

  RNA_FilteredBaseBaseInfo( RNA_FilteredBaseBaseInfo const & src );

  basic::datacache::CacheableDataOP
  clone() const
  {
    return new RNA_FilteredBaseBaseInfo( *this );
  }

  Size
  size() const {
    return filtered_base_pair_array_.size1();
  }

  void
  resize( Size const & total_residue );

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
  set_calculated( bool const & setting )
  {
    calculated_ = setting;
  }

	void
	carry_out_filtering( RNA_RawBaseBaseInfo const & raw_base_base_info );

  ObjexxFCL::FArray2D < Real > & filtered_base_pair_array()  { return  filtered_base_pair_array_; }
  ObjexxFCL::FArray2D < Real > & filtered_base_axis_array()  { return  filtered_base_axis_array_; }
  ObjexxFCL::FArray2D < Real > & filtered_base_stagger_array()  { return  filtered_base_stagger_array_; }
  ObjexxFCL::FArray2D < Real > & filtered_base_stack_array()  { return  filtered_base_stack_array_; }
  ObjexxFCL::FArray2D < Real > & filtered_base_stack_axis_array()  { return  filtered_base_stack_axis_array_; }

  ObjexxFCL::FArray2D < Real > const & filtered_base_pair_array() const { return  filtered_base_pair_array_; }
  ObjexxFCL::FArray2D < Real > const & filtered_base_axis_array() const { return  filtered_base_axis_array_; }
  ObjexxFCL::FArray2D < Real > const & filtered_base_stagger_array() const { return  filtered_base_stagger_array_; }
  ObjexxFCL::FArray2D < Real > const & filtered_base_stack_array() const { return  filtered_base_stack_array_; }
  ObjexxFCL::FArray2D < Real > const & filtered_base_stack_axis_array() const { return  filtered_base_stack_axis_array_; }

	Real const & get_total_base_pair_score() const { return total_base_pair_score_; }
	Real const & get_total_base_axis_score() const { return total_base_axis_score_; }
	Real const & get_total_base_stagger_score() const { return total_base_stagger_score_; }
	Real const & get_total_base_stack_score() const { return total_base_stack_score_; }
	Real const & get_total_base_stack_axis_score() const { return total_base_stack_axis_score_; }

	bool const & scale_axis_stagger() const { return scale_axis_stagger_; }
	Real const & basepair_axis_stagger_scaling() const { return basepair_axis_stagger_scaling_; }
	Real const & basestack_axis_scaling() const { return basestack_axis_scaling_; }

	Energy_base_pair_list const scored_base_pair_list() const{ return scored_base_pair_list_; }
	Energy_base_stack_list const scored_base_stack_list() const{ return scored_base_stack_list_; }

	Real get_data_score( RNA_DataInfo const & rna_data_info ) const;

private:

	void
	figure_out_rna_base_pairs_to_score( RNA_RawBaseBaseInfo const & raw_base_base_info );

	void
	figure_out_rna_base_stacks_to_score( RNA_RawBaseBaseInfo const & raw_base_base_info );

	// data

  // For now, direct copy of what was in rosetta++.
  //The third dimension here refers to the edge of the filtered_base that is pairing.
  // Note that these are "scratch" arrays, not filtered to avoid,
  // e.g. one filtered_base edge forming multiple filtered_base pairs.
  ObjexxFCL::FArray2D < Real > filtered_base_pair_array_;
  ObjexxFCL::FArray2D < Real > filtered_base_axis_array_;
  ObjexxFCL::FArray2D < Real > filtered_base_stagger_array_;
  ObjexxFCL::FArray2D < Real > filtered_base_stack_array_;
  ObjexxFCL::FArray2D < Real > filtered_base_stack_axis_array_;

	Energy_base_pair_list  scored_base_pair_list_;
	Energy_base_stack_list scored_base_stack_list_;

	Real total_base_pair_score_;
	Real total_base_axis_score_;
	Real total_base_stagger_score_;
	Real total_base_stack_score_;
	Real total_base_stack_axis_score_;

	bool scale_axis_stagger_;
	Real basepair_axis_stagger_scaling_;
	Real basestack_axis_scaling_;
	bool include_neighbor_base_stacks_;

  bool calculated_;

	bool rna_verbose_;
};

} //rna
} //scoring
} //core

#endif
