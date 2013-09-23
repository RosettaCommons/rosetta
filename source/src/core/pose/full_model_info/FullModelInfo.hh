// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/FullModelInfo.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_pose_full_model_info_FullModelInfo_hh
#define INCLUDED_core_pose_full_model_info_FullModelInfo_hh

#include <core/types.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/pose/full_model_info/FullModelInfo.fwd.hh>

// C++
#include <string>
#include <map>

namespace core {
namespace pose {
namespace full_model_info {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of all information related to how a subpose 'fits in' to global modeling scheme.
class FullModelInfo: public basic::datacache::CacheableData  {

public:

  FullModelInfo(); //empty...

  FullModelInfo( std::string  const full_sequence );

	FullModelInfo( pose::Pose & pose,
								 std::string const & full_sequence,
								 utility::vector1< Size > const & cutpoint_open_in_full_model,
								 utility::vector1< Size > const & res_numbers_in_pose );

	FullModelInfo( pose::Pose & pose );

  FullModelInfo( FullModelInfo const & src );

  ~FullModelInfo();

	basic::datacache::CacheableDataOP
  clone() const
  {
    return new FullModelInfo( *this );
  }

	FullModelInfoOP
  clone_info() const
  {
    return new FullModelInfo( *this );
  }

	// properties of full model.
	std::string const & full_sequence() const{ return full_sequence_;}

	utility::vector1< Size > const & cutpoint_open_in_full_model() const { return cutpoint_open_in_full_model_;}

	utility::vector1< Size > const & fixed_domain_map() const { return fixed_domain_map_;}

	void clear_other_pose_list();

	utility::vector1< core::pose::PoseOP > const & other_pose_list() const { return other_pose_list_; }

	utility::vector1< Size > const & res_list() const { return res_list_; }

	utility::vector1< Size > full_to_sub( utility::vector1< Size > const & res_in_full_model_numbering ) const;

	void add_other_pose( core::pose::PoseOP & pose );

	// set properties of full model.
	void set_full_sequence( std::string const & setting ) { full_sequence_ = setting;}

	void set_cutpoint_open_in_full_model( utility::vector1< Size > const & setting ){ cutpoint_open_in_full_model_ = setting;}

	void set_fixed_domain_map( utility::vector1< Size > const & setting ){ fixed_domain_map_ = setting;}

	void set_other_pose_list( utility::vector1< pose::PoseOP > const & setting );

	void set_res_list( utility::vector1< Size > const & res_list ){ res_list_ = res_list; }

	void remove_other_pose_at_idx( Size const idx );

	Size find_index_in_other_pose_list( pose::Pose const & pose ) const;

	Size get_idx_for_other_pose_with_residue( Size const input_res ) const;

	utility::vector1< Size > chains_in_full_model() const;

	utility::vector1< Size > moving_res_in_full_model() const;

	Size size(){ return full_sequence_.size(); }

private:

	utility::vector1< Size >
	get_res_num_from_pdb_info( pose::Pose const & pose ) const;

	utility::vector1< Size >
	get_cutpoint_open_from_pdb_info( pose::Pose const & pose ) const;

	std::string
	get_sequence_with_gaps_filled_with_n( pose::Pose const & pose ) const;

private:

	// properties of full model.
	std::string full_sequence_;
	utility::vector1< Size > cutpoint_open_in_full_model_;
	utility::vector1< Size > fixed_domain_map_; //in case model involves several fixed regions.

	// residues that go with pose. In principle, this is redundant with PDBInfo,
	// but not in eventual case where user has favorite numbering/chain scheme.
	utility::vector1< Size > res_list_;

	// what's known about this pose and any neighbors in a "PoseTree"
	utility::vector1< core::pose::PoseOP > other_pose_list_;

	// Following would be great, but would cause a headache in cloning --
	//  the FullModelInfo doesn't know what pose it is inside.
	//  instead, we'll probably need to create a lightweight PoseTree object that
	//  infers daughters and parents from FullModelInfo when it is needed.
	//	PoseOP pose_tree_parent;


};

FullModelInfo const &
const_full_model_info( pose::Pose const & pose );

FullModelInfo &
nonconst_full_model_info( core::pose::Pose & pose );


}
}
}
#endif
