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
#include <utility/vector1.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <core/pose/full_model_info/FullModelInfo.fwd.hh>

// C++
#include <string>
#include <map>

namespace core {
namespace pose {
namespace full_model_info {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA centroid information inside the pose.
//// Rhiju move this to its own namespace!
class FullModelInfo: public basic::datacache::CacheableData  {

public:

  FullModelInfo(); //empty...

  FullModelInfo( 	utility::vector1< Size > const & sub_to_full,
									utility::vector1< Size > const & moving_res_list,
									std::string  const full_sequence,
									utility::vector1< Size >  const & cutpoint_open_in_full_model); //proper constructor

	FullModelInfo( pose::Pose const & pose );

  FullModelInfo( FullModelInfo const & src );

  basic::datacache::CacheableDataOP
  clone() const
  {
    return new FullModelInfo( *this );
  }

	// properties of current, working pose
	utility::vector1< Size > sub_to_full() const { return sub_to_full_;} //mapping from working pose to full pose.

	utility::vector1< Size > moving_res_list() const{ return moving_res_list_;}

	// properties of full model.
	std::string full_sequence() const{ return full_sequence_;}

	// set properties of full model.
	void set_full_sequence( std::string const & setting ) { full_sequence_ = setting;}

	utility::vector1< Size > cutpoint_open_in_full_model() const { return cutpoint_open_in_full_model_;}

	// properties of current, working pose
	void set_sub_to_full( utility::vector1< Size > const & setting ){ sub_to_full_ = setting;} //mapping from working pose to full pose.

	void set_moving_res_list( utility::vector1< Size > const & setting ){ moving_res_list_ = setting;}

	utility::vector1< Size >
	get_res_num_from_pdb_info( pose::Pose const & pose ) const;

	utility::vector1< Size >
	get_cutpoint_open_from_pdb_info( pose::Pose const & pose ) const;

	std::string
	get_sequence_with_gaps_filled_with_n( pose::Pose const & pose ) const;

private:

	// properties of current, working pose
	utility::vector1< Size > sub_to_full_; //mapping from working pose to full pose.
	utility::vector1< Size > moving_res_list_;

	// properties of full model.
	std::string full_sequence_;
	utility::vector1< Size > cutpoint_open_in_full_model_;

};

FullModelInfo const &
const_full_model_info_from_pose( pose::Pose const & pose );

FullModelInfo &
nonconst_full_model_info_from_pose( core::pose::Pose & pose );


}
}
}
#endif
