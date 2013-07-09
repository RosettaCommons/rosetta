// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SubToFullInfo.hh
/// @brief  Statistically derived rotamer pair potential class implementation
/// @author Phil Bradley
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_protocols_swa_monte_carlo_SubToFullInfo_hh
#define INCLUDED_protocols_swa_monte_carlo_SubToFullInfo_hh

#include <core/types.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <basic/datacache/CacheableData.hh>
#include <protocols/swa/monte_carlo/SubToFullInfo.fwd.hh>

// C++
#include <string>
#include <map>

namespace protocols {
namespace swa {
namespace monte_carlo {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Keep track of RNA centroid information inside the pose.
//// Rhiju move this to its own namespace!
class SubToFullInfo: public basic::datacache::CacheableData  {

public:

  SubToFullInfo(); //empty...

  SubToFullInfo( 	std::map< Size, Size >  sub_to_full,
									utility::vector1< Size >  moving_res_list,
									std::string  full_sequence,
									utility::vector1< Size >  cutpoints_in_full_pose); //proper constructor

  SubToFullInfo( SubToFullInfo const & src );

  basic::datacache::CacheableDataOP
  clone() const
  {
    return new SubToFullInfo( *this );
  }

	// properties of current, working pose
	std::map< Size, Size > sub_to_full(){ return sub_to_full_;} //mapping from working pose to full pose.

	utility::vector1< Size > moving_res_list() const{ return moving_res_list_;}

	// properties of full model.
	std::string full_sequence() const{ return full_sequence_;}

	utility::vector1< Size > cutpoints_in_full_pose() const{ return cutpoints_in_full_pose_;}

	// properties of current, working pose
	void set_sub_to_full( std::map< Size, Size > const & setting ){ sub_to_full_ = setting;} //mapping from working pose to full pose.

	void set_moving_res_list( utility::vector1< Size > const & setting ){ moving_res_list_ = setting;}

private:

	// properties of current, working pose
	std::map< Size, Size > sub_to_full_; //mapping from working pose to full pose.
	utility::vector1< Size > moving_res_list_;

	// properties of full model.
	std::string full_sequence_;
	utility::vector1< Size > cutpoints_in_full_pose_;

};

// Undefinded, commenting out to fix PyRosetta build  SubToFullInfo const & sub_to_full_info_from_pose( core::pose::Pose const & pose );

SubToFullInfo &
nonconst_sub_to_full_info_from_pose( core::pose::Pose & pose );

}
}
}
#endif
