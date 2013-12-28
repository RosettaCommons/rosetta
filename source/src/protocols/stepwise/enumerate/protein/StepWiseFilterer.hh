// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_Filterer.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_StepWiseFilterer_hh
#define INCLUDED_protocols_stepwise_StepWiseFilterer_hh

#include <core/pose/Pose.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

#include <map>

#include <utility/vector1.hh>

//Auto Headers



namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {

	typedef std::map< std::string, core::pose::PoseOP > PoseList;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseFilterer: public utility::pointer::ReferenceCount {
  public:

    //constructor!
    StepWiseFilterer();

    //destructor -- necessary?
    ~StepWiseFilterer();

    /// @brief Filter a list of poses by score.
    void filter( PoseList & pose_list, PoseList & minimize_pose_list ) const;

		void set_final_number( core::Size const & setting ){ final_number_ = setting; }

  private:

		core::Size final_number_;

  };

} //protein
} //enumerate
} //stepwise
} //protocols

#endif
