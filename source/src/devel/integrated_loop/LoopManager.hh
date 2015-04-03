// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///	  Contains currently: LoopModeler
///
///
/// @author Vatsan Raman

#ifndef INCLUDED_devel_integrated_loop_LoopManager_hh
#define INCLUDED_devel_integrated_loop_LoopManager_hh

//devel headers
//v#include <devel/IntegratedLoop/Loops.fwd.hh>

//core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>

// protocols header

//C++ headers

//Auto Headers
#include <protocols/loops/Loops.hh>


namespace protocols {
namespace moves {

//////////////////////////////////////////////////////////
/// @brief manages loop business
///
/////////////////////////////////////////////////////////

class LoopManager {

public:
  LoopManager(){}

  LoopManager(
    core::pose::Pose & pose,
    protocols::loops::Loops LoopList
  ) :
    pose_( pose ),
    LoopList_( LoopList )
  {
    ReorderLoops();
  }

  //  Size NumSingleLoops();
  protocols::loops::Loops LoopsToPerturb();
  void ReorderLoops();

private:
  core::pose::Pose pose_;
  protocols::loops::Loops LoopList_;

  bool IsFirstLoop( protocols::loops::Loop const & ThisLoop );
  bool IsLastLoop( protocols::loops::Loop const & ThisLoop );

  protocols::loops::Loop PreviousLoop( protocols::loops::Loop const & ThisLoop );
  protocols::loops::Loop NextLoop( protocols::loops::Loop const & ThisLoop );

  bool IsNtermLoop( protocols::loops::Loop const & ThisLoop );
  bool IsCtermLoop( protocols::loops::Loop const & ThisLoop );

  protocols::loops::Loop VaryStems( protocols::loops::Loop const & ThisLoop );
  protocols::loops::Loop VaryCutpoint( protocols::loops::Loop const & ThisLoop );
};

} //moves
}//protocols


#endif
