// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/InvKinLigLoopDesign/Cloner.hh
///
/// @brief
/// @author

#ifndef DEVEL_INVKINLIGLOOPDESIGN_CLONER_HH
#define DEVEL_INVKINLIGLOOPDESIGN_CLONER_HH

#include <vector>
#include <map>


#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <devel/inv_kin_lig_loop_design/ResID.hh>
#include <devel/inv_kin_lig_loop_design/Loop.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace devel {

namespace inv_kin_lig_loop_design {

  using namespace std;
  using core::conformation::Residue;
  using platform::Size;

  // ================================================
  // ==================== Cloner ====================
  // ================================================

  typedef map< Residue*, Residue* > clones_type;
  typedef vector<Segment> segments_type;
  typedef vector<Loop> loops_type;

  struct Cloner {

    Cloner(TagCOP const& tag0, core::pose::PoseOP pose);

    core::pose::PoseOP clone();
    core::kinematics::FoldTree getFoldTree();
    void setInitialConfig();
    loops_type getLoops();

  private:

    TagCOP tag0;

    core::pose::PoseOP pose0;
    core::pose::PoseOP pose1;

    clones_type clones;
    resids_type resids;
    segments_type indels,segments,loops;

  }; // class Cloner

} // namespace LoopDesign

} // namespace devel

#endif // DEVEL_LOOPDESIGN_CLONER_HH
