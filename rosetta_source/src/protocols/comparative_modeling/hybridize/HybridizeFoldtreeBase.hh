// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/comparative_modeling/hybridize/HybridizeFoldtreeBase.hh
/// @author Yifan Song

#ifndef PROTOCOLS_COMPARATIVE_MODELING_HYBRIDIZE_HYBRIDIZEFOLDTREEBASE_HH_
#define PROTOCOLS_COMPARATIVE_MODELING_HYBRIDIZE_HYBRIDIZEFOLDTREEBASE_HH_

// Unit headers
#include <protocols/comparative_modeling/hybridize/HybridizeFoldtreeBase.fwd.hh>

// C/C++ headers
#include <string>

// Project headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <protocols/loops/Loops.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {
namespace hybridize {

class HybridizeFoldtreeBase {
 public:
  HybridizeFoldtreeBase();

    // initialize pose fold tree with chunks, add virtual residue
    void initialize(core::pose::Pose & pose);
    // update pose fold tree with chunks, no adding virtual residue
    void update(core::pose::Pose & pose);
    
    void save_foldtree(core::pose::Pose const & pose);

    /// @brief Restore to the saved foldtree, remove the virtual residue added to the end of the pose
    void restore_foldtree(core::pose::Pose & pose);

    void set_chunks(const protocols::loops::Loops & chunks,
                    const utility::vector1 < core::Size > & anchor_positions);
    void update(core::pose::Pose const & pose);

protected:
  /// @brief Stochastically selects an anchor position
  //core::Size choose_anchor_position(const protocols::loops::Loop& chunk) const;

private:
    /// @brief Index of the virtual residue we added to the pose in set_up()
    int virtual_res_;
    core::Size num_nonvirt_residues_;
    
    protocols::loops::Loops chunks_last_; 
    utility::vector1 < core::Size > anchor_positions_last_;
    
    protocols::loops::Loops chunks_; 
    utility::vector1 < core::Size > anchor_positions_;
    
    // backup original info
    core::kinematics::FoldTree saved_ft_;
    core::Size saved_n_residue_;

};

}  //  namespace comparative_modeling
}  //  namespace hybridize
}  //  namespace protocols

#endif  // PROTOCOLS_COMPARATIVE_MODELING_HYBRIDIZE_HYBRIDIZEFOLDTREEBASE_HH_
