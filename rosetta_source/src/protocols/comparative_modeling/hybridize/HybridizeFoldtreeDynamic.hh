// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/comparative_modeling/hybridize/HybridizeFoldtreeDynamic.hh
/// @author Yifan Song

#ifndef PROTOCOLS_COMPARATIVE_MODELING_HYBRIDIZE_HYBRIDIZEFOLDTREEDYNAMIC_HH_
#define PROTOCOLS_COMPARATIVE_MODELING_HYBRIDIZE_HYBRIDIZEFOLDTREEDYNAMIC_HH_

// Unit headers
#include <protocols/comparative_modeling/hybridize/HybridizeFoldtreeDynamic.fwd.hh>
#include <protocols/comparative_modeling/hybridize/HybridizeFoldtreeBase.hh>

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

class HybridizeFoldtreeDynamic:  public HybridizeFoldtreeBase {
public:
    HybridizeFoldtreeDynamic(                                                   utility::vector1 < protocols::loops::Loops > template_contigs
                             );
    
    void set_core_chunks(const protocols::loops::Loops & chunks);

    protocols::loops::Loops make_complete_chunks( utility::vector1 < core::Size > cut_positions );

    void set_complete_chunks(const protocols::loops::Loops & chunks);

    void initialize(
                    core::pose::Pose & pose,
                    protocols::loops::Loops const & core_chunks
                    );
    
    utility::vector1 < core::Size > decide_cuts(core::Size n_residues);
    void choose_anchors();
    protocols::loops::Loops make_complete_chunks(utility::vector1 < core::Size > cut_positions,
                              core::Size n_residues);

protected:
    /// @brief Stochastically selects an anchor position
    core::Size choose_anchor_position(const protocols::loops::Loop& chunk) const;
    
private:
    utility::vector1 < protocols::loops::Loops > template_chunks_;
	utility::vector1 < protocols::loops::Loops > template_contigs_;

    utility::vector1 < core::Size > anchor_positions_;

    protocols::loops::Loops complete_chunks_last_; 
    protocols::loops::Loops complete_chunks_; 

    protocols::loops::Loops core_chunks_last_; 
    protocols::loops::Loops core_chunks_; 
    
};

}  //  namespace comparative_modeling
}  //  namespace hybridize
}  //  namespace protocols

#endif  // PROTOCOLS_COMPARATIVE_MODELING_HYBRIDIZE_HYBRIDIZEFOLDTREEDYNAMIC_HH_
