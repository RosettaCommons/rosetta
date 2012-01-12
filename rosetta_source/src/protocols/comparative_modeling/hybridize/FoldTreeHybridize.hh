// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief Align a random jump to template
/// @detailed
/// @author Yifan Song

#ifndef INCLUDED_protocols_comparative_modeling_hybridize_FoldTreeHybridize_hh
#define INCLUDED_protocols_comparative_modeling_hybridize_FoldTreeHybridize_hh

#include <protocols/comparative_modeling/hybridize/InsertChunkMover.hh>
#include <protocols/comparative_modeling/hybridize/FoldTreeHybridize.fwd.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/util/kinematics_util.hh>
#include <core/fragment/FragSet.fwd.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/scoring/ScoreFunction.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/rigid.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace comparative_modeling {
namespace hybridize {

using namespace core;
using namespace protocols::moves;
using namespace protocols::loops;
	
class FoldTreeHybridize: public protocols::moves::Mover
{
	
public:

	FoldTreeHybridize(core::Size const initial_template_index,
                      utility::vector1 < core::pose::PoseCOP > const & template_poses,
                      utility::vector1 < protocols::loops::Loops > const & template_chunks,
                      utility::vector1 < protocols::loops::Loops > const & template_contigs,
                      utility::vector1 < core::fragment::FragSetOP > & frag_libs);

	void revert_loops_to_original(core::pose::Pose & pose, Loops loops);

	void set_loops_to_virt_ala(core::pose::Pose & pose, Loops loops);

	Real
	gap_distance(Size Seq_gap);

	void add_gap_constraints_to_pose(core::pose::Pose & pose, Loops const & chunks, int gap_edge_shift=0, Real stdev=0.1);

    void backup_original_foldtree(core::pose::Pose const & pose);
    void restore_original_foldtree(core::pose::Pose & pose);

	void
	setup_foldtree(core::pose::Pose & pose);

    core::Size choose_anchor_position(const protocols::loops::Loop & chunk) const;

	numeric::xyzVector<Real> center_of_mass(core::pose::Pose const & pose);

	void translate_virt_to_CoM(core::pose::Pose & pose);

	utility::vector1< core::Real > get_residue_weights_from_loops(core::pose::Pose & pose);
    
    protocols::loops::Loops renumber_template_chunks(
                                                     protocols::loops::Loops & template_chunk,
                                                     core::pose::PoseCOP template_pose
                                                     );

	Loops loops();
    inline void set_scorefunction(core::scoring::ScoreFunctionOP const scorefxn) {
        scorefxn_ = scorefxn;
    }
    
	void apply(core::pose::Pose & pose);

	std::string	get_name() const;
	
private:
    core::Size initial_template_index_;
    core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1 < core::pose::PoseCOP > template_poses_;
	utility::vector1 < protocols::loops::Loops > template_chunks_;
	utility::vector1 < protocols::loops::Loops > template_contigs_;
	utility::vector1 < core::fragment::FragSetOP > frag_libs_;

	Loops ss_chunks_pose_;
	//Loops loops_pose_;
	
    // backup original info
    core::kinematics::FoldTree orig_ft_;
    Size orig_n_residue_;
}; //class FoldTreeHybridize
	
} // hybridize 
} // comparative_modeling 
} // protocols

#endif
