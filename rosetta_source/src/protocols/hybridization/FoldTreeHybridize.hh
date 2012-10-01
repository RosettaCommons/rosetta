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

#ifndef INCLUDED_protocols_hybridization_FoldTreeHybridize_hh
#define INCLUDED_protocols_hybridization_FoldTreeHybridize_hh

#include <protocols/hybridization/InsertChunkMover.hh>
#include <protocols/hybridization/FoldTreeHybridize.fwd.hh>
#include <protocols/hybridization/HybridizeFoldtreeDynamic.hh>

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
//namespace comparative_modeling {
namespace hybridization {

using namespace core;
using namespace protocols::moves;
using namespace protocols::loops;

class FoldTreeHybridize: public protocols::moves::Mover {

public:
	FoldTreeHybridize();

	FoldTreeHybridize(
		core::Size const initial_template_index,
		utility::vector1 < core::pose::PoseCOP > const & template_poses,
		utility::vector1 < core::Real > const & template_weights,
		utility::vector1 < protocols::loops::Loops > const & template_chunks,
		utility::vector1 < protocols::loops::Loops > const & template_contigs,
		core::fragment::FragSetOP fragments3_in,
		core::fragment::FragSetOP fragments9_in );

	// initialize options to defaults
	void init();

	void revert_loops_to_original(core::pose::Pose & pose, Loops loops);

	void set_loops_to_virt_ala(core::pose::Pose & pose, Loops loops);

	Real gap_distance(Size Seq_gap);

	void add_gap_constraints_to_pose(core::pose::Pose & pose, Loops const & chunks, int gap_edge_shift=0, Real stdev=0.1);

	void backup_original_foldtree(core::pose::Pose const & pose);
	void restore_original_foldtree(core::pose::Pose & pose);

	void setup_foldtree(core::pose::Pose & pose);

	numeric::xyzVector<Real> center_of_mass(core::pose::Pose const & pose);

	void translate_virt_to_CoM(core::pose::Pose & pose);

	utility::vector1< core::Real > get_residue_weights_from_loops(core::pose::Pose & pose);

	protocols::loops::Loops renumber_template_chunks(
	          protocols::loops::Loops & template_chunk,
	          core::pose::PoseCOP template_pose);

	void set_constraint_file(std::string cst_file_in) { cst_file_=cst_file_in; }
	void set_increase_cycles(core::Real increase_cycles_in) { increase_cycles_=increase_cycles_in; }
	void set_add_non_init_chunks(bool add_non_init_chunks_in) { add_non_init_chunks_=add_non_init_chunks_in; }
	void set_domain_assembly(bool domain_assembly_in) { domain_assembly_=domain_assembly_in; }
	void set_add_hetatm(bool add_hetatm_in, core::Real hetatm_cst_weight_in) { add_hetatm_=add_hetatm_in; hetatm_cst_weight_=hetatm_cst_weight_in; }
	void set_frag_insertion_weight(core::Real frag_insertion_weight_in) { frag_insertion_weight_=frag_insertion_weight_in; }
	void set_frag_weight_aligned(core::Real frag_weight_aligned_in) { frag_weight_aligned_=frag_weight_aligned_in; }
	void set_max_registry_shift(core::Size max_registry_shift_in) { max_registry_shift_=max_registry_shift_in; }
	
	inline void set_scorefunction(core::scoring::ScoreFunctionOP const scorefxn) { scorefxn_ = scorefxn; }

	void setup_scorefunctions( 
		core::scoring::ScoreFunctionOP score0,
		core::scoring::ScoreFunctionOP score1,
		core::scoring::ScoreFunctionOP score2,
		core::scoring::ScoreFunctionOP score5,
		core::scoring::ScoreFunctionOP score3);

	void apply(core::pose::Pose & pose);

	std::string	get_name() const;

private:
	core::Real increase_cycles_;
	core::Real frag_insertion_weight_; // fragment insertion weight, vs. chunk insertion
	bool add_non_init_chunks_;
	bool domain_assembly_;
	bool add_hetatm_;
	core::Real hetatm_cst_weight_;
	core::Real frag_weight_aligned_; // fragment insertion to the aligned region, vs. unaligned region
	core::Size max_registry_shift_;
    std::string cst_file_;

	core::Size initial_template_index_;
	core::scoring::ScoreFunctionOP scorefxn_;
	utility::vector1 < core::Real > template_wts_;
	utility::vector1 < core::pose::PoseCOP > template_poses_;
	utility::vector1 < protocols::loops::Loops > template_chunks_;
	utility::vector1 < protocols::loops::Loops > template_contigs_;
	utility::vector1 < core::fragment::FragSetOP > frag_libs_;
	utility::vector1 < core::fragment::FragSetOP > frag_libs9_;
	utility::vector1 < core::fragment::FragSetOP > frag_libs3_;

	Loops ss_chunks_pose_;

	// backup original info
	HybridizeFoldtreeDynamic foldtree_mover_;
	//core::kinematics::FoldTree orig_ft_;
	//Size orig_n_residue_;
}; //class FoldTreeHybridize

} // hybridize
//} // comparative_modeling
} // protocols

#endif
