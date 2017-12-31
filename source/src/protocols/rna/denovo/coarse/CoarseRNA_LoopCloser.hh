// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file loopRNA_minimizer.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_coarse_rna_CoarseRNA_LoopCloser_HH
#define INCLUDED_protocols_coarse_rna_CoarseRNA_LoopCloser_HH

#include <protocols/moves/Mover.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <utility/fixedsizearray1.fwd.hh>
#include <core/types.hh>
#include <ObjexxFCL/FArray1D.hh>

//// C++ headers
#include <string>


namespace protocols {
namespace rna {
namespace denovo {
namespace coarse {

/// @brief The RNA de novo structure modeling protocol
class CoarseRNA_LoopCloser: public protocols::moves::Mover {
public:
	/// @brief Construct the protocol object
	CoarseRNA_LoopCloser();

	/// @brief Clone this object
	protocols::moves::MoverOP clone() const override {
		return protocols::moves::MoverOP( new CoarseRNA_LoopCloser(*this) );
	}

	/// @brief Apply the loop-rebuild protocol to the input pose
	using protocols::moves::Mover::apply;

	void apply( core::pose::Pose & pose ) override;

	/// @brief Apply the loop-rebuild protocol to the input pose
	bool
	apply( core::pose::Pose & pose, Size const & seqpos_moved );

	bool
	apply_after_jump_change( core::pose::Pose & pose, Size const & jumpno );

	std::string get_name() const override;

	void
	choose_best_solution_based_on_score_function( core::scoring::ScoreFunctionOP scorefxn );

	void
	choose_least_perturb_solution();

	// Undefined, commenting out to fix PyRosetta build  void choose_random_solution();

	void
	get_all_solutions( core::pose::Pose & pose,
		utility::vector1< core::pose::PoseOP > & pose_list );

	void
	set_atom_level_domain_map( core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map );

	core::pose::toolbox::AtomLevelDomainMapOP  atom_level_domain_map();

	Size nsol(){ return nsol_; }

	utility::vector1< core::Size > const & cutpos_list(){ return cutpos_list_; }

	ObjexxFCL::FArray1D< bool > const & partition_definition(){ return partition_definition_; }

private:

	bool
	close_at_all_cutpoints( core::pose::Pose & pose );

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	figure_out_which_cutpoints_were_affected( core::pose::Pose const & pose );

	////////////////////////////////////////////////////////////////////////////////////////
	bool
	figure_out_pivots( core::pose::Pose const & pose );

	////////////////////////////////////////////////////////////////////////////////////////
	void
	remove_res( utility::vector1< core::Size > & res_vector,
		Size const & res );

	void
	backtrack( core::kinematics::tree::Atom const * current_atom,
		utility::vector1< core::Size > & upstream_res,
		utility::vector1< bool > & is_upstream_res,
		core::pose::Pose const & pose );

	void
	figure_out_forward_backward_res_by_backtracking( core::pose::Pose const & pose );

	void
	output_forward_backward_res();

	void
	filter_path( utility::vector1< core::Size > & upstream_res,
		utility::vector1< bool > const & is_filter_res,
		core::pose::Pose const & pose );

	void
	figure_out_pivot_res_and_scratch_res();

	void
	close_the_loop( core::pose::Pose & pose );

	void
	figure_out_dof_ids_and_offsets( core::pose::Pose const & pose,
		utility::vector1<core::Real> const & dt_ang );


	void
	figure_out_offset(
		core::pose::Pose const & pose,
		core::Size const & pivot,
		core::id::DOF_ID const & dof_id,
		core::Real const & original_torsion_value,
		utility::vector1< core::Real > & offset_save );

	void
	apply_solutions( core::pose::Pose & pose );

	void
	fill_solution( core::pose::Pose & pose,
		Size const n ) const;

	void
	output_chainTORS( utility::vector1< core::Real > const & dt_ang,
		utility::vector1< core::Real > const & db_ang,
		utility::vector1< core::Real > const & db_len ) const;

	void
	fill_chainTORS(
		core::pose::Pose const & pose,
		utility::vector1< core::id::NamedAtomID> const & atom_ids,
		utility::vector1<utility::fixedsizearray1<core::Real,3> > & atoms,
		utility::vector1<core::Real> & dt_ang,
		utility::vector1<core::Real> & db_ang,
		utility::vector1<core::Real> & db_len) const;

private:

	bool const a_little_verbose_;
	bool const verbose_;
	core::Size seqpos_moved_;
	core::Size cutpos_;
	utility::vector1< core::Size > cutpos_list_;
	int nsol_;

	core::pose::toolbox::AtomLevelDomainMapOP atom_level_domain_map_;
	utility::vector1< core::Size > backward_res_, forward_res_, pivot_res_, scratch_res_, pivot_to_scratch_res_;
	utility::vector1< bool > is_backward_res_, is_forward_res_, is_pivot_res_, is_scratch_res_;
	core::Size which_scratch_res_is_cut_;

	core::scoring::ScoreFunctionOP scorefxn_;

	utility::vector1< core::id::NamedAtomID > atom_ids_;
	utility::vector1< core::Real > offset_save1_, offset_save2_;
	utility::vector1< core::id::DOF_ID > dof_ids1_, dof_ids2_;

	utility::vector1<utility::vector1<core::Real> > t_ang_, b_ang_, b_len_;

	bool choose_least_perturb_solution_;
	bool choose_best_solution_;
	bool choose_random_solution_;

	ObjexxFCL::FArray1D< bool > partition_definition_;

}; // class CoarseRNA_LoopCloser


} //coarse
} //denovo
} //rna
} //protocols

#endif
