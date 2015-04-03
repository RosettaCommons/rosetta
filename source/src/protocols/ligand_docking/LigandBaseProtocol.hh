// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/LigandBaseProtocol.hh
///
/// @brief
/// @author Ian W. Davis


#ifndef INCLUDED_protocols_ligand_docking_LigandBaseProtocol_hh
#define INCLUDED_protocols_ligand_docking_LigandBaseProtocol_hh

#include <core/types.hh>
#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <protocols/moves/Mover.hh>


#include <set>

#include <core/conformation/Residue.fwd.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>


namespace protocols {
namespace ligand_docking {


/// @brief Convenience wrapper: selects the best ligand docking results
/// from a silent file and appends their tags to the supplied set.
void select_best_poses(
	core::import_pose::atom_tree_diffs::AtomTreeDiff const & atdiff,
	std::set< std::string > & tags_out
);

/// @brief Selects the best ligand docking results from a silent file
/// and appends their scores to the supplied list.
void select_best_poses(
	core::import_pose::atom_tree_diffs::AtomTreeDiff const & atdiff,
	core::import_pose::atom_tree_diffs::ScoresPairList & scores_out,
	core::Real to_keep = 0.05
);

/// @brief Trims scores_in based on ligand_is_touching (if present) and
/// then by total_score.
void select_best_poses(
	core::import_pose::atom_tree_diffs::ScoresPairList const & scores_in,
	core::import_pose::atom_tree_diffs::ScoresPairList & scores_out,
	core::Real to_keep = 0.05
);

/// @brief Without superimposing, automorphically computes the fraction of atoms
/// in these residues that are within the given cutoff(s) of each other.
void frac_atoms_within(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	utility::vector1< core::Real > const & cutoffs,
	utility::vector1< core::Real > & fractions_out
);


class LigandBaseProtocol; // fwd declaration
typedef utility::pointer::shared_ptr< LigandBaseProtocol > LigandBaseProtocolOP;
typedef utility::pointer::shared_ptr< LigandBaseProtocol const > LigandBaseProtocolCOP;

/// @brief Shared functionality for protocols that dock ligands.
///
/// @details Includes score function setup, interface definitions, and ligand flexibility.
///
class LigandBaseProtocol : virtual public protocols::moves::Mover
{
public:

	LigandBaseProtocol();

	virtual ~LigandBaseProtocol();

	core::scoring::ScoreFunctionOP  get_scorefxn();
	core::scoring::ScoreFunctionCOP get_scorefxn() const;

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	core::Size
	get_ligand_jump_id(
		core::pose::Pose const & pose
	) const;

	core::Size
	get_ligand_id(
		core::pose::Pose const & pose
	) const;

  core::Size
  get_ligand_id(
    core::pose::Pose const & pose,
    core::Size jump_id
  ) const;

    void
  restrain_protein_Calphas(
    core::pose::Pose & pose,
    utility::vector1< bool > const & is_restrained,
      //core::Real stddev_Angstroms,
    core::scoring::func::FuncOP restr_func
  ) const;

    void
  reorder_foldtree_around_mobile_regions(
    core::pose::Pose & pose,
    core::Size const & jump_id,
    utility::vector1< bool > const & mobile_bb,
    core::Size const & lig_id
  ) const;

	void
	get_non_bb_clashing_rotamers(
		core::pose::Pose const & pose,
		core::Size seqpos,
		core::scoring::ScoreFunctionCOP scofx,
		utility::vector1< core::conformation::ResidueCOP > & accepted_rotamers
	) const;

protected:

	core::scoring::ScoreFunctionOP
	make_tweaked_scorefxn(
		std::string const & weights_tag,
		bool estat_exclude_protein,
		bool estat_upweight,
		bool hbonds_downweight
	);

	core::Vector choose_desired_centroid(
		core::pose::Pose const & pose,
		core::Size jump_id,
		utility::vector1< core::Vector >
	);

	void move_ligand_to_desired_centroid(
		core::pose::Pose & pose,
		core::Size jump_id,
		utility::vector1< core::Vector >  start_from_pts
	);

	void move_ligand_to_desired_centroid(
		core::pose::Pose & pose,
		core::Size jump_id,
		core::Vector desired_centroid
	);

	core::kinematics::MoveMapOP
	make_movemap(
		core::pose::Pose const & pose,
		core::Size jump_id,
		core::Real sc_padding,
		bool include_all_rsds,
		bool include_backbone,
		bool include_ligands,
		bool include_water
	) const;

	/// @brief Shared machinery for the next two
	core::pack::task::PackerTaskOP
	make_packer_task(
		core::pose::Pose const & pose,
		ObjexxFCL::FArray1D_bool const & allow_repack,
		bool ligand_protonation
	) const;

	/// @brief Receptor (interface?) plus ligand
	core::pack::task::PackerTaskOP
	make_packer_task(
		core::pose::Pose const & pose,
		int jump_id,
		core::Real sc_padding,
		bool include_all_rsds,
		bool ligand_protonation
	) const;

	/// @brief Just ligand, not the receptor
	core::pack::task::PackerTaskOP
	make_packer_task_ligand_only(
		core::pose::Pose const & pose,
		int jump_id,
		bool ligand_protonation
	) const;

	void
	find_interface_rsds(
		core::pose::Pose const & pose,
		int jump_id,
		core::Real padding,
		ObjexxFCL::FArray1D_bool & is_interface //< output
	) const;

	void
	find_interface_backbone(
		core::pose::Pose const & pose,
		int jump_id,
		core::Real cutoff_dist,
		utility::vector1< bool > & is_interface, //< output
		utility::vector1< bool > & is_around_interface //< output
	) const;

	core::scoring::constraints::ConstraintOP
	restrain_ligand_nbr_atom(
		core::pose::Pose & pose,
		core::Size lig_id,
		core::Real stddev_Angstroms
	) const;


	void
	setup_bbmin_foldtree(
		core::pose::Pose & pose,
		core::Size const & jump_id,
		core::Real cutoff_dist,
		core::Real stddev_Angstroms
	);


protected:

	bool use_soft_rep_;
	core::scoring::ScoreFunctionOP scorefxn_, hard_scorefxn_, soft_scorefxn_;
	core::Real sc_interface_padding_; //< if using subset of residues, extra distance to use
	core::Real bb_interface_cutoff_; //< if using subset of residues, absolute distance to use
	core::pack::rotamer_set::UnboundRotamersOperationOP unboundrot_;

}; // class LigandBaseProtocol


} // namespace ligand_docking
} // namespace protocols

#endif // INCLUDED_protocols_ligand_docking_LigandBaseProtocol_HH
