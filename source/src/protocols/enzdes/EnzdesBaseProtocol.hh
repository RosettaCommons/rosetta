// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/enzdes/EnzdesBaseProtocol.hh
///
/// @brief
/// @author Florian Richter




#ifndef INCLUDED_protocols_enzdes_EnzdesBaseProtocol_hh
#define INCLUDED_protocols_enzdes_EnzdesBaseProtocol_hh


#include <protocols/enzdes/EnzdesBaseProtocol.fwd.hh>
// AUTO-REMOVED #include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh>

#include <protocols/ligand_docking/LigandBaseProtocol.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <utility/vector1.hh>



namespace protocols {
namespace enzdes {


class EnzdesBaseProtocol : public protocols::ligand_docking::LigandBaseProtocol
{

public:

	EnzdesBaseProtocol();
	virtual std::string get_name() const;

	//virtual void apply( core::pose::Pose & pose) const = 0;

	//toolbox::match_enzdes_util::EnzConstraintIOOP cst_io(){ return cst_io_; }

 //catalytic res INCLUDING all ligands in pose numbering in a particular order
	utility::vector1<Size> catalytic_res( core::pose::Pose const & pose) const;

	std::set< Size > const & design_targets( core::pose::Pose const & pose ) const;

	bool
	is_catalytic_position( core::pose::Pose const & pose, core::Size const seqpos ) const;

	static void register_options();

	core::chemical::ResidueTypeSetCAP
	restype_set() const {
		return restype_set_; }

	void
	generate_explicit_ligand_rotamer_poses(
		core::pose::Pose const & orig_pose,
		utility::vector1< core::pose::PoseOP > & ligrot_poses,
		core::scoring::ScoreFunctionCOP scofx
	);

	core::scoring::ScoreFunctionCOP
	reduced_scorefxn() const;

	core::scoring::ScoreFunctionOP
  reduced_scorefxn();

	core::Real
	design_targets_score(
		core::pose::Pose const & pose
	) const;


	virtual
	void
	remap_resid(
		core::pose::Pose const & pose,
		core::id::SequenceMapping const & smap
	);

	void
	remove_enzdes_constraints(
		core::pose::Pose & pose,
		bool keep_covalent
	) const;

	void
	add_pregenerated_enzdes_constraints(
		core::pose::Pose & pose
	) const;

	void
	cst_minimize(
		core::pose::Pose & pose,
		core::pack::task::PackerTaskCOP task,
		bool cst_opt = false
	) const;

	core::pack::task::PackerTaskOP
	create_enzdes_pack_task(
		core::pose::Pose & pose,
		bool design = true
	); //the task

	void
	setup_sequence_recovery_cache(
		core::pose::Pose & pose,
		core::pack::task::PackerTask const & task
	) const;

	void
	set_all_jumps_minimizable( bool const & setting ){ min_all_jumps_ = setting; }

	void
	set_minimize_options(
		bool const & min_sc,
		bool const & min_bb,
		bool const & min_rb,
		bool const & min_lig,
		bool backrub = false){
	chi_min_ = min_sc; bb_min_= min_bb; rb_min_=min_rb; minimize_ligand_torsions_=min_lig; bb_backrub_ =backrub;
	}

	void
	set_fix_cataa(
		bool const & setting){ fix_catalytic_aa_ = setting;}

	core::kinematics::MoveMapOP
	create_enzdes_movemap(
		core::pose::Pose & pose,
		core::pack::task::PackerTaskCOP task,
		bool min_all_jumps= false
	) const;

	void
	enzdes_pack(
		core::pose::Pose & pose,
		core::pack::task::PackerTaskCOP,
		core::scoring::ScoreFunctionCOP scorefxn,
		core::Size cycles,
		bool minimize_after_packing,
		bool pack_unconstrained,
		bool favor_native
	) const;

	void set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn );

	void
	setup_bbmin_ft_and_csts(
		core::pose::Pose & pose,
		utility::vector1< bool > allow_move_bb,
		core::Size jump_id
	) const;
protected:
	void
	setup_enzdes_constraints(
		core::pose::Pose & pose,
		bool allow_missing_remark_blocks
	) const;

	void enable_constraint_scoreterms();

	/// @brief function to disable constraint scoring terms:
	/// @NOTE: this will leave eventual covalent connections set up by EnzConstraintIO untouched.
	void disable_constraint_scoreterms();

	bool
	exchange_ligands_in_pose(
		core::pose::Pose & pose,
		bool check_bb_clashes,
		core::scoring::ScoreFunctionCOP scofx
	);

	core::scoring::ScoreFunctionOP reduced_sfxn_;

	mutable std::set< core::Size > design_targets_;

	bool include_all_design_targets_in_design_interface_;

private:

	void
	read_ligand_superposition_file( std::string filename );

	//data

	//toolbox::match_enzdes_util::EnzConstraintIOOP cst_io_;

	core::chemical::ResidueTypeSetCAP restype_set_;
	core::scoring::EnergyMap constraint_weights_;

	core::Real bb_min_allowed_dev_;
	core::Real loop_bb_min_allowed_dev_;

	//stuff for optional ligand superposition
	bool lig_superposition_file_read_;
	std::pair< std::string, std::string > res_to_superimpose_;
	utility::vector1< std::pair< std::string, std::string > > atoms_to_superimpose_on_;
	bool chi_min_, bb_min_, rb_min_, minimize_ligand_torsions_, bb_backrub_;
	bool min_all_jumps_; // minimize protein jumps too?
	bool minimize_all_ligand_torsions_;
	core::Real lig_min_stddev_;

	bool exclude_protein_protein_fa_elec_, fix_catalytic_aa_;

}; //class EnzdesBaseProtocol


} //namespace enzdes
} //namespace protocols




#endif // INCLUDED_protocols_enzdes_EnzdesBaseProtocol_HH
