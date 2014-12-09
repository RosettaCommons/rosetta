// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AntibodyModeler
/// @brief Build a homology model of an antibody
/// @detailed
///
///
/// @author Aroop Sircar


#ifndef INCLUDED_protocols_antibody_legacy_AntibodyModeler_hh
#define INCLUDED_protocols_antibody_legacy_AntibodyModeler_hh

#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#include <protocols/antibody_legacy/AntibodyClass.hh>
#include <protocols/antibody_legacy/AntibodyModeler.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace antibody_legacy {

class AntibodyModeler: public moves::Mover {
public:

	// default constructor
	AntibodyModeler();

	// default destructor
	~AntibodyModeler();

	void set_default();

	void init_from_options();

	virtual void apply( core::pose::Pose & pose_in );
	virtual std::string get_name() const;

	virtual protocols::moves::MoverOP clone() const;

	virtual	void init_on_new_input();

	void set_snugdock_foldtree( core::pose::Pose & pose_in );

	void setup_simple_fold_tree(
	    core::Size jumppoint1,
	    core::Size cutpoint,
	    core::Size jumppoint2,
	    core::Size nres,
	    core::pose::Pose & pose_in );

	void relax_cdrs();

	void all_cdr_VL_VH_fold_tree(
	    core::pose::Pose & pose_in,
	    const loops::Loops & loops );

	void repulsive_ramp(
	    core::pose::Pose & pose_in,
	    loops::Loops loops_in );

	void snugfit_MC_min (
	    core::pose::Pose & pose_in,
	    core::kinematics::MoveMapOP cdr_dock_map,
	    core::Size cycles,
	    core::Real minimization_threshold,
	    core::scoring::ScoreFunctionOP scorefxn,
	    core::scoring::ScoreFunctionOP pack_scorefxn,
	    utility::vector1< bool> is_flexible );

	void snugfit_mcm_protocol(
	    core::pose::Pose & pose_in,
	    loops::Loops loops_in );

	void setup_packer_task( core::pose::Pose & pose_in );

	core::Real global_loop_rmsd (
	    const core::pose::Pose & pose_in,
	    const core::pose::Pose & native_pose,
	    std::string cdr_type );

	void read_and_store_fragments();

	void display_constraint_residues();


private:

	// Modeling H3 options
	bool model_h3_;
	bool snugfit_;
	bool native_present_;
	bool graft_l1_;
	bool graft_l2_;
	bool graft_l3_;
	bool graft_h1_;
	bool graft_h2_;
	bool graft_h3_;
	bool camelid_;
	bool camelid_constraints_;

	// Benchmark mode for shorter_cycles
	bool benchmark_;

	// flag for one time fragment initialization
	bool init_for_input_yet_;

	// used as a flag to enable reading in of cst files
	core::Real cst_weight_;

	std::string native_filename_;
	core::pose::Pose start_pose_;
	core::pose::Pose native_pose_;
	Antibody antibody_in_;
	utility::vector1< core::fragment::FragData > H3_base_library_;
	utility::vector1< core::fragment::FragSetOP > offset_frags_;

	//packer task
	core::pack::task::TaskFactoryOP tf_;
	core::pack::task::TaskFactoryOP init_task_factory_;

	std::map < std::string, core::Real > score_map_;

}; // class AntibodyModeler

} // namespace antibody
} // namespace protocols
#endif
