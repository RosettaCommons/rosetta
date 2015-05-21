// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody.design/AntibodyDesignModeler.hh
/// @brief Handles modeling of the antibody.  Before, after, and during design
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_AntibodyDesignModeler_hh
#define INCLUDED_protocols_antibody_design_AntibodyDesignModeler_hh

#include <protocols/antibody/design/AntibodyDesignModeler.fwd.hh>
#include <protocols/antibody/design/AntibodySeqDesignTFCreator.fwd.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.fwd.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/toolbox/task_operations/RestrictToLoopsAndNeighbors.hh>

#include <utility/vector1.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace antibody {
namespace design {
	
using namespace core::scoring;
using namespace core::pose;
using namespace core::pack::task;
using std::string;
using core::Real;
using utility::vector1;

/// @brief Basic Class that can run modeling for various protocols.
/// Does not do any design.
/// Mainly deprecated by AntibodyDesignMoverGenerator for the new class.
class AntibodyDesignModeler : public utility::pointer::ReferenceCount {
            
public:
	AntibodyDesignModeler(AntibodyInfoOP ab_info);
	virtual ~AntibodyDesignModeler();


	
	
	/// @brief Set to model CDRs.  Default is all of them false.
	void
	set_cdr(CDRNameEnum const cdr, bool setting);
	
	/// @brief Set to model only one cdr, or all others but one.
	void
	set_cdr_only(CDRNameEnum const cdr, bool setting);
	
	/// @brief Set a range of CDRs.
	void
	set_cdr_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool setting);
	
	void
	set_cdrs(utility::vector1<bool> cdrs);
	
	void
	set_scorefunction(ScoreFunctionCOP scorefxn);
	
	
	////////////////////////////////////////////////////////////////////////
	
	/// @brief How many residues on either side of the CDR to include in CDR - modeling? (Default is 2)
	void
	cdr_overhang(CDRNameEnum const cdr, core::Size const overhang);
	
	/// @brief Set the ab_dock chains, which will be used for any docking or minimization on the jump.  (LH_A, L_H, etc. - A for antigen)
	void
	ab_dock_chains(std::string ab_dock_chains);	
	
	
	/// @brief Set a value for interface detection across the class.
	void
	interface_detection_dis(core::Real interface_distance);
	
	core::Real
	interface_detection_dis() const {
		return interface_dis_;
	}
	
	void
	neighbor_detection_dis(core::Real neighbor_distance);
	
	core::Real
	neighbor_detection_dis() const {
		return neighbor_dis_;
	}
	
	///////////////////////////////////////////////////////////////////////////////////////
	// Minimization
	//
	//
	
	/// @brief Repack the CDRs given.  Nothing special.
	void
	repack_cdrs(Pose & pose, bool include_neighbor_sc = true);
	
	/// @brief Vanilla minimizer using dfpmin_armijo_nonmonotone at .001 tolerance (or lfbgs for cartmin).
	void
	minimize_cdrs(Pose & pose,
		bool min_sc = true,
		bool include_neighbor_sc = true,
		bool min_interface = false,
		bool cartmin = false,
		bool use_talaris_cartmin = false) const;
	

	/// @brief Relax CDRs using FastRelax.  CentroidRelax unfortunately does not work well.  Perhaps with centroid rotamers..
	/// @details Cluster Constraints should already be set. Optionally use start_coordinate constraints.  
	/// All coordinate constraints on the pose will then be removed.
	void
	relax_cdrs(Pose & pose, bool include_neighbor_sc = true,  bool starting_coordinate_constraints = false,  bool min_interface = false, bool dualspace = false) const;
	
	
	void
	minimize_interface(Pose & pose, bool min_interface_sc = true) const;
	
	//@brief Relax torsion (jump) between chains in dock_chains, optionally 
	void
	relax_interface(Pose & pose, bool min_interface_sc = true) const;
	

	/// @brief Repack the interface between Antibody and Antigen. 
	void
	repack_antigen_ab_interface(Pose & pose) const;

	/// @brief Repack the interface between antibody and antigen, but only pack antigen.
	void
	repack_antigen_interface(Pose & pose) const;
	
	/// @brief Repack the interface between antibody and antigen, but only pack antibody.
	void
	repack_antibody_interface(Pose & pose) const;
	
	
	///////////////////////////////////////////////////////////////////////////////////////
	// Docking
	//
	//
	
	/// @brief Dock the antibody. Using A for designation of Antigen chains. Not full antigen chains. (L_A, LH_A, L_H, etc.)
	void
	dock_low_res(Pose & pose, bool pack_interface = false) const;
	
	void
	dock_high_res(Pose & pose, core::Size first_cycle=4, core::Size second_cycle = 45) const;
	
	
	void
	set_dock_low_res_scorefunction(ScoreFunctionCOP scorefxn);
	
	void
	set_dock_high_res_scorefunction(ScoreFunctionCOP scorefxn);
	

	///////////////////////////////////////////////////////////////////////////////////////
	// Benchmarking
	//

	/// @brief Randomizes the dihedrals of the CDR.
	void
	extend_CDR(Pose & pose, CDRNameEnum cdr) const;
	
	
	
public:
	
	///////////////////////////////////////////////////////////////////////////////////////
	// Loops, Foldtrees, and Movemaps
	//
	//
	
	/// @brief Applies a docking foldtree..  Forked from AbInfo - Its version uses hardcoded order of the PDB.  However, this calls docking FT, which is hard coded as well for now.
	/// Fix Me.
	void
	apply_LH_A_foldtree(core::pose::Pose & pose) const;
	
	
	/// @brief Get CDR loops with cutpoint in the middle.
	protocols::loops::LoopsOP
	get_cdr_loops(Pose const & pose) const;
	
	/// @brief Get CDR loops with set overhang with cutpoint at the stop position +1.
	protocols::loops::LoopsOP
	get_cdr_loops_with_overhang(Pose const & pose) const ;
	
	protocols::loops::Loop
	get_cdr_loop_with_overhang(Pose const & pose, CDRNameEnum cdr) const;
	
	
	core::kinematics::MoveMapOP
	get_cdrs_movemap_with_overhang(Pose const & pose, bool min_bb = true, bool min_sc = true, bool include_neighbor_sc = true, bool include_neighbor_bb = false) const;
		
	core::kinematics::MoveMapOP
	get_movemap_from_task(Pose const & pose, core::pack::task::PackerTaskCOP task) const;
	
private:
	void
	set_defaults();
	
	void
	read_command_line_options();

	void
	setup_scorefxns();
	
	void
	setup_task_operations();
	
	AntibodyInfoOP ab_info_;
	ScoreFunctionOP scorefxn_;
	ScoreFunctionOP docking_scorefxn_low_;
	ScoreFunctionOP docking_scorefxn_high_;
	
	//TaskFactories
	TaskFactoryOP interface_tf_; //Interface TF at jump 1  - Used as baseline for interface repacking/design
	TaskFactoryOP tf_; //Basic TF that gets cleared upon use..
	
	
	//TaskOperations.  So that they are not constructed every graft.
	protocols::toolbox::task_operations::RestrictToLoopsAndNeighborsOP loops_operation_;
	operation::InitializeFromCommandlineOP cmd_line_operation_;
	operation::RestrictToRepackingOP restrict_design_operation_;
	
	core::Real interface_dis_;
	core::Real neighbor_dis_;
	
	std::string ab_dock_chains_; //L_H, LH_A, etc.  Antigen as A, not full antigen chain strings
	vector1 <bool> model_cdrs_; //These cdrs are cdr's to work on.  Whether copying, refining, designing, etc.etc.
	vector1 <core::Size> overhangs_;
            
	//bool design_;
	//bool design_neighbors_;

};
}
}
}


#endif //INCLUDED_protocols_antibody_design_AntibodyDesignMover_hh
