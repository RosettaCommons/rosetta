// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody.design/AntibodyDesignMoverGenerator.hh
/// @brief Handles modeling of the antibody.  Before, after, and during design
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_AntibodyDesignMoverGenerator_hh
#define INCLUDED_protocols_antibody_design_AntibodyDesignMoverGenerator_hh

#include <protocols/antibody/design/AntibodyDesignMoverGenerator.fwd.hh>
#include <protocols/antibody/design/AntibodySeqDesignTFCreator.fwd.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.fwd.hh>
#include <protocols/antibody/AntibodyInfo.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>


#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.fwd.hh>


#include <protocols/toolbox/task_operations/RestrictToLoopsAndNeighbors.hh>
#include <protocols/simple_moves/ChangeAndResetFoldTreeMover.fwd.hh>
#include <protocols/moves/ChangeFoldTreeMover.hh>
#include <protocols/moves/MoverApplyingMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <protocols/docking/DockMCMProtocol.hh>
#include <protocols/docking/DockingLowRes.hh>

#include <utility/vector1.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace antibody {
namespace design {
	

///@brief Generates pre-configured general MoverOPs for use in AntibodyDesign and modeling
/// Helper functions for setting up FoldTrees.
/// Can set the generated mover to apply, which will include the needed FoldTree and any changing of typeset.
/// FoldTree stuff will eventually be refactored in favor of the TopologyBroker/Environment system.
///
class AntibodyDesignMoverGenerator : public protocols::moves::MoverApplyingMover {
	
public:
	
	///Default constructor for RosettaScripts.  Do not use.
	AntibodyDesignMoverGenerator();
	
	AntibodyDesignMoverGenerator(AntibodyInfoCOP ab_info);
	
	virtual ~AntibodyDesignMoverGenerator();

	///@brief Set the last setup mover as the current mover to call during apply.
	/// Default is True.  See apply for more details.
	void
	set_as_mover( bool setting );
	
	///@brief Apply the last mover generated or run the set mover.
	/// Will do the following:
	/// 1) Run any set ChangeFoldTreeMover.  This is set by any generated class if needed.
	/// 2) Will call set_cutpoint_variants.  If there are none, no harm done.
	/// 3) Call the set mover's apply function.
	/// 4) Remove cutpoint variants and any coord csts left over from [ relax ].
	/// Foldtree calls will eventually be replaced with the Environment/Topology Broker system.
	///
	virtual void
	apply( core::pose::Pose & pose );
	
	///@brief Generate a ChangeAndResetFoldTreeMover from any set mover, FoldTree mover, and ScoreFunction.
	/// Be careful of needed cutpoint variants.  If needed, generate everything manually!!
	simple_moves::ChangeAndResetFoldTreeMoverOP
	generate_protocol_mover() const;
	
	
	void
	set_defaults();
	
	
	///@brief Set to model CDRs.  Default is all of them false.
	void
	set_cdr( CDRNameEnum const cdr, bool setting );
	
	///@brief Set to model only one cdr, or all others but one.
	void
	set_cdr_only( CDRNameEnum const cdr, bool setting );
	
	void
	set_cdrs( utility::vector1<bool> cdrs );
	
	
	void
	set_scorefunction( core::scoring::ScoreFunctionCOP scorefxn );
	
	
	///@brief Set the mover for apply
	/// Generator functions will set this themselves, but can be changed before apply.
	virtual void
	set_mover( protocols::moves::MoverOP mover );
	
	///@brief Get the currently set mover
	virtual protocols::moves::MoverOP
	mover() const;
	
	
	///@brief Set the particular ChangeFoldTreeMover to use.
	/// Generator functions will set this themselves or to NULL, but can be changed before apply.
	void
	set_ft_mover( protocols::moves::ChangeFoldTreeMoverOP ft_mover );
	
	protocols::moves::ChangeFoldTreeMoverCOP
	ft_mover() const;
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////
	// Design Settings -  Use AntibodySeqDesignTFCreator for now.
	//
	//
	
	///@brief
	/// Override the TaskFactory that will be created and use this instead (for everything that uses a task!)
	///@details
	/// Hack for graft cdrs extremely specialized tf for design and to limit the amount of casting that would be done otherwise.
	/// Use of AntibodySeqDesignTFCreator or helper utility functions for TF creation are recommended.
	/// A better way here is preferred.
	void
	set_task_factory( core::pack::task::TaskFactoryOP tf );
	
	///////////////////////////////////////////////////////////////////////////////////////
	// Distance Detection
	//
	//
	
	///@brief Set a value for interface detection across the class.
	void
	interface_detection_dis( core::Real interface_distance );
	
	core::Real
	interface_detection_dis() const {
		return interface_dis_;
	}
	
	///@brief Set a value for neighbor detection across the class.
	void
	neighbor_detection_dis( core::Real neighbor_distance );
	
	core::Real
	neighbor_detection_dis() const {
		return neighbor_dis_;
	}
	
	
	///////////////////////////////////////////////////////////////////////////////////////
	// Protocol Level options
	//
	//
	
	///@brief How many residues on either side of the CDR to include in CDR - modeling? (2 is a good number - more if not North)
	void
	stem_size( core::Size overhang );
	
	///@brief Set the ab_dock chains, which will be used for any docking or minimization on the jump.  (LH_A, L_H, etc. - A for antigen)
	void
	ab_dock_chains( std::string ab_dock_chains );	
	
	
	///////////////////////////////////////////////////////////////////////////////////////
	// Additional Protocol Level options
	//
	//
	
	///@brief Set to minimize the set interface (during minimization/relax)
	/// Not implemented while minimizing CDRs as well.
	void
	set_min_interface( bool min_interface );
	
	///@brief Set to minimize sidechains (add to movemap)
	void
	set_min_sc( bool min_sc );
	
	///@brief Set to include neighbor side chains during any regional repacking
	void
	set_include_neighbor_sc( bool include_neighbor_sc );
	
	///@brief Use cartesian minimization for MinMovers
	/// Setup appropriate scorefunction and min settings
	void
	set_cartmin( bool cartmin );
	
	///@brief Use DualSpace relax during FastRelax
	/// Setup appropriate scorefunction and settings
	void
	set_dualspace( bool dualspace );
	
	///@brief Set to use starting coordinate constraints during relax
	void
	set_start_coord_csts( bool coord_csts );
	
	
	void
	set_dock_low_res_scorefunction( core::scoring::ScoreFunctionCOP scorefxn );
	
	void
	set_dock_high_res_scorefunction( core::scoring::ScoreFunctionCOP scorefxn );
	
	///@brief Set the outer cycles for DockMCM.
	/// Normal in docking is 4.
	void
	set_high_res_dock_outer_cycles( core::Size first_cycle = 3 );
	
	///@brief Set the inner cycles for DockMCM.
	/// Normal in docking is 45.
	void
	set_high_res_dock_inner_cycles( core::Size second_cycle = 10 );
	
	
	///////////////////////////////////////////////////////////////////////////////////////
	// FoldTree setup - will be deprecated in favor of Topology Broker
	//
	//
	
	void
	setup_general_min_foldtree(core::pose::Pose const & pose, protocols::moves::ChangeFoldTreeMover & ft_mover);
	
	void
	setup_dock_foldtree(core::pose::Pose const & pose, protocols::moves::ChangeFoldTreeMover & ft_mover);
	
	
	///////////////////////////////////////////////////////////////////////////////////////
	// Minimization
	//
	//
	
	
	///@brief Generate a packer and optionally set as the current mover
	simple_moves::PackRotamersMoverOP
	generate_repack_cdrs( core::pose::Pose const & pose );
	
	//@brief Setup a packer to repack CDRs set.
	void
	setup_repack_cdrs( core::pose::Pose const & pose, simple_moves::PackRotamersMoverOP packer );

	
	///@brief Get a pre-configured Packer to repack the interface between Antibody and Antigen.
	/// If set_as_mover is true, will set the needed ab-ag foldtree to this class which the task op will use.
	protocols::simple_moves::PackRotamersMoverOP
	generate_repack_antigen_ab_interface( Pose const & pose );
	
	///@brief Configure a packer to repack the interface between Antibody and Antigen
	/// If set_as_mover is true, will set the needed ab-ag foldtree to this class which the task op will use.
	void
	setup_repack_antigen_ab_interface( core::pose::Pose const & pose, simple_moves::PackRotamersMoverOP packer );
	
	
	///@brief Vanilla minimizer using dfpmin_armijo_nonmonotone at .001 tolerance (or lfbgs for cartmin).
	/// If set_as_mover is true, will set the appropriate FoldTree to this class.
	simple_moves::MinMoverOP
	generate_minimizer( core::pose::Pose const & pose );
	
	///@brief Vanilla minimizer using dfpmin_armijo_nonmonotone at .001 tolerance (or lfbgs for cartmin).
	/// If set_as_mover is true, will set the appropriate FoldTree to this class.
	void
	setup_minimizer( core::pose::Pose const & pose, simple_moves::MinMoverOP min_mover);
	
	
	///@brief Generate FastRelax of CDRs +/or Interface.  
	///@details Any Cluster Constraints should already be set. Optionally use start_coordinate constraints.  
	/// All coordinate constraints on the pose will then be removed in this classes apply.
	relax::FastRelaxOP
	generate_relax( core::pose::Pose const & pose );
	
	///@brief Setup FastRelax of CDRs +/or Interface.  
	///@details Any Cluster Constraints should already be set. Optionally use start_coordinate constraints.  
	/// All coordinate constraints on the pose will then be removed in this classes apply.
	void
	setup_relax( core::pose::Pose const & pose, relax::FastRelaxOP rel);
	
	
	///////////////////////////////////////////////////////////////////////////////////////
	// Docking
	//
	//
	
	///@brief Configure a low-res docking mover.  No repacking, but will switch chains.
	/// Using A for designation of Antigen chains. Not full antigen chains. (L_A, LH_A, L_H, etc.)
	docking::DockingLowResOP
	generate_dock_low_res( core::pose::Pose const & pose );
	
	void
	setup_dock_low_res( core::pose::Pose const & pose, docking::DockingLowResOP docker );
	
	
	///@brief Configure a high-res, DockMCMProtocol for docking.
	docking::DockMCMProtocolOP
	generate_dock_high_res( core::pose::Pose const & pose );
	
	void
	setup_dock_high_res( core::pose::Pose const & pose, docking::DockMCMProtocolOP docker );
	
		
public:
	
	///////////////////////////////////////////////////////////////////////////////////////
	// Loops, Foldtrees, and Movemaps
	//
	//
	
	///@brief Set a range of CDRs.
	void
	set_cdr_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool setting);
	
	core::kinematics::MoveMapOP
	get_cdrs_movemap_with_overhang(core::pose::Pose const & pose, bool min_bb = true, bool min_sc = true, bool include_neighbor_sc = true, bool include_neighbor_bb = false) const;
		
	core::kinematics::MoveMapOP
	get_movemap_from_task(core::pose::Pose const & pose, core::pack::task::PackerTaskCOP task) const;
	
	
public:
	virtual std::string
	get_name() const {
		return "AntibodyDesignMoverGenerator";
	}
	
private:
	
	void
	read_command_line_options();

	void
	setup_scorefxns();
	
	void
	setup_task_operations();
	
	core::pack::task::TaskFactoryOP
	get_dock_tf();
	
	AntibodyInfoCOP ab_info_;
	protocols::moves::MoverOP mover_;
	protocols::moves::ChangeFoldTreeMoverOP ft_mover_;
	
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreFunctionOP docking_scorefxn_low_;
	core::scoring::ScoreFunctionOP docking_scorefxn_high_;
	
	//TaskFactories
	core::pack::task::TaskFactoryOP interface_tf_; //Interface TF at jump 1  - Used as baseline for interface repacking/design
	core::pack::task::TaskFactoryOP tf_; //Basic TF that gets cleared upon use..
	core::pack::task::TaskFactoryOP set_tf_; //Set task factory to use - mainly for design.  Not ideal. Change.
	
	
	//TaskOperations.  So that they are not constructed every graft.
	protocols::toolbox::task_operations::RestrictToLoopsAndNeighborsOP loops_operation_;
	core::pack::task::operation::InitializeFromCommandlineOP cmd_line_operation_;
	core::pack::task::operation::RestrictToRepackingOP restrict_design_operation_;
	
	core::Real interface_dis_;
	core::Real neighbor_dis_;
	
	std::string ab_dock_chains_; //L_H, LH_A, etc.  Antigen as A, not full antigen chain strings
	utility::vector1 <bool> model_cdrs_; //These cdrs are cdr's to work on.  Whether copying, refining, designing, etc.etc.
	core::Size overhang_;
            
	//bool design_; Unused at the moment
	//bool design_neighbors_;
	//bool design_neighbor_cdrs_;
	//bool design_neighbor_framework_;
	

	
	bool min_interface_;
	bool min_sc_;
	bool include_neighbor_sc_;
	bool cartmin_;
	bool dualspace_;
	bool start_coord_csts_;
	
	core::Size high_res_dock_outer_cycles_;
	core::Size high_res_dock_inner_cycles_;
	
	bool set_as_mover_;
};
}
}
}



#endif //INCLUDED_protocols_antibody_design_AntibodyDesignMover_hh
