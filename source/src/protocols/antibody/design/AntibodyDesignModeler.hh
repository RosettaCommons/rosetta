// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody.design/AntibodyDesignModeler.hh
/// @brief Handles modeling of the antibody.  Before, after, and during design
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_AntibodyDesignModeler_hh
#define INCLUDED_protocols_antibody_design_AntibodyDesignModeler_hh

#include <protocols/antibody/design/AntibodyDesignModeler.fwd.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.fwd.hh>

#include <protocols/toolbox/task_operations/RestrictToLoopsAndNeighbors.hh>

#include <utility/vector1.hh>
#include <core/scoring/ScoreFunction.hh>
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

///@brief Class that runs the modeling for AntibodyDesign protocols.  Could be moved to antibody.
class AntibodyDesignModeler : public utility::pointer::ReferenceCount {
            
public:
	AntibodyDesignModeler(AntibodyInfoOP ab_info);
	virtual ~AntibodyDesignModeler();

	///@brief Set to model CDRs.  Default is all of them false.
	void
	set_cdr(CDRNameEnum const cdr, bool setting);
	
	///@brief Set to model only one cdr, or all others but one.
	void
	set_cdr_only(CDRNameEnum const cdr, bool setting);
	
	///@brief Set a range of CDRs.
	void
	set_cdr_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool setting);
	
	
	void
	set_scorefunction(ScoreFunctionCOP scorefxn);
	
	
	///@brief Set a value for interface detection across the class.
	void
	set_interface_detection_dis(core::Real interface_distance);
	
	void
	set_neighbor_detection_dis(core::Real neighbor_distance);
	
	
	///////////////////////////////////////////////////////////////////////////////////////
	// Minimization
	//
	//
	
	///@brief Relax CDRs only either using FastRelax or CentroidRelax
	///@details Cluster Constraints should already be set. Optionally use start_coordinate constraints.  
	/// All coordinate constraints on the pose will then be removed.
	void
	relax_cdrs(Pose & pose, bool centroid_mode, bool starting_coordinate_constraints = false) const;
	
	void
	relax_cdrs_and_neighbor_sc(Pose & pose, bool starting_coordinate_constraints = false) const;
	
	
	///@brief Vanilla minimizer using dfpmin_armijo_nonmonotone at .01 tolerance.
	void
	minimize_cdrs(Pose & pose) const;
	
	void
	minimize_cdrs_and_neighbor_sc(Pose & pose) const;

	
	///@brief Repack the interface between Antibody and Antigen.  Foldtree for docking (LH_A/etc.) must be set.
	void
	repack_antigen_ab_interface(Pose & pose) const;

	///@brief Repack the interface between antibody and antigen, but only pack antigen.
	void
	repack_antigen_interface(Pose & pose) const;
	
	///@brief Repack the interface between antibody and antigen, but only pack antibody.
	void
	repack_antibody_interface(Pose & pose) const;
	
	
	///@brief Repack the CDRs given.  Nothing special.
	void
	repack_CDRs(Pose & pose);

	///@brief Repack CDRs and their neighbors within distance.
	void
	repack_CDRs_and_neighbors(Pose & pose);

	
	///////////////////////////////////////////////////////////////////////////////////////
	// Docking
	//
	//
	
	void
	dock_LH_A_low_res(Pose & pose, bool pack_interface = false) const;
	
	void
	dock_LH_A_high_res(Pose & pose, int first_cycle=4, int second_cycle=45) const;
		
	
	void
	dock_A_LH_low_res(Pose & pose, bool pack_interface = false) const;

	void
	dock_A_LH_high_res(Pose & pose, int first_cycle=4, int second_cycle=45) const;
	
	
	void
	dock_L_H_low_res(Pose & pose, bool pack_interface = false) const;

	void
	dock_L_H_high_res(Pose & pose, int first_cycle=4, int second_cycle=45) const;
	
	
	void
	dock_H_L_low_res(Pose & pose, bool pack_interface = false) const;
	
	void
	dock_H_L_high_res(Pose & pose, int first_cycle=4, int second_cycle=45) const;
	

	

	

	
	///////////////////////////////////////////////////////////////////////////////////////
	// Design
	//
	//
	
protected:
            
	void
	read_command_line_options();

	///@brief Applies an LH-A Foldtree.  Forked from AbInfo - Its version uses hardcoded order of the PDB.  Need to refactor that using this.
	/// Need to find out if they are exactly the same first.  
	void
	apply_LH_A_foldtree(core::pose::Pose & pose) const;
	
	///@brief AntibodyInfo has this function.  Need to refactor more to get single loop on the fly
	protocols::loops::LoopsOP
	get_cdr_loops(Pose & pose) const;
	
	AntibodyInfoOP ab_info_;
	ScoreFunctionOP scorefxn_;
	ScoreFunctionOP docking_scorefxn_high_;
	
	//TaskFactories
	TaskFactoryOP antigen_interface_tf_; //Antigen-Antibody interface.  Used as baseline for interface repacking/design
	TaskFactoryOP tf_; //Basic TF that gets cleared upon use..
	
	
	//TaskOperations.  So that they are not constructed every graft.
	protocols::toolbox::task_operations::RestrictToLoopsAndNeighborsOP loops_operation_;
	operation::InitializeFromCommandlineOP cmd_line_operation_;
	operation::RestrictToRepackingOP restrict_design_operation_;
	
	core::Real interface_dis_;
	core::Real neighbor_dis_;
	
	vector1 <bool> cdrs_; //These cdrs are cdr's to work on.  Whether copying, refining, designing, etc.etc.

            

            
};
}
}
}



#endif //INCLUDED_protocols_antibody_design_AntibodyDesignMover_hh
