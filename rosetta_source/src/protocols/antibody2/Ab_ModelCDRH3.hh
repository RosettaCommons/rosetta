// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file Ab_ModelCDRH3
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody2_Ab_ModelCDRH3_hh
#define INCLUDED_protocols_antibody2_Ab_ModelCDRH3_hh

#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#include <protocols/antibody2/Ab_Info.hh>
#include <protocols/antibody2/Ab_ModelCDRH3.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/antibody2/CDRH3Modeler2.fwd.hh>


#include <utility/vector1.hh>

namespace protocols {
namespace antibody2 {

class Ab_ModelCDRH3: public moves::Mover {
public:

	// default constructor
	Ab_ModelCDRH3();

	// default destructor
	~Ab_ModelCDRH3();

	virtual protocols::moves::MoverOP clone() const;

	/// @brief Assigns default values to primitive members
	void set_default();

	/// @brief Instantiates non-primitive members based on the value of the primitive members
	void sync_objects_with_flags();

	virtual void apply( core::pose::Pose & pose );

	// simple inline setters
	void set_model_h3( bool model_h3 ) { model_h3_ = model_h3; }
	void set_snugfit( bool snugfit ) { snugfit_ = snugfit; }
	void set_graft_l1( bool graft_l1 ) { graft_l1_ = graft_l1; }
	void set_graft_l2( bool graft_l2 ) { graft_l2_ = graft_l2; }
	void set_graft_l3( bool graft_l3 ) { graft_l3_ = graft_l3; }
	void set_graft_h1( bool graft_h1 ) { graft_h1_ = graft_h1; }
	void set_graft_h2( bool graft_h2 ) { graft_h2_ = graft_h2; }
	void set_graft_h3( bool graft_h3 ) { graft_h3_ = graft_h3; }
	void set_camelid( bool camelid ) { camelid_ = camelid; }
	void set_camelid_constraints( bool camelid_constraints ) { camelid_constraints_ = camelid_constraints; }
	void set_benchmark( bool benchmark ) { benchmark_ = benchmark; }

	virtual std::string get_name() const;

	void setup_simple_fold_tree(
		core::Size jumppoint1,
		core::Size cutpoint,
		core::Size jumppoint2,
		core::Size nres,
		core::pose::Pose & pose_in );

	void relax_cdrs( core::pose::Pose & pose );

	void all_cdr_VL_VH_fold_tree( core::pose::Pose & pose_in, const loops::Loops & loops );

	void repulsive_ramp( core::pose::Pose & pose_in, loops::Loops loops_in );

	void snugfit_MC_min (
		core::pose::Pose & pose_in,
		core::kinematics::MoveMapOP cdr_dock_map,
		core::Size cycles,
		core::Real minimization_threshold,
		core::scoring::ScoreFunctionOP scorefxn,
	  core::scoring::ScoreFunctionOP pack_scorefxn,
		utility::vector1< bool> is_flexible );

	void snugfit_mcm_protocol(core::pose::Pose & pose_in, loops::Loops loops_in );

	void setup_packer_task( core::pose::Pose & pose_in );

	core::Real global_loop_rmsd ( const core::pose::Pose & pose_in, const core::pose::Pose & native_pose, std::string cdr_type );

	void read_and_store_fragments( core::pose::Pose & pose );

	void display_constraint_residues( core::pose::Pose & pose );

    
    
    void show( std::ostream & out=std::cout );
    friend std::ostream & operator<<(std::ostream& out, const Ab_ModelCDRH3 & ab_m_2 );
    
    
    /// @brief Associates relevant options with the AntibodyModeler class
    static void register_options();
    


public:

	// Modeling H3 options
	bool model_h3_;
	bool snugfit_;
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

	bool user_defined_; // for constructor options passed to init

	// flag for one time fragment initialization
	bool flags_and_objects_are_in_sync_;
	bool first_apply_with_current_setup_;

	// used as a flag to enable reading in of cst files
	core::Real cst_weight_;

	// score functions
	core::scoring::ScoreFunctionOP scorefxn_;

	// external objects
	Ab_Info ab_info_;
	utility::vector1< core::fragment::FragSetOP > offset_frags_;

	//packer task
	core::pack::task::TaskFactoryOP tf_;
	core::pack::task::TaskFactoryOP init_task_factory_;

	// movers

	protocols::antibody2::CDRH3Modeler2OP model_cdrh3_;

	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();


	/// @brief Performs the portion of setup of non-primitive members that requires a pose - called on apply
	void finalize_setup( core::pose::Pose & pose );

	/// @brief Sets up the instance of AntibodyModeler and initializes all members based on values passed in at construction
	///		or via the command line.
	void init();

	void setup_objects();

}; // class Ab_ModelCDRH3

    
    
    
    
} // namespace antibody2
} // namespace protocols

#endif
