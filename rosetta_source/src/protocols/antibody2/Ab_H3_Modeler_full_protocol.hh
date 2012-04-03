// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/Ab_H3_Modeler_full_protocol.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#ifndef INCLUDED_protocols_antibody2_Ab_H3_Modeler_full_protocol_hh
#define INCLUDED_protocols_antibody2_Ab_H3_Modeler_full_protocol_hh

#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/PyMolMover.fwd.hh>

#include <protocols/antibody2/Ab_H3_Model_CDR_H3.fwd.hh>
#include <protocols/antibody2/Ab_Info.hh>
#include <protocols/antibody2/Ab_H3_Modeler_full_protocol.fwd.hh>


#include <utility/vector1.hh>

namespace protocols {
namespace antibody2 {

class Ab_H3_Modeler_full_protocol: public moves::Mover {
public:

	// default constructor
	Ab_H3_Modeler_full_protocol();

	// default destructor
	~Ab_H3_Modeler_full_protocol();

	virtual protocols::moves::MoverOP clone() const;

	/// @brief Assigns default values to primitive members
	void set_default();

	/// @brief Instantiates non-primitive members based on the value of the primitive members
	void sync_objects_with_flags();

	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;
    
    /// @brief Associates relevant options with the AntibodyModeler class
    static void register_options();
    
	// simple inline setters
    void set_h3modeler( bool model_h3 ){ 
        model_h3_ = model_h3; 
    }
	void set_snugfit( bool snugfit ) { 
        snugfit_ = snugfit; 
    }
	void set_camelid( bool camelid ) { 
        camelid_ = camelid; 
    }
	void set_camelid_constraints(bool camelid_constraints ) { 
        camelid_constraints_ = camelid_constraints; 
    }
	void set_benchmark( bool benchmark ) { 
        benchmark_ = benchmark; 
    }
    void set_cst_weight( core::Real const cst_weight){
        cst_weight_=cst_weight;
    }
    void set_H3_filter(bool H3_filter){
        H3_filter_ = H3_filter;
    }
    void set_cter_insert(bool cter_insert){
        cter_insert_ = cter_insert;
    }

	void relax_cdrs( core::pose::Pose & pose );

	void all_cdr_VL_VH_fold_tree( core::pose::Pose & pose_in, 
                                  const loops::Loops & loops );


	core::Real global_loop_rmsd ( const core::pose::Pose & pose_in, 
                                  const core::pose::Pose & native_pose, 
                                  std::string cdr_type );

	void display_constraint_residues( core::pose::Pose & pose );

        
    void show( std::ostream & out=std::cout );
    friend std::ostream & operator<<(std::ostream& out, const Ab_H3_Modeler_full_protocol & ab_m_2 );
    
    

private:
    bool model_h3_;
    bool extreme_repacking_;
	bool snugfit_;
    bool refine_h3_;
	bool camelid_;
	bool camelid_constraints_;
    bool H3_filter_;
    bool cter_insert_;
    core::pose::Pose start_pose_;
    
    /// @brief refine H3 only
	bool antibody_refine_;
    core::Real high_cst_;

    moves::PyMolMoverOP pymol_;
    bool use_pymol_diy_;

    
    
	// Benchmark mode for shorter_cycles
	bool benchmark_;

	bool user_defined_; // for constructor options passed to init

	bool flags_and_objects_are_in_sync_;
	bool first_apply_with_current_setup_;

	// used as a flag to enable reading in of cst files
	core::Real cst_weight_;

	// score functions
	core::scoring::ScoreFunctionOP scorefxn_;
    core::scoring::ScoreFunctionOP highres_scorefxn_;


	// external objects
	Ab_InfoOP ab_info_;

	//packer task
	core::pack::task::TaskFactoryOP tf_;

	// movers
	protocols::antibody2::Ab_H3_Model_CDR_H3OP model_cdrh3_;

	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();

	/// @brief Performs the portion of setup of non-primitive members that requires a pose - called on apply
	void finalize_setup( core::pose::Pose & pose );

	/// @brief Sets up the instance of AntibodyModeler and initializes all members based on values passed in at construction
	///		or via the command line.
	void init();

	void setup_objects();

}; // class Ab_H3_Modeler_full_protocol

    
    
    
    
} // namespace antibody2
} // namespace protocols

#endif

