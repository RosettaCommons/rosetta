// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/Ab_ModelCDRH3.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#ifndef INCLUDED_protocols_antibody2_Ab_ModelCDRH3_hh
#define INCLUDED_protocols_antibody2_Ab_ModelCDRH3_hh

#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
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

	virtual std::string get_name() const;
    
    /// @brief Associates relevant options with the AntibodyModeler class
    static void register_options();
    
	// simple inline setters
    void set_h3modeler( bool model_h3 ){ model_h3_ = model_h3; }
	void set_snugfit( bool snugfit ) { snugfit_ = snugfit; }
	void set_camelid( bool camelid ) { camelid_ = camelid; }
	void set_camelid_constraints( bool camelid_constraints ) 
        { camelid_constraints_ = camelid_constraints; }
	void set_benchmark( bool benchmark ) { benchmark_ = benchmark; }
    void set_cst_weight( core::Real const cst_weight){cst_weight_=cst_weight;}


	void relax_cdrs( core::pose::Pose & pose );

	void all_cdr_VL_VH_fold_tree( core::pose::Pose & pose_in, 
                                  const loops::Loops & loops );


	core::Real global_loop_rmsd ( const core::pose::Pose & pose_in, 
                                  const core::pose::Pose & native_pose, 
                                  std::string cdr_type );

	void display_constraint_residues( core::pose::Pose & pose );

        
    void show( std::ostream & out=std::cout );
    friend std::ostream & operator<<(std::ostream& out, const Ab_ModelCDRH3 & ab_m_2 );
    
    

private:
    bool model_h3_;
	bool snugfit_;
	bool camelid_;
	bool camelid_constraints_;
    core::pose::Pose start_pose_;
    /// @brief refine H3 only
	bool antibody_refine_;
    bool is_camelid_;
    core::scoring::ScoreFunctionOP highres_scorefxn_;
    core::Real high_cst_;

	// Benchmark mode for shorter_cycles
	bool benchmark_;

	bool user_defined_; // for constructor options passed to init

	bool flags_and_objects_are_in_sync_;
	bool first_apply_with_current_setup_;

	// used as a flag to enable reading in of cst files
	core::Real cst_weight_;

	// score functions
	core::scoring::ScoreFunctionOP scorefxn_;

	// external objects
	Ab_InfoOP ab_info_;


	//packer task
	core::pack::task::TaskFactoryOP tf_;

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
