// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/AntibodyModelerProtocol.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu ( xubest@gmail.com )


#ifndef INCLUDED_protocols_antibody2_AntibodyModelerProtocol_hh
#define INCLUDED_protocols_antibody2_AntibodyModelerProtocol_hh

#include <utility/vector1.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/moves/PyMolMover.fwd.hh>
#include <protocols/antibody2/AntibodyInfo.hh>
#include <protocols/antibody2/AntibodyModelerProtocol.fwd.hh>




using namespace core;
namespace protocols {
namespace antibody2 {

class AntibodyModelerProtocol: public moves::Mover {
public:

	// default constructor
	AntibodyModelerProtocol();

	// default destructor
	~AntibodyModelerProtocol();

	virtual protocols::moves::MoverOP clone() const;

	/// @brief Assigns default values to primitive members
	void set_default();

	/// @brief Instantiates non-primitive members based on the value of the primitive members
	void sync_objects_with_flags();

	virtual void apply( pose::Pose & pose );

	virtual std::string get_name() const;
    
    /// @brief Associates relevant options with the AntibodyModeler class
    static void register_options();
    
	// simple inline setters
    void set_BenchMark  ( bool benchmark  ) {benchmark_   = benchmark;}
    void set_ModelH3    ( bool model_h3   ) {model_h3_    = model_h3; }
	void set_SnugFit    ( bool snugfit    ) {snugfit_     = snugfit;  }
    void set_H3Filter   ( bool H3_filter  ) {H3_filter_   = H3_filter;}
    void set_CterInsert ( bool cter_insert) {cter_insert_ = cter_insert;}
    void set_cst_weight ( core::Real const cst_weight){cst_weight_=cst_weight;}
    void set_camelid( bool camelid )  { camelid_ = camelid; }
	void set_camelid_constraints(bool camelid_constraints ) { 
        camelid_constraints_ = camelid_constraints; 
    }
    void set_sc_min (bool scmin) {sc_min_ = scmin ;}
    void set_rt_min (bool rtmin) {rt_min_ = rtmin ;}

	void display_constraint_residues( pose::Pose & pose );
        
    void show( std::ostream & out=std::cout );
    friend std::ostream & operator<<(std::ostream& out, const AntibodyModelerProtocol & ab_m );
    
    

private:
    bool model_h3_;
	 bool snugfit_;
    bool refine_h3_;
    bool H3_filter_;
    bool cter_insert_;
    bool LH_repulsive_ramp_;
    bool sc_min_;
    bool rt_min_;
    bool camelid_;
	 bool camelid_constraints_;
    
    /// @brief refine H3 only
    core::Real cen_cst_, high_cst_;
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
    core::scoring::ScoreFunctionOP loop_scorefxn_highres_;
    core::scoring::ScoreFunctionOP loop_scorefxn_centroid_;
    core::scoring::ScoreFunctionOP dock_scorefxn_highres_;
    core::scoring::ScoreFunctionOP pack_scorefxn_;
    
	// external objects
	AntibodyInfoOP ab_info_;

	//packer task
	pack::task::TaskFactoryOP tf_;

	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();

	/// @brief Performs the portion of setup of non-primitive members that requires a pose - called on apply
	void finalize_setup( core::pose::Pose & pose );

	/// @brief Sets up the instance of AntibodyModeler and initializes all members based on values passed in at construction
	///		or via the command line.
	void init();

	void setup_objects();

}; // class
    
    
    
    
} // namespace antibody2
} // namespace protocols

#endif

