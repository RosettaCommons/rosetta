// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/LHRepulsiveRamp.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#ifndef INCLUDED_protocols_antibody2_LHRepulsiveRamp_hh
#define INCLUDED_protocols_antibody2_LHRepulsiveRamp_hh





#include <protocols/antibody2/LHRepulsiveRamp.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/antibody2/AntibodyInfo.hh>
#include <core/scoring/ScoreFunction.fwd.hh>




using namespace core;
namespace protocols {
namespace antibody2 {
        
class LHRepulsiveRamp: public moves::Mover {
            
            
public:
    
    /// @brief default constructor
	LHRepulsiveRamp();
    
	/// @brief constructor with arguments
    LHRepulsiveRamp(loops::Loops loops_in );
    LHRepulsiveRamp(antibody2::AntibodyInfoOP antibody_in );
	LHRepulsiveRamp(antibody2::AntibodyInfoOP antibody_in, bool camelid );
        
    virtual protocols::moves::MoverOP clone() const;
    
	/// @brief default destructor
	~LHRepulsiveRamp();
    
    void set_default();
    

    virtual void apply( core::pose::Pose & pose );
    
    virtual std::string get_name() const;
    
    
    void set_task_factory(pack::task::TaskFactoryCOP tf){
        tf_ = new pack::task::TaskFactory(*tf);
    }
    
private:

    AntibodyInfoOP ab_info_;
    
    bool user_defined_;
    bool benchmark_;
    bool is_camelid_;
    loops::Loops all_loops_; 
    Size nres_;
    kinematics::MoveMapOP cdr_dock_map_;
    Size rep_ramp_cycles_;
    std::string min_type_;
    Real rot_mag_;
    Real trans_mag_;
    Real temperature_;
    Real min_threshold_;
    Size cycles_;

    
    scoring::ScoreFunctionOP dock_scorefxn_;
    scoring::ScoreFunctionOP pack_scorefxn_;
    
    void init(loops::Loops loops_in, bool camelid);
    
    void setup_objects();
    
    void finalize_setup(pose::Pose & pose );


    
	void repulsive_ramp( core::pose::Pose & pose_in, loops::Loops loops_in );
    
	void snugfit_MC_min (core::pose::Pose & pose_in);
    
    
	//packer task
    pack::task::TaskFactoryOP tf_;

};
    
    
    
    
    
    
    
} // namespace antibody2
} // namespace protocols

#endif








