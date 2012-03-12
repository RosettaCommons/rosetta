// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/Ab_LH_RepulsiveRamp_Mover.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#ifndef INCLUDED_protocols_antibody2_Ab_LH_RepulsiveRamp_Mover_hh
#define INCLUDED_protocols_antibody2_Ab_LH_RepulsiveRamp_Mover_hh






#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/antibody2/Ab_Info.hh>

#include <protocols/antibody2/Ab_LH_RepulsiveRamp_Mover.fwd.hh>




namespace protocols {
namespace antibody2 {
        
class Ab_LH_RepulsiveRamp_Mover: public moves::Mover {
            
            
public:
    
    /// @brief default constructor
	Ab_LH_RepulsiveRamp_Mover();
    
	/// @brief constructor with arguments
    Ab_LH_RepulsiveRamp_Mover(loops::Loops loops_in );
    Ab_LH_RepulsiveRamp_Mover(antibody2::Ab_Info & antibody_in );
	Ab_LH_RepulsiveRamp_Mover(antibody2::Ab_Info & antibody_in, bool camelid );
        
    virtual protocols::moves::MoverOP clone() const;
    
	/// @brief default destructor
	~Ab_LH_RepulsiveRamp_Mover();
    
    void set_default();
    

    virtual void apply( core::pose::Pose & pose );
    
    virtual std::string get_name() const;
    
    
private:

    Ab_Info ab_info_;
    
    bool user_defined_;
    bool benchmark_;
    bool is_camelid_;
    loops::Loops all_loops_; 
    
    
    
    void init(loops::Loops loops_in, bool camelid);
    
    void setup_objects();
    
    
    
    void setup_packer_task(core::pose::Pose & Pose_in);
    
	void repulsive_ramp( core::pose::Pose & pose_in, loops::Loops loops_in );
    
	void snugfit_MC_min (
                         core::pose::Pose & pose_in,
                         core::kinematics::MoveMapOP cdr_dock_map,
                         core::Size cycles,
                         core::Real minimization_threshold,
                         core::scoring::ScoreFunctionOP scorefxn,
                         core::scoring::ScoreFunctionOP pack_scorefxn,
                         utility::vector1< bool> is_flexible );
    
    
	//packer task
	core::pack::task::TaskFactoryOP tf_;
	core::pack::task::TaskFactoryOP init_task_factory_;

};
    
    
    
    
    
    
    
} // namespace antibody2
} // namespace protocols

#endif








