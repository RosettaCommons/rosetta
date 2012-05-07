// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/RefineBetaBarrel.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)


#ifndef INCLUDED_protocols_antibody2_RefineBetaBarrel_hh
#define INCLUDED_protocols_antibody2_RefineBetaBarrel_hh

#include <protocols/moves/Mover.hh>
#include <protocols/antibody2/RefineBetaBarrel.fwd.hh>
#include <protocols/antibody2/AntibodyInfo.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/antibody2/LHRepulsiveRamp.fwd.hh>
#include <protocols/antibody2/LHSnugFitLegacy.fwd.hh>
#include <protocols/docking/DockMCMProtocol.fwd.hh>
#include <protocols/moves/PyMolMover.fwd.hh>

#ifdef PYROSETTA
#include <protocols/moves/PyMolMover.hh>
#endif


using namespace core;
namespace protocols {
namespace antibody2 {
        
        
        
class RefineBetaBarrel: public moves::Mover {
            
            
public:
    /// @brief default constructor
    RefineBetaBarrel();
        
    /// @brief default destructor
    ~RefineBetaBarrel();
        
    RefineBetaBarrel(AntibodyInfoOP antibody_info);
    
    RefineBetaBarrel(AntibodyInfoOP antibody_info,
                    core::scoring::ScoreFunctionCOP dock_scorefxn,
                    core::scoring::ScoreFunctionCOP pack_scorefxn);
    
    virtual void apply( core::pose::Pose & pose_in );
    virtual std::string get_name() const;
    
    void set_task_factory(core::pack::task::TaskFactoryCOP tf);
    
    
    void set_dock_score_func(core::scoring::ScoreFunctionCOP dock_scorefxn ){
        dock_scorefxn_ = new core::scoring::ScoreFunction(*dock_scorefxn);
    }
    
    void set_pack_score_func(core::scoring::ScoreFunctionCOP pack_scorefxn){
        pack_scorefxn_ = new core::scoring::ScoreFunction(*pack_scorefxn);
    }
    
    void turn_off_repulsive_ramp(){
        repulsive_ramp_ = false;
    }
    
    void turn_on_and_pass_the_pymol(moves::PyMolMoverOP pymol){
        use_pymol_diy_ = true;
        pymol_ = pymol;
    }
    
    void set_sc_min(bool sc_min){
        sc_min_ = sc_min;
    }
    
    void set_rt_min(bool rt_min){
        rt_min_ = rt_min;
    }
            
private:
    bool sc_min_;
    bool rt_min_;
    bool user_defined_;
    bool repulsive_ramp_;
    AntibodyInfoOP ab_info_;
    pack::task::TaskFactoryOP tf_;
    
    moves::PyMolMoverOP pymol_;
    bool use_pymol_diy_;
        
    void init( );
    void finalize_setup(core::pose::Pose & pose_in );
    
    LHRepulsiveRampOP lh_repulsive_ramp_;
    LHSnugFitLegacyOP lh_snugfit_;
    docking::DockMCMProtocolOP dock_mcm_protocol_;
    
    core::scoring::ScoreFunctionOP dock_scorefxn_;
    core::scoring::ScoreFunctionOP pack_scorefxn_;
    
    kinematics::MoveMapOP cdr_dock_map_;
    
    
};
    
    
}// namespace antibody2
}// namespace protocols




#endif


