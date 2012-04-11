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
#include <core/pack/task/TaskFactory.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/antibody2/LHRepulsiveRamp.fwd.hh>
#include <protocols/antibody2/LHSnugFitLegacy.fwd.hh>
#include <protocols/docking/DockMCMProtocol.fwd.hh>




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
    
        void set_default();
        virtual void apply( core::pose::Pose & pose_in );
        virtual std::string get_name() const;
    
        void set_task_factory(core::pack::task::TaskFactoryCOP tf){
            tf_ = new pack::task::TaskFactory(*tf);
        }
    
        void turn_off_repulsive_ramp(){
            repulsive_ramp_ = false;
        }
            
    private:
        bool user_defined_;
        bool repulsive_ramp_;
        AntibodyInfoOP ab_info_;
    	pack::task::TaskFactoryOP tf_;
        
        void init(AntibodyInfoOP antibody_info);

    
        LHRepulsiveRampOP lh_repulsive_ramp_;
        LHSnugFitLegacyOP lh_snugfit_;
        docking::DockMCMProtocolOP dock_mcm_protocol_;
    
    
};
    
    
}// namespace antibody2
}// namespace protocols




#endif


