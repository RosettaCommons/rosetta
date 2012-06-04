// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer, email:license@u.washington.edu

/// @file protocols/antibody2/RefineCDRH3HighRes.hh
/// @brief Build a homology model of an antibody2
/// @detailed
///
///
/// @author Jianqing Xu (xubest@gmail.com)



#ifndef INCLUDED_protocols_antibody2_RefineCDRH3HighRes_hh
#define INCLUDED_protocols_antibody2_RefineCDRH3HighRes_hh


#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/antibody2/AntibodyInfo.fwd.hh>
#include <protocols/antibody2/RefineCDRH3HighRes.fwd.hh>



using namespace core;
namespace protocols {
namespace antibody2 {

class RefineCDRH3HighRes: public moves::Mover {


public:

    /// @brief default constructor
	RefineCDRH3HighRes();
    
    /// @brief constructor with arguments
    RefineCDRH3HighRes( AntibodyInfoOP antibody_info);

    /// @brief constructor with arguments
    RefineCDRH3HighRes(AntibodyInfoOP antibody_info, std::string refine_mode);
    
    /// @brief constructor with arguments
    RefineCDRH3HighRes(AntibodyInfoOP antibody_info, std::string refine_mode, 
                       core::scoring::ScoreFunctionCOP highres_scorefxn );


	/// @brief default destructor
	~RefineCDRH3HighRes();

    
    void turn_on_benchmark(){
        benchmark_=true;
    }
    void set_h3_filter(bool setting){
        H3_filter_=setting;
    }
    void set_num_filter_tries(core::Size setting){
        num_filter_tries_=setting;
    }
    void set_flank_relax(bool setting){
        flank_relax_=setting;
    }
    void set_flank_size(core::Size setting){
        flank_size_=setting;
    }
    void set_refine_mode(std::string refine_mode){
        refine_mode_ = refine_mode;
    }
    virtual void apply( core::pose::Pose & pose );

    
	void set_highres_score_func(core::scoring::ScoreFunctionCOP highres_scorefxn){
        highres_scorefxn_ = new core::scoring::ScoreFunction(*highres_scorefxn);
    }

    void pass_start_pose(core::pose::Pose & start_pose);
    virtual std::string get_name() const;
    virtual protocols::moves::MoverOP clone() const;

    
private:

    AntibodyInfoOP ab_info_;
    bool user_defined_;
    bool benchmark_;
	core::scoring::ScoreFunctionOP highres_scorefxn_;
	core::Size flank_size_;
	bool flank_relax_;
    core::pose::Pose start_pose_;
	core::pack::task::TaskFactoryOP tf_;
    bool H3_filter_;
    core::Size num_filter_tries_;
    std::string refine_mode_;
    core::Real high_cst_;
    
    void set_default();
    void init();
    void finalize_setup( core::pose::Pose & pose );
    
};







} // namespace antibody2
} // namespace protocols

#endif


