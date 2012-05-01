// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody2/ModelCDRH3.hh
/// @brief
///// @author Jianqing Xu ( xubest@gmail.com )
//


#ifndef INCLUDED_protocols_antibody2_ModelCDRH3_hh
#define INCLUDED_protocols_antibody2_ModelCDRH3_hh

#include <protocols/antibody2/ModelCDRH3.fwd.hh>

#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/moves/Mover.hh>
#include <utility/vector1.fwd.hh>
#include <protocols/antibody2/AntibodyInfo.fwd.hh>
#include <protocols/antibody2/H3PerturbCCD.fwd.hh>
#include <protocols/antibody2/H3CterInsert.fwd.hh>
#include <protocols/moves/PyMolMover.fwd.hh>
#include <basic/Tracer.hh>





using namespace core;

namespace protocols {
namespace antibody2 {

//////////////////////////////////////////////////////////////////////////
/// @brief Ab initio modeling of CDR H3 loop
/// @details
class ModelCDRH3 : public protocols::moves::Mover {
public:
	/// @brief default constructor
	ModelCDRH3();

    /// @brief constructor with arguments
	ModelCDRH3(  AntibodyInfoOP antibody_info );
    
	/// @brief constructor with arguments
	ModelCDRH3(  AntibodyInfoOP antibody_info ,                         
                core::scoring::ScoreFunctionCOP lowres_scorefxn,
                core::scoring::ScoreFunctionCOP highres_scorefxn);

	/// @brief default destructor
	~ModelCDRH3();
    
	void set_default();
	virtual void apply( pose::Pose & pose_in );
	virtual std::string get_name() const;

	/// @brief enable benchmark mode
	inline void enable_benchmark_mode( bool setting ) {
		benchmark_ = setting;
	}

	/// @brief enable camelid modeling mode
	inline void set_camelid( bool setting ) {
		is_camelid_ = setting;
	}


	/// @brief set scorefunction for low resolution of CDR H3 modeling
	void set_lowres_score_func(scoring::ScoreFunctionCOP lowres_scorefxn );

	/// @brief set scorefunction for high resolution of CDR H3 modeling
	void set_highres_score_func(scoring::ScoreFunctionCOP highres_scorefxn);

	/// @brief set task factory 
    void set_task_factory(pack::task::TaskFactoryCOP tf);
    
    void turn_off_H3_filter();
    
    void turn_off_cter_insert(){
        do_cter_insert_ = false;
    }
    
    void turn_on_and_pass_the_pymol(moves::PyMolMoverOP pymol);
    
private:
	bool user_defined_;
	pose::Pose start_pose_;
    
    bool do_cter_insert_;

	// constraints
	Real cen_cst_;
	Real high_cst_;
    pose::Pose hfr_pose_;
    
    /// @brief size of loop above which 9mer frags are used
	Size cutoff_9_; // default 16
    
	/// @brief size of loop above which 3mer frags are used
	Size cutoff_3_; // default 6
    
    /// @brief Number of ADDITIONAL residues modeled from H3_CTERM
	///        These residues range from H:n-2,n-1,n,n+1 of H3
	Size c_ter_stem_;
    
    Size max_cycle_;
    
    bool use_pymol_diy_;
    moves::PyMolMoverOP pymol_;

	// score functions
	scoring::ScoreFunctionOP lowres_scorefxn_;
	scoring::ScoreFunctionOP highres_scorefxn_;

	/// @brief benchmark flag
	bool benchmark_;
    
	/// @brief flag indicating that current loop being modeled is CDR H3
	bool current_loop_is_H3_;

    
	/// @brief loop_building in docking
	bool loops_flag_;
    
	/// @brief insert fragment in docking
	bool dle_flag_;
    


	/// @brief is camelid antibody without light chain
	bool is_camelid_;

	antibody2::AntibodyInfoOP ab_info_;

	//packer task
	pack::task::TaskFactoryOP tf_;

	void init();

    H3CterInsertOP h3_cter_insert_mover_;
    H3PerturbCCDOP h3_perturb_ccd_build_;
    
}; // class ModelCDRH3


    


// JQX: make a class, which inherits from abstract class "loop_mover"
class my_LoopMover: public protocols::loops::loop_mover::LoopMover {
protected:
    virtual basic::Tracer & tr() const;
};
    



} // antibody2
} // protocols


#endif
