// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody2/Ab_H3_Model_CDR_H3.hh
/// @brief
///// @author Jianqing Xu ( xubest@gmail.com )
//


#ifndef INCLUDED_protocols_antibody2_Ab_H3_Model_CDR_H3_hh
#define INCLUDED_protocols_antibody2_Ab_H3_Model_CDR_H3_hh

#include <protocols/antibody2/Ab_H3_Model_CDR_H3.fwd.hh>

#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/moves/Mover.hh>
#include <utility/vector1.fwd.hh>
#include <protocols/antibody2/Ab_Info.fwd.hh>
#include <protocols/antibody2/Ab_H3_perturb_ccd_build.fwd.hh>
#include <protocols/antibody2/Ab_H3_cter_insert_mover.fwd.hh>
#include <protocols/moves/PyMolMover.fwd.hh>




using namespace core;

namespace protocols {
namespace antibody2 {

//////////////////////////////////////////////////////////////////////////
/// @brief Ab initio modeling of CDR H3 loop
/// @details
class Ab_H3_Model_CDR_H3 : public protocols::moves::Mover {
public:
	/// @brief default constructor
	Ab_H3_Model_CDR_H3();

	/// @brief constructor with arguments
	Ab_H3_Model_CDR_H3( bool camelid, bool benchmark, Ab_InfoOP antibody_info );

	/// @brief default destructor
	~Ab_H3_Model_CDR_H3();
    
	void set_default();
	virtual void apply( core::pose::Pose & pose_in );
	virtual std::string get_name() const;

	/// @brief enable benchmark mode
	inline void enable_benchmark_mode( bool setting ) {
		benchmark_ = setting;
	}

	/// @brief enable camelid modeling mode
	inline void set_camelid( bool setting ) {
		is_camelid_ = setting;
		if( is_camelid_ ) {
			snug_fit_ = false;
		}
	}


	/// @brief set scorefunction for low resolution of CDR H3 modeling
	void set_lowres_score_func(core::scoring::ScoreFunctionOP lowres_scorefxn );

	/// @brief set scorefunction for high resolution of CDR H3 modeling
	void set_highres_score_func(core::scoring::ScoreFunctionOP highres_scorefxn);

    
	void loop_centroid_relax(
		core::pose::Pose & pose_in,
		core::Size const loop_begin,
		core::Size const loop_end );

    void set_task_factory(pack::task::TaskFactoryCOP tf){
        tf_ = new core::pack::task::TaskFactory(*tf);
    }
    
    void turn_off_H3_filter();
    
    void turn_off_cter_insert(){
        do_cter_insert_ = false;
    }
    
    void turn_on_and_pass_the_pymol(moves::PyMolMoverOP pymol);
    
private:
	bool user_defined_;
	core::pose::Pose start_pose_;
    
    bool do_cter_insert_;

	// constraints
	core::Real cen_cst_;
	core::Real high_cst_;
    core::pose::Pose hfr_pose_;
    
    /// @brief size of loop above which 9mer frags are used
	core::Size cutoff_9_; // default 16
    
	/// @brief size of loop above which 3mer frags are used
	core::Size cutoff_3_; // default 6
    
    /// @brief Number of ADDITIONAL residues modeled from H3_CTERM
	///        These residues range from H:n-2,n-1,n,n+1 of H3
	core::Size c_ter_stem_;
    
    core::Size max_cycle_;
    
    bool use_pymol_diy_;
    moves::PyMolMoverOP pymol_;

	// score functions
	core::scoring::ScoreFunctionOP lowres_scorefxn_;
	core::scoring::ScoreFunctionOP highres_scorefxn_;

	/// @brief benchmark flag
	bool benchmark_;
    
	/// @brief flag indicating that current loop being modeled is CDR H3
	bool current_loop_is_H3_;


	/// @brief refine H3 only
	bool antibody_refine_;

	/// @brief enable docking local refine of LH chains & simultaneous H3 min
	bool snug_fit_;
    
	/// @brief loop_building in docking
	bool loops_flag_;
    
	/// @brief insert fragment in docking
	bool dle_flag_;
    
	/// @brief just refine input loop
	bool refine_input_loop_;

	/// @brief is camelid antibody without light chain
	bool is_camelid_;

	antibody2::Ab_InfoOP ab_info_;

	//packer task
	core::pack::task::TaskFactoryOP tf_;

	void init(bool camelid, bool benchmark, Ab_InfoOP antibody_info);

    void setup_objects();

    Ab_H3_cter_insert_moverOP h3_cter_insert_mover_;
    Ab_H3_perturb_ccd_buildOP h3_perturb_ccd_build_;
    
}; // class Ab_H3_Model_CDR_H3


    
/*
// JQX: make a class, which inherits from abstract class "loop_mover"
class my_LoopMover: public protocols::loops::loop_mover::LoopMover {
public:
    my_LoopMover();
    ~my_LoopMover();
    virtual void apply( core::pose::Pose & pose_in );
	virtual std::string get_name() const;
    virtual void set_extended_torsions( core::pose::Pose & pose, loops::Loop const & loop );
protected:
    virtual basic::Tracer & tr() const;
private:
    
};
*/  
    



} // antibody2
} // protocols


#endif
