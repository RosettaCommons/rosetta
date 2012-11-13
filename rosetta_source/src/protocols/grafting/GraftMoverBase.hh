// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/grafting/GraftMoverBase.hh
/// @brief   Base class for graftmovers
/// @author  Jared Adolf-Bryfogle

#ifndef INCLUDED_protocols_grafting_GraftMoverBase_HH
#define INCLUDED_protocols_grafting_GraftMoverBase_HH

// Unit header
#include <protocols/grafting/GraftMoverBase.fwd.hh>
#include <protocols/moves/Mover.hh>


// Project headers
#include <core/scoring/ScoreFunction.hh>
#include <protocols/loops/Loop.hh>
#include <core/pose/Pose.hh>



namespace protocols {
namespace grafting {
	using namespace core;
	using namespace core::pose;
	using protocols::loops::Loop;
	using core::kinematics::MoveMapOP;
///@brief Base class for GraftMovers.  Created for simplicity and control for C++ programmer, as well as PyRosetta user.
///
/// Feel free to add.
///
class GraftMoverBase: public moves::Mover {

public:
	///@brief Start and end are the residue numbers you want your insert to go between.  start->Insert<-end
	GraftMoverBase(Size const start, Size const end, std::string mover_name);

	virtual ~GraftMoverBase();

	//@brief copy ctor
	//GraftMoverBase( GraftMoverBase const & rhs);

    ///@brief Sets the piece that will be inserted, and any overhang residues.
    ///@details Overhang will be deleted upon insertion.  They are used for the base class function superimpose_overhangs.  Useful if using Double loop methods which keep the insertion frozen in cartesian space.
	virtual void 
    set_piece(Pose & piece, Size Nter_overhang, Size Cter_overhang);

	///@brief Set the region if changed since construction.
	virtual void 
    set_insert_region(Size const start, Size const end);

    virtual void
    set_cen_scorefunction(core::scoring::ScoreFunctionOP score);
    
	virtual void 
    set_fa_scorefunction(core::scoring::ScoreFunctionOP score);
    
	///@brief Grabs the scorefunction from command line options.
	virtual void 
    set_default_fa_scorefunction();
	
    ///@brief Grabs the cen scorefunction set in command line loops option group.
    virtual void
    set_default_cen_scorefunction();

	///@brief uses rms_util to superimpose piece onto pose.  Not run during apply for added control to user.
	void 
    superimpose_overhangs_heavy(Pose const & pose, bool ca_only, bool silence_rms);

	/// @brief  Return the name of the Mover.
	virtual std::string 
    get_name() const;

    




protected:
    
	/// @brief Steven Lewis' insertion method from insert_pose_into_pose. Wrapper to his function, using variables defined in this baseclass.
	/// @details Need to set piece to use. Deletes any overhang in piece.  Deletes any region from start to end in pose. Updates end_.
	/// Recommended use is within apply method.
	///
	Pose 
    insert_piece(Pose const & pose);

	/// @brief Uses a small mover at high KT to perturb residues in the movemap for testing.
	virtual core::Real 
    perturb_backbone_for_test(Pose & pose, MoveMapOP mm);
        
	///@brief combines the two main movemaps to use after the insertion.
	///@details Pose piece must be set.
	MoveMapOP 
    combine_movemaps(MoveMapOP & scaffold_mm, MoveMapOP & insert_mm);
	
	/// @brief deletes overhang residues of the pose piece set.
	/// Recommended use is within apply method
	void 
    delete_overhang_residues();
    
	///@brief Set overhang residues
	void 
    set_overhang(Size Nter_overhang, Size Cter_overhang);
    
    //Reference of the pose piece.  Should be changed to local copy, but I'm not sure how to do that.
	Pose piece_;

	///@brief Residue insertion will start from
	Size start_;
	///@brief Residue insertion will end before here. Updates after insertion.
	Size end_;
	///@brief some functions need to only work on the original numbers. (combine movemaps)
	Size original_end_;

	Size insertion_length_;
	core::scoring::ScoreFunctionOP cen_scorefxn_;
    core::scoring::ScoreFunctionOP fa_scorefxn_;
	///@brief Number of overhang residues on N terminus.  Updates on delete_overhang_residues
	Size Nter_overhang_;
	///@brief Number of overhang residues on C terminus.  Updates on delete_overhang_residues
	Size Cter_overhang_;

protected:
	//////////////////////////////////////////////////////////////////////////////////////////////////////
    ///MOVEMAP and REGION SETUP
    ///Note: Only these functions interact with class variables defined here. To be used optionally in apply.
    //////////////////////////////////////////////////////////////////////////////////////////////////////

    ///@brief sets up either the default movemap or a new combined movemap at apply time.  Updates regions as needed.
	virtual void 
    setup_movemap_and_regions(Pose & pose);
	
	///@brief Sets up the default movemap
	virtual void 
    set_default_movemap();

    ///@brief Sets up the regions at apply using insert variables and flexibility.
	virtual void 
    set_regions_from_flexibility();
    
	///@brief Sets up region variables from the class movemap for the combined pose.  
	virtual void 
    set_regions_from_movemap(Pose & pose);

    MoveMapOP movemap_;
	MoveMapOP scaffold_movemap_;
	MoveMapOP insert_movemap_;
    
    bool use_default_movemap_; //Instructs setup_movemap_and_regions what to do.
	Size Nter_scaffold_flexibility_;
	Size Nter_insert_flexibility_;
	Size Nter_loop_start_;//First flexible residue
	Size Nter_loop_end_;
	
	Size Cter_scaffold_flexibility_;
	Size Cter_insert_flexibility_;
	Size Cter_loop_start_;
	Size Cter_loop_end_; //Last flexible residue

protected:
    /////////////////////////////////////////////////////////////////////////////////////////////////
	///FOLDTREE SETUP.  options depending on how you want your graft algorithm to work!
	///---Indicates Flexible regions, | indicates cutpoint. Arrows are direction of ARMs used to close the loop in conjunction with algorithm. (CCD, KIC)
	///
	/////////////////////////////////////////////////////////////////////////////////////////////////
    

	/////////////////////////////////////////////////////////////////
	/// @brief ****Nter_loop_start---->Piece----> | Cter_loop_end****
	/// Insert will move in cartesian space
	/// @params lower_cutpoint for CCD and loops is Cter_loop_end-1
	///
	virtual void 
    setup_single_loop_single_arm_remodeling_foldtree(Pose & pose, Size const Nter_loop_start, Size const Cter_loop_end, bool loop_modeling=false);

	//////////////////////////////////////////////////////////////////
	/// @brief ****Nter_loop_start---->Piece | <----Nter_loop_end****
	/// Insert will move in cartesian space
	/// @params lower_cutpoint for CCD and loops is end_-1
	///
	virtual void 
    setup_single_loop_double_arm_remodeling_foldtree(Pose & pose, Size const Nter_loop_start, Size const Cter_loop_end, bool loop_modeling=false);

	/////////////////////////////////////////////////////////////////////
	/// ****Nter_loop_start-----> | Piece | <----Nter_loop_end****
	/// or
	///   ****Nter_loop_start----> | <---- Piece -----> | <------ Cter_loop_end****
	/// Insert is fixed
	/// Two loops.  Cutpoint vs loop def determines if one arm or two arms on each side.
	///Use FoldTreeFromLoopsWrapper




};  // class GraftMoverBase

}  // namespace grafting
}  // namespace protocols

#endif  // INCLUDED_protocols_grafting_GraftMoverBase_HH
