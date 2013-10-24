// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/grafting/AnchoredGraftMover.hh
/// @brief   Class to graft a piece into a pose.
/// @Author  Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author  Original algorithm - Steven Lewis (smlewi@gmail.com)

#ifndef INCLUDED_protocols_grafting_AnchoredGraftMover_HH
#define INCLUDED_protocols_grafting_AnchoredGraftMover_HH


//Unit Headers
#include <protocols/grafting/GraftMoverBase.hh>
#include <protocols/grafting/AnchoredGraftMover.fwd.hh>
#include <protocols/moves/Mover.hh>
//Core
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>




namespace protocols {
namespace grafting {
	using namespace core::pose;

///@brief Grafting class adapted from Steven Lewis' pose_into_pose algorithm.  Basic, and quick, but with many options.
///
/// example:
///    mover = AnchoredGraftMover(start, end)
///    mover.set_piece(piece, cter_overhang, nter_overhang)
///    mover.apply(pose)
///
/// see also: grafting/util.hh.
///
///@details Uses a single loop and a single arm to close the loop by default.
/// ****Nter_loop_start---->Piece----> | Cter_loop_end****
/// Default movemap keeps insert Frozen in dihedral angle space, But think of the insert as part of a giant arm.
/// Default flexibility on Nter and Cter is only two residues (--> part of diagram).
/// Will delete any residues between start and end, and any overhang residues from the insert.
///
/// Algorithm originally from pose_into_pose:
/// The insert will be left unchanged in internal-coordinate space except for the phi on the first residue, and the psi/omega on the last residue, and atoms whose bonding partners change as a result of the insertion.
/// Internally, apply performs the insertion, idealizes the loop residues (omegas to 180, peptide bonds idealized) and the newly made polymer connections at the insert point, and then attempts to close the loop(s).
/// It is intended, but not guaranteed, to produce a loop with good rama, omega, and chainbreak/peptide_bond scores.  It does NOT attempt to give a loop with good sidechains (it does not repack at all) or worry overmuch about van der Waals
///
class AnchoredGraftMover : public protocols::grafting::GraftMoverBase {
public:
    
    ///@brief Start and end are the residue numbers you want your insert to go between.  start->Insert<-end
    AnchoredGraftMover(Size const start, Size const end);

    virtual ~AnchoredGraftMover();
    
    void 
    set_defaults();
    
	///@brief Sets scaffold flexiblity on either end of scaffold
	void 
    set_scaffold_flexibility(Size const Nter_scaffold_flexibility, Size const Cter_scaffold_flexibility);

	Size
	get_nterm_scaffold_flexibility();
	
	Size
	get_cterm_scaffold_flexibility();
	
	///@brief Sets insert flexibility on either end of insert
	void 
    set_insert_flexibility(Size const Nter_insert_flexibility, Size const Cter_insert_flexibility);
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///@brief Advanced way to set options.  Use if you know what you are doing.
	///@details Will combine the movemaps for apply, and renumber everything. Flexible residues in multiple chains not recommended.
	/// single_loop_single_arm (Default) From the first bb flexible residue to the last will act as one loop.  This can be amazing as you can use loop regions in various parts of your protein to help the insertion.
	/// single_loop_double_arm: One arm will be from first flexible N terminal residue to end_; Second arm will be from there to last flexible residue
	/// double_loop_double_arm: One arm Will go from first flexible N terminal residue to after start_+any contiguous residues on in the movemap from there. Opposite for Cter side.
	/// double_loop_quad_arm:   same as double, but cutpoints will be at each new peptide bond.
	///
	/// Note: Will disregard flexibility settings, as the movemaps will be used as primary way to define flexibility.
	void 
    set_movemaps(MoveMapOP const scaffold_mm, MoveMapOP const insert_mm);

	///@brief Instructs the class to ignore any movemaps from set_movemaps.
	void 
    set_use_default_movemap_from_flexibility(bool def);
	
	void 
    set_cycles(Size cycles);
    
    ///@brief override of base class
    void 
    set_default_cen_scorefunction();

    ///@brif override of base class
    void
    set_cen_scorefunction(core::scoring::ScoreFunctionOP scorefxn_);
    
    ///@brief Sets the mintype for the MinMover
    void 
    set_mintype(std::string mintype);

    ///@brief Sets the mover to skip the small mover sampling step.
    void 
    set_skip_sampling(bool skip_sampling);

    ///@brief sets up the smooth centroid scorefunction + any changes to VDW.  if false, switches back to default scorefunction.
    void 
    set_use_smooth_centroid_settings(bool use_smooth);
    
    ///@brief Uses a single loop, two arm loop closer.
    ///@details ****Nter_loop_start---->Piece | <----Nter_loop_end****
    void 
    set_use_single_loop_double_CCD_arms(bool single_loop_double_arm);

    ///@brief Keeps the insert frozen in cartesian space. Sets the algorithm to use two loops on either side in order to close the graft. (Please superimpose piece onto scaffold first)
    ///@details ****Nter_loop_start-----> | Piece | <----Nter_loop_end****
    void 
    set_use_double_loop_double_CCD_arms(bool double_loop_double_arm);

    ///@brief Keeps the insert frozen in cartesian space. Sets the algorithm to use two loops on either side in order to close the graft. (Please superimpose piece onto scaffold first)
    ///@details ****Nter_loop_start----> | <---- Piece -----> | <------ Cter_loop_end**** Use only if you have continuous flexibility in your insert from start and or end.
    void
    set_use_double_loop_quad_CCD_arms(bool double_loop_quad_arm);

    ///@returns the Cterminal loop end (Last flexible residue).  Useful to use after insertion.
    Size 
    get_Cter_loop_end();

	///@brief TESTING ONLY Sets the protocol to 'randomize' the flexible residues before trying to graft.  This is used to test the protocol by grafting a piece of a protein back onto itself and looking at RMSD.
	void 
    set_test_control_mode(bool test_control_mode);

    ///@brief Grafts the piece into the pose, uses CCD to close the connection.  Insert does not change dihedral space, but DOES change cartesian space by default.
    ///Pose is returned without repacking any sidechains.  Use repack_connection after apply method.
    ///Deletes overhang and region between start and end if residues are present.
    virtual void 
    apply(Pose & pose);

public:
    
    ///@brief convenience function for AFTER apply method.
    ///@details flexible Nter and Cter residues plus the first and last residue of the insert.
    /// If passing movemaps, will respect those movemaps.
    virtual void 
    repack_connection_and_residues_in_movemap(Pose & pose, core::scoring::ScoreFunctionOP fa_scorefxn);

    ///@brief convenience function for AFTER apply method.
    ///@details flexible Nter and Cter residues plus the entire insert.
    /// If passing movemaps, will respect those movemaps.
    virtual void 
    repack_connection_and_residues_in_movemap_and_piece(Pose & pose, core::scoring::ScoreFunctionOP fa_scorefxn);

private:

	Size cycles_;
    
	std::string mintype_;
	bool skip_sampling_;//Option to skip the small mover sampling step.
	bool single_loop_double_arm_; //Option to use two CCD arms to close the loop.
	bool double_loop_double_arm_;//Option to freeze the insert in cartesian space.  Uses two arms, one on either side of insert to do the graft.
	bool double_loop_quad_arm_;
	bool test_control_mode_;//TESTING ONLY  Set to randomize the flexible residues before trying to graft.

}; //Class AnchoredGraftMover


}// namespace grafting
}// namespace protocols

#endif  // INCLUDED_protocols_grafting_AnchoredGraftMover_HH
