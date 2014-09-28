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


//Protocols
#include <protocols/simple_moves/MinMover.fwd.hh>
#include <protocols/simple_moves/BackboneMover.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <utility/tag/Tag.hh>

namespace protocols {
namespace grafting {
	using core::pose::Pose;
	using core::kinematics::MoveMapCOP;
	using core::kinematics::MoveMapOP;
	using core::scoring::ScoreFunctionCOP;
	using core::scoring::ScoreFunctionOP;
	using protocols::simple_moves::MinMoverOP;
	using protocols::simple_moves::SmallMoverOP;
	
///@brief Grafting class adapted from Steven Lewis' pose_into_pose algorithm.
///
/// 1) Inserts the pose piece into the scaffold pose, deleting any overhang residues or residues in the region the insertion will occur between.
/// 2) Connects the left side of the piece to the scaffold by inserting ideal geometry 
///     (aka superposition of terminal overhangs is not needed)
/// 3) Cycles of:
///      a) SmallMover for sampling
///      b) CCD to close right (Nter) side of the graft.  
///        The insert is included as dihedral-inflexible residues in the CCD arm.
///        Residues on either side of the scaffold are used by default to close the loop. 
///        A Moevemap can be given to include more residues for CCD to use, such as loops connecting secondary structural elements
///      c) MinMover to help close the graft through chainbreak terms
///      d) Closure check - will stop if closed (This can be overridden)
///      d) MonteCarlo Boltzmann criterion
///
/// 5) Repack flexible residues, insert residues, connecting residues between scaffold and insert, and neighbors
///
/// example:
///    mover = AnchoredGraftMover(start, end, piece, cter_overhang, nter_overhang);
///    mover.apply(pose);
///
/// see also: grafting/util.hh.
///
///@details Uses a single arm to close the loop by default.
/// ****Nter_loop_start---->Piece----> | Cter_loop_end****
/// Default movemap keeps insert Frozen in dihedral angle space, But think of the insert as part of a CCD giant arm.
/// Default flexibility on Nter and Cter is only two residues (--> part of diagram).
/// Will delete any residues between start and end, and any overhang residues from the insert.
///
/// Algorithm originally from pose_into_pose:
/// The insert will be left unchanged in internal-coordinate space except for the phi on the first residue, and the psi/omega on the last residue, and atoms whose bonding partners change as a result of the insertion.
/// Internally, apply performs the insertion, idealizes the loop residues (omegas to 180, peptide bonds idealized) and the newly made polymer connections at the insert point, and then attempts to close the loop(s).
/// It is intended, but not guaranteed, to produce a graft with good rama, omega, and chainbreak/peptide_bond scores.  All-atom minimization of graft or pose after insertion is recommended.
///
class AnchoredGraftMover : public protocols::grafting::GraftMoverBase {
public:
    
	AnchoredGraftMover();
	
	///@brief Start and end are the residue numbers you want your insert to go between.  start->Insert<-end
	AnchoredGraftMover(Size const start, Size const end, bool copy_pdbinfo = false);

	AnchoredGraftMover(
		Size const start, Size const end, 
		core::pose::Pose const & piece, Size Nter_overhang_length=0, Size Cter_overhang_length=0, bool copy_pdbinfo = false);
    
	AnchoredGraftMover(AnchoredGraftMover const & src);
	
	virtual ~AnchoredGraftMover();

	void 
	set_cycles(Size cycles);
	
	virtual void 
	set_defaults();
	
	///@brief Grafts the piece into the pose, uses CCD to close the connection.  Insert does not change dihedral space, but DOES change cartesian space by default.
	///Does not repack any sidechains.
	///Deletes overhang and region between start and end if residues are present.
	virtual void 
	apply(Pose & pose);
	
public:
	
	protocols::moves::MoverOP
	clone() const;
	
	protocols::moves::MoverOP
	fresh_instance() const;
	
	
	virtual void 
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		moves::Movers_map const & movers,
		Pose const & pose
	);
	
	
public:
	///@brief Stop at closure of both ends? 
	/// Default True.
	void stop_at_closure(bool stop_at_closure);
	
	///@brief Stop at closure of both ends?
	bool stop_at_closure();
	
	///@brief Pack sidechains of flexible residues, insert and neighbors at end?
	/// Default True.
	void final_repack(bool final_repack);
	
	///@brief Pack sidechains of flexible residues, insert and neighbors at end?
	bool final_repack();
	
	void idealize_insert(bool idealize);
	
	bool idealize_insert() const {
		return idealize_insert_;
	}
	
public:
    
	virtual void set_cen_scorefunction(core::scoring::ScoreFunctionCOP score);
    
	virtual void set_fa_scorefunction(core::scoring::ScoreFunctionCOP score);
	
	///@brief Sets the mintype for the MinMover
	void set_mintype(std::string mintype);

	///@brief Sets the mover to skip the small mover sampling step.
	void set_skip_sampling(bool skip_sampling);
	
public:

	///@brief Sets scaffold flexiblity on either end of scaffold
	virtual void set_scaffold_flexibility(Size const Nter_scaffold_flexibility, Size const Cter_scaffold_flexibility);
	
	///@brief Sets insert flexibility on either end of insert
	virtual void set_insert_flexibility(Size const Nter_insert_flexibility, Size const Cter_insert_flexibility);
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///@brief Advanced way to set flexibility and residues used for CCD
	///@details Will combine the movemaps for apply, and renumber everything. Flexible residues in multiple chains not recommended.
	/// From the first bb flexible residue to the last will act as one giant CCD arm, 
	///  using those residues for minimization and CCD that have bb true in movemap. 
	///  This can be amazing as you can use loop regions in various parts of your protein to help the graft complete.
	/// Note: Will disregard flexibility settings, as the movemaps will be used as primary way to define flexibility.
	/// May want to consider turning off the sampling step when passing crazy movemaps.
	
	virtual void set_movemaps(MoveMapCOP const scaffold_mm, MoveMapCOP const insert_mm);

	///@brief Neighbor distance for any repacking of side-chains.
	void neighbor_dis(core::Real dis);
	core::Real neighbor_dis() const;
	
public:
	Size get_nterm_scaffold_flexibility();
	
	Size get_cterm_scaffold_flexibility() ;

	///@returns the Cterminal loop end (Last flexible residue).  Useful to use after insertion.
	Size get_Cter_loop_end();

	protocols::loops::Loops
	get_loops() {
		return *loops_;
	};
	
protected:
	
	//ScoreFunction Setup
	void
	setup_scorefunction();
    
    	virtual void 
	set_default_fa_scorefunction();
	
	///@brief Smooth version of Steven's original scorefunction used for grafting
	virtual void
	set_default_cen_scorefunction();
	
	void
	set_loops(protocols::loops::LoopsOP loops);
	
protected:
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	// MOVEMAP + REGION +SCORE SETUP
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

	
protected:
	//Testing
		
	///@brief TESTING ONLY Sets the protocol to 'randomize' the flexible residues before trying to graft.  This is used to test the protocol by grafting a piece of a protein back onto itself and looking at RMSD.
	void 
	set_test_control_mode(bool test_control_mode);
	
protected:
	virtual MinMoverOP
	setup_default_min_mover();
	
	virtual SmallMoverOP
	setup_default_small_mover();
	
protected:
	//Accessors and Mutators of private data for derived classes
	void movemap(MoveMapOP movemap);
	MoveMapOP movemap() const;
	
	void scaffold_movemap(MoveMapCOP scaffold_movemap);
	MoveMapCOP scaffold_movemap() const;
	
	void insert_movemap(MoveMapCOP insert_movemap);
	MoveMapCOP insert_movemap() const;
	
	ScoreFunctionOP cen_scorefxn() const;
	ScoreFunctionOP fa_scorefxn() const;
	
	void use_default_movemap(bool use_default_movemap);
	bool use_default_movemap() const;
	
	bool test_control_mode() const;
	
	Size Nter_scaffold_flexibility() const;
	Size Cter_scaffold_flexibility() const;
	Size Nter_insert_flexibility() const;
	Size Cter_insert_flexibility() const;
	Size Nter_loop_start() const;
	Size Nter_loop_end() const;
	Size Cter_loop_start() const;
	Size Cter_loop_end() const;
	
	std::string mintype() const;
	Size cycles() const;
	bool skip_sampling() const;
	
	TagCOP tag() const;
	
private:
	Size cycles_;
    
	std::string mintype_;
	bool skip_sampling_;//Option to skip the small mover sampling step.
	bool test_control_mode_;//TESTING ONLY  Set to randomize the flexible residues before trying to graft.

	MoveMapOP movemap_;
	MoveMapCOP scaffold_movemap_;
	MoveMapCOP insert_movemap_;
	
	core::scoring::ScoreFunctionOP cen_scorefxn_;
	core::scoring::ScoreFunctionOP fa_scorefxn_;
	
	bool use_default_movemap_; //Instructs setup_movemap_and_regions what to do.
	Size Nter_scaffold_flexibility_;
	Size Nter_insert_flexibility_;
	Size Nter_loop_start_;//First flexible residue
	Size Nter_loop_end_;
	
	Size Cter_scaffold_flexibility_;
	Size Cter_insert_flexibility_;
	Size Cter_loop_start_;
	Size Cter_loop_end_; //Last flexible residue
	
	bool stop_at_closure_;
	bool final_repack_;
	
	core::Real neighbor_dis_;
	TagCOP tag_; //This is so pdb_num can be parsed at apply time instead of construction time.
	bool idealize_insert_;
	
	protocols::loops::LoopsOP loops_;
	
}; //Class AnchoredGraftMover


}// namespace grafting
}// namespace protocols

#endif  // INCLUDED_protocols_grafting_AnchoredGraftMover_HH
