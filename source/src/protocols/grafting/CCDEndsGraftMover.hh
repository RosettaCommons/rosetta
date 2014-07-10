// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/grafting/CCDEndsGraftMover.hh
/// @brief   Class to graft a piece into a pose.
/// @Author  Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_grafting_CCDEndsGraftMover_HH
#define INCLUDED_protocols_grafting_CCDEndsGraftMover_HH


//Unit Headers
#include <protocols/grafting/GraftMoverBase.hh>
#include <protocols/grafting/AnchoredGraftMover.hh>
#include <protocols/grafting/CCDEndsGraftMover.fwd.hh>
#include <protocols/moves/Mover.hh>

//Core
#include <core/pose/Pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>

//Protocols
#include <protocols/simple_moves/BackboneMover.fwd.hh>

namespace protocols {
namespace grafting {
	
	using protocols::simple_moves::SmallMoverOP;
	
///@brief General purpose Grafting class which:
/// 1) superimposes the insert onto the scaffold using any overhang residues,
/// 2) Inserts the pose piece into the scaffold pose, deleting any overhang residues or residues in the region the isertion will occur between.
/// 3) Cycles of:
///      a) SmallMover for sampling (Can be disabled))
///      b) CCDs both terminal ends of the scaffold using flexibility settings or movemap
///          To close the graft.  
///      c) MinMover to help close the graft through chainbreak terms
///      d) Closure check - will return if closed (Can be disabled)
///      d) MonteCarlo Boltzmann criterion
///
/// 5) Repacks flexible residues and insert residues
///
///
/// Differs from AnchoredGraftMover in that the insert is held fixed in space after the superposition
///
/// example:
///    mover = CCDEndsGraftMover(start, end, piece, cter_overhang, nter_overhang)
///    mover.apply(pose)
///
/// see also: grafting/util.hh.
///
///@details Uses two CCD arms on either side of the scaffold to close graft.
/// Insert is held fixed in cartesian space after superposition of any overhang residues by default.
///
/// ****Nter_loop_start-----> | Piece | <----Nter_loop_end****
/// Default flexibility on Nter and Cter is only two residues (--> part of diagram).
/// Will delete any residues between start and end of the scaffold, and any overhang residues from the insert.
///
/// Internally, apply performs the insertion, idealizes the loop residues (omegas to 180, peptide bonds idealized) and the newly made polymer connections at the insert point, and then attempts to close the loop(s).
/// Repacks sidechains of insert and overhang residues after insertion is complete
/// By default, will stop apply once graft is closed.  Closure is measured as all bond lengths and angles of loop being acceptable.
///
///
class CCDEndsGraftMover : public protocols::grafting::AnchoredGraftMover {
	
public:
    
	CCDEndsGraftMover();
	
	///@brief Start and end are the residue numbers you want your insert to go between.  start->Insert<-end
	CCDEndsGraftMover(Size const start, Size const end, bool copy_pdbinfo = false);

	CCDEndsGraftMover(
		Size const start, Size const end, 
		core::pose::Pose const & piece, Size Nter_overhang=2, Size Cter_overhang=2, bool copy_pdbinfo = false);
    
	CCDEndsGraftMover(CCDEndsGraftMover const & src);

	virtual ~CCDEndsGraftMover();
    
	
	virtual void
	apply(core::pose::Pose & pose);
	
public:
	
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	///@brief Advanced way to set flexibility. 
	///@details Will combine the movemaps for apply, and renumber everything. Flexible residues in multiple chains not recommended.
	/// One arm Will go from first flexible N terminal residue to after start_+any contiguous residues on in the movemap from there. Opposite for Cter side.
	/// This way any loop regions within a chain on either side can be used as flexible residues to close full graft.
	///
	/// Note: Will disregard flexibility settings, as the movemaps will be used as primary way to define flexibility.
	/// May want to consider turning off the sampling step when passing crazy movemaps.
	///
	virtual void 
	set_movemaps(MoveMapCOP const scaffold_mm, MoveMapCOP const insert_mm);
	

public:
	
	std::string
	get_name();
	
	protocols::moves::MoverOP
	clone() const;
	
	protocols::moves::MoverOP
	fresh_instance() const;
	

private:
	
	virtual SmallMoverOP
	setup_default_small_mover();
	
}; //Class CCDEndsGraftMover


}// namespace grafting
}// namespace protocols

#endif  // INCLUDED_protocols_grafting_CCDEndsGraftMover_HH


