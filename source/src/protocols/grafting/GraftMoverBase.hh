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
	using core::Size;
	using core::pose::Pose;
	using core::pose::PoseOP;
	using protocols::loops::Loop;
	using core::kinematics::MoveMapOP;
	using core::kinematics::MoveMap;
	
///@brief Fairly abstract base class for GraftMover classes
class GraftMoverBase: public moves::Mover {

public:
	///@brief Start and end are the residue numbers you want your insert to go between.  start->Insert<-end
	GraftMoverBase(Size const start, Size const end, std::string mover_name);

	GraftMoverBase(Size const start, Size const end, std::string mover_name,
		core::pose::Pose const & piece, Size Nter_overhang_length=0, Size Cter_overhang_length=0);
        
	virtual ~GraftMoverBase();

	//@brief copy ctor
	//GraftMoverBase( GraftMoverBase const & rhs);

	///@brief Sets the piece that will be inserted, and any overhang residues.
	///@details Overhang residues are residues not being inserted into the scaffold. 
	/// These residues are deleted before insertion and are used by classes usually for superposition or
	/// Initial orientation of the insert relative to the scaffold.  
	virtual void 
	set_piece(Pose const & piece, Size Nter_overhang_length, Size Cter_overhang_length);

	virtual void
	set_insert_region(Size const start, Size const end);
	
public:
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
	
protected:
	///Setters and accessors of private data
	void original_end(Size original_end);
	Size original_end();
	
	void insertion_length(Size insertion_length);
	Size insertion_length();
	
	Size start();
	Size end();
	Size Nter_overhang_length();
	Size Cter_overhang_length();
	PoseOP piece();
	
private:
   	//Reference of the pose piece.  Should be changed to local copy, but I'm not sure how to do that.
	PoseOP piece_;

	///@brief Residue insertion will start from
	Size start_;
	///@brief Residue insertion will end before here. Updates after insertion.
	Size end_;
	///@brief some functions need to only work on the original numbers. (combine movemaps)
	Size original_end_;

	Size insertion_length_;
	
	///@brief Number of overhang residues on N terminus.  Updates on delete_overhang_residues
	Size Nter_overhang_length_;
	///@brief Number of overhang residues on C terminus.  Updates on delete_overhang_residues
	Size Cter_overhang_length_;
	
};  // class GraftMoverBase

}  // namespace grafting
}  // namespace protocols

#endif  // INCLUDED_protocols_grafting_GraftMoverBase_HH
