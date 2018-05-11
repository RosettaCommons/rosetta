// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/pose/Pose.fwd.hh>


namespace protocols {
namespace grafting {


/// @brief Fairly abstract base class for GraftMover classes
class GraftMoverBase: public moves::Mover {

public:

	/// @brief Default constructor for RosettaScripts
	GraftMoverBase(std::string mover_name);

	/// @brief Start and end are the residue numbers you want your insert to go between.  start->Insert<-end
	GraftMoverBase(core::Size const start, core::Size const end, std::string mover_name);

	GraftMoverBase(core::Size const start, core::Size const end, std::string mover_name,
		core::pose::Pose const & piece, core::Size Nter_overhang_length=0, core::Size Cter_overhang_length=0);

	GraftMoverBase(GraftMoverBase const & src);

	~GraftMoverBase() override;

	//@brief copy ctor
	//GraftMoverBase( GraftMoverBase const & rhs);

	/// @brief Sets the piece that will be inserted, and any overhang residues.
	/// @details Overhang residues are residues not being inserted into the scaffold.
	/// These residues are deleted before insertion and are used by classes usually for superposition or
	/// Initial orientation of the insert relative to the scaffold.
	virtual void
	set_piece(core::pose::Pose const & piece, core::Size Nter_overhang_length, core::Size Cter_overhang_length);

	virtual void
	set_insert_region(core::Size const start, core::Size const end);

	/// @brief Copy PDBInfo from the pose piece into pose during graft.  If false(default), PDBInfo will be obsoleted.
	void
	copy_pdbinfo(bool copy_pdbinfo);

	bool
	copy_pdbinfo() const {
		return copy_pdbinfo_;
	}


public:
	/// @brief  Return the name of the Mover.
	std::string
	get_name() const override;

	Size start();
	Size original_end();
	Size end();
	Size insertion_length();

protected:

	/// @brief Steven Lewis' insertion method from insert_pose_into_pose. Wrapper to his function, using variables defined in this baseclass.
	/// @details Need to set piece to use. Deletes any overhang in piece.  Deletes any region from start to end in pose. Updates end_.
	/// Recommended use is within apply method.
	///
	core::pose::Pose
	insert_piece(core::pose::Pose const & pose, core::pose::Pose & piece);

protected:
	///Setters and accessors of private data

	void original_end(core::Size original_end);
	void insertion_length(core::Size insertion_length);
	void start(core::Size start);
	void end(core::Size end);

	core::Size Nter_overhang_length();
	void Nter_overhang_length(core::Size overhang);

	core::Size Cter_overhang_length();
	void Cter_overhang_length(core::Size overhang);

	core::pose::PoseCOP piece();
	void piece(core::pose::PoseCOP piece);

private:

	/// @brief COP of the set pose.  This is to have better support of the SavePoseMover in RosettaScripts
	core::pose::PoseCOP piece_ = nullptr;

	/// @brief Residue insertion will start from
	core::Size start_;

	/// @brief Residue insertion will end before here. Updates after insertion.
	core::Size end_;

	/// @brief some functions need to only work on the original numbers. (combine movemaps)
	core::Size original_end_;

	core::Size insertion_length_;

	/// @brief Number of overhang residues on N terminus.  Updates on delete_overhang_residues
	core::Size Nter_overhang_length_;

	/// @brief Number of overhang residues on C terminus.  Updates on delete_overhang_residues
	core::Size Cter_overhang_length_;

	bool copy_pdbinfo_;

};  // class GraftMoverBase

}  // namespace grafting
}  // namespace protocols

#endif  // INCLUDED_protocols_grafting_GraftMoverBase_HH
