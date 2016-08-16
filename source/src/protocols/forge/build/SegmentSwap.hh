// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/SegmentSwap.hh
/// @brief instruction to swap a segment with an external segment
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_SegmentSwap_hh
#define INCLUDED_protocols_forge_build_SegmentSwap_hh

// unit headers
#include <protocols/forge/build/SegmentSwap.fwd.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>

// project headers
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief instruction to swap a segment with an external segment
class SegmentSwap : public BuildInstruction {


private: // typedefs


	typedef BuildInstruction Super;


public: // typedefs


	typedef Super::Size Size;

	typedef Super::ResidueTypeSetCAP ResidueTypeSetCAP;
	typedef Super::LengthEvent LengthEvent;
	typedef Super::Pose Pose;

	typedef Super::Positions Positions;
	typedef Super::String String;

	typedef core::kinematics::MoveMap MoveMap;


public: // construct/destruct


	/// @brief default constructor
	SegmentSwap();


	/// @brief constructor
	/// @param[in] interval swap out this range of residues
	/// @param[in] move_map fixed backbone residues in this movemap will be used for new jumps
	/// @param[in] swap_in swap in this pose
	SegmentSwap(
		Interval const & i,
		MoveMap const & swap_in_movemap,
		Pose const & swap_in
	);


	/// @brief copy constructor
	SegmentSwap( SegmentSwap const & rval );


	/// @brief default destructor
	virtual
	~SegmentSwap();


public: // assignment


	/// @brief copy assignment
	SegmentSwap & operator =( SegmentSwap const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	BuildInstructionOP clone() const;


public: // accessors


	/// @brief fixed backbone residues in this movemap will be used for new jumps
	/// @remarks Procedure will attempt to honor this movemap as much as it can.
	///  The caveat is that sequences of calls to some FoldTree routines may shift
	///  the jumps internally in a way that is not easily predictable.  If the
	///  procedure cannot  find an allowed residue for a jump, it will make a jump
	///  to the (lower) median residue in the disconnected fold tree interval.
	MoveMap const & swap_in_movemap() const;


	/// @brief the pose to swap in
	Pose const & swap_in() const;


public: // virtual accessors


	/// @brief is the original interval storing valid information, or is empty
	///  or being used for something else?
	/// @return true, stores valid interval
	inline
	virtual
	bool original_interval_valid() const {
		return true;
	}


	/// @brief a copy of the working range of residues specifying the swapped region
	/// @details This residue range can change wrt length changes in Pose /Conformation
	///  being watched.
	virtual
	Interval interval() const;


	/// @brief return a copy of the set of positions within the new region
	///  that were pre-existing in the original Pose prior to modify()
	/// @return An empty set -- no positions are pre-existing.
	virtual
	Positions preexisting_positions() const;


	/// @brief return a copy of the set of positions that are "new" and did
	///  not exist in the original Pose.
	/// @return A set of positions spanning the interval -- all positions are
	///  are defined.
	virtual
	Positions new_positions() const;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has a defined conformation.  E.g. existing or copied residues.
	/// @return A set of positions spanning the interval -- all positions are
	///  are defined.
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	Positions defined_positions() const;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has an undefined conformation.  E.g. newly created residues.
	/// @return An empty set -- no undefined positions.
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	Positions undefined_positions() const;


	/// @brief return a copy of the MoveMap that defines the moveable/fixed
	///  positions/dofs for this instruction
	/// @return a MoveMap with [interval.left, interval.right] bb set to false
	///  at the MoveMapTorsionID level
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	MoveMap movemap() const;


public: // virtual Conformation observer interface


	/// @brief update indexing on residue append
	virtual
	void on_residue_append( LengthEvent const & event );


	/// @brief update indexing on residue prepend
	virtual
	void on_residue_prepend( LengthEvent const & event );


	/// @brief update indexing on residue delete
	virtual
	void on_residue_delete( LengthEvent const & event );


public: // original positions


	/// @brief return the set of positions within the original interval that
	///  will be kept in this BuildInstruction
	/// @return An empty set -- no positions are kept.
	virtual
	Positions original_kept_positions() const;


	/// @brief return set of positions within the original interval that will
	///  be deleted in this BuildInstruction
	/// @return An empty set -- no positions are deleted.
	virtual
	Positions original_deleted_positions() const;


public: // instruction comparison


	/// @brief return set of any fixed positions necessary with respect to the original
	///  interval and original Pose numbering
	/// @remarks Used for ensuring build regions for instructions do not overlap and
	///  so that jumps may be placed correctly.
	/// @return empty set if no fixed positions necessary
	virtual
	Positions original_fixed_positions() const;


	/// @brief return set of any mutable positions necessary with respect to the original
	///  interval and original Pose numbering
	/// @remarks Used for ensuring build regions for instructions do not overlap and
	///  so that jumps may be placed correctly.
	virtual
	Positions original_mutable_positions() const;


public: // virtual object descriptor


	/// @brief does this object create undefined backbone in the modified region?
	inline
	virtual
	bool creates_undefined_backbone() const {
		return false;
	}


protected: // virtual Pose modification methods


	/// @brief are dependencies satisfied so that modify_impl() can complete
	///  successfully?
	/// @return always True, this BuildInstruction has no dependencies
	inline
	virtual
	bool dependencies_satisfied() const {
		return true;
	}


	/// @brief do the actual work of modifying the Pose
	virtual
	void modify_impl( Pose & pose );


protected: // virtual mutators


	/// @brief do the actual reset of intervals, positions, etc to initial state
	virtual
	void reset_accounting_impl();


private: // init


	/// @brief init to be called during non-default constructors
	void init();


private: // data


	/// @brief range of residues to swap out
	/// @remarks this range can shift if listening to a Pose/Conformation and the number
	///  of residues changes
	Interval interval_;


	/// @brief fixed backbone residues in this movemap will be used for new jumps
	/// @remarks Procedure will attempt to honor this movemap as much as it can.
	///  The caveat is that sequences of calls to some FoldTree routines may shift
	///  the jumps internally in a way that is not easily predictable.  If the
	///  procedure cannot  find an allowed residue for a jump, it will make a jump
	///  to the (lower) median residue in the disconnected fold tree interval.
	MoveMap swap_in_movemap_;


	/// @brief swap in this Pose
	Pose swap_in_;


};


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_SegmentSwap_HH */
