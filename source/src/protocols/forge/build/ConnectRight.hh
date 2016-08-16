// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/build/ConnectRight.hh
/// @brief instruction to connect one Pose onto the right side of another
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_ConnectRight_hh
#define INCLUDED_protocols_forge_build_ConnectRight_hh

// unit headers
#include <protocols/forge/build/ConnectRight.fwd.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>

// project headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/id/NamedStubID.hh>
#include <core/kinematics/RT.hh>
#include <core/pose/Pose.hh>

// utility headers

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief instruction to connect one Pose onto the right side of another
/// @details Denote pose_left = [a,b] and pose_right = [c,d] and the newly
///  connected pose_left + pose_right = pose_total = [a,d].  If 'b' of
///  pose_left is an upper terminus or 'c' of pose_right is a lower terminus,
///  then ConnectRight will start a new chain in the Conformation at 'c' when
///  constructing pose_total, otherwise it will continue the chain at 'b'.
class ConnectRight : public BuildInstruction {


private: // typedefs


	typedef BuildInstruction Super;


public: // typedefs


	typedef Super::Size Size;

	typedef Super::ResidueTypeSetCAP ResidueTypeSetCAP;
	typedef Super::LengthEvent LengthEvent;
	typedef Super::MoveMap MoveMap;
	typedef Super::Pose Pose;

	typedef Super::Positions Positions;
	typedef Super::String String;

	typedef core::kinematics::RT RT;
	typedef core::id::NamedStubID NamedStubID;
	typedef NamedStubID::AtomList AtomNameList;


public: // construct/destruct


	/// @brief default constructor
	ConnectRight();


	/// @brief position to position jump constructor
	/// @param[in] left_position connect at this position on 'pose_left'
	///  passed into modify()
	/// @param[in] right_position connect at this position on 'pose_right'
	/// @param[in] pose_right connect this pose to the right of pose_left when
	///  modify( pose_left ) is called
	ConnectRight(
		Size const left_position,
		Size const right_position,
		Pose const & pose_right
	);


	/// @brief copy constructor
	ConnectRight( ConnectRight const & rval );


	/// @brief default destructor
	virtual
	~ConnectRight();


public: // assignment


	/// @brief copy assignment
	ConnectRight & operator =( ConnectRight const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	BuildInstructionOP clone() const;


public: // accessors


	/// @brief connect this pose to the right of pose_left when modify( pose_left )
	///  is called
	Pose const & pose_right();


	/// @brief apply the transform over the jump connection between pose_left and
	///  pose_right?  default False
	inline
	bool use_rt() const {
		return use_rt_;
	}


	/// @brief Use these atoms to compute the stub on pose_left.
	///  Default ["CA", "N", "CA", "C"].
	/// @details This is an array with 4 elements in the same order as Stub
	///  initialization: 'center atom', 'point1', 'point2', 'point3'.  See
	///  Stub docs for more details.
	inline
	AtomNameList const & left_stub_atoms() const {
		return left_stub_atoms_;
	}


	/// @brief Use these atoms to compute the stub on pose_right.
	///  Default ["CA", "N", "CA", "C"].
	/// @details This is an array with 4 elements in the same order as Stub
	///  initialization: 'center atom', 'point1', 'point2', 'point3'.  See
	///  Stub docs for more details.
	inline
	AtomNameList const & right_stub_atoms() const {
		return right_stub_atoms_;
	}


	/// @brief the rotation-translation representing an explicit transform from
	///  the stub on pose_left to the stub on pose_right, i.e. the representation
	///  of the "jump"
	inline
	RT const & rt() {
		return rt_;
	}


public: // virtual accessors


	/// @brief is the original interval storing valid information, or is empty
	///  or being used for something else?
	/// @return false, stores invalid interval
	inline
	virtual
	bool original_interval_valid() const {
		return false;
	}


	/// @brief a copy of the working range of residues specifying the modified region
	/// @details before modify() holds [1, pose_right.n_residue()]; after
	///  modify() holds [1 + pose_left.n_residue(), pose_right.n_residue + pose_left.n_residue()]
	/// @remarks This can change if listening to Conformation LengthEvents
	inline
	virtual
	Interval interval() const {
		return interval_;
	}


	/// @brief return a copy of the set of positions within the new region
	///  that were pre-existing in the original Pose prior to modify()
	/// @return An empty set -- no positions are pre-existing.
	virtual
	Positions preexisting_positions() const;


	/// @brief return a copy of the set of positions that are "new" and did
	///  not exist in the original Pose.
	/// @return A set spanning the entire interval -- all positions are new.
	virtual
	Positions new_positions() const;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has a defined conformation.  E.g. existing or copied residues.
	/// @return A set spanning the entire interval -- all positions are defined.
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


public: // mutators


	/// @brief apply the transform over the jump connection between pose_left and
	///  pose_right?  default False
	inline
	void use_rt( bool const flag ) {
		use_rt_ = flag;
	}


	/// @brief extract appropriately computed transform between two stubs that
	///  represent a jump in the given Pose between the two residues and set it
	///  as the rt for this ConnectRight
	/// @param[in] pose The Pose to use.
	/// @param[in] jump_start_residue The starting residue of the jump.
	///  This residue should be the equivalent of the jump position on pose_left.
	/// @param[in] jump_stop-residue The stopping residue of the jump.
	///  This residue should be the equivalent of the jump position on pose_right.
	/// @remarks Uses left_stub_atoms() and right_stub_atoms() as the stub atoms
	///  to compute the transform.  Remember to set use_rt() to True after
	///  calling this function if you want to actually use the transform
	///  during modify().
	void extract_rt(
		Pose const & pose,
		Size const jump_start_residue,
		Size const jump_stop_residue
	);


	/// @brief Use these atoms to compute the stub on pose_left.
	///  Default ["CA", "N", "CA", "C"].
	/// @details This is an array with 4 elements in the same order as Stub
	///  initialization: 'center atom', 'point1', 'point2', 'point3'.  See
	///  Stub docs for more details.
	inline
	void left_stub_atoms( AtomNameList const & atoms ) {
		assert( atoms.size() == 4 );
		left_stub_atoms_ = atoms;
	}


	/// @brief Use these atoms to compute the stub on pose_right.
	///  Default ["CA", "N", "CA", "C"].
	/// @details This is an array with 4 elements in the same order as Stub
	///  initialization: 'center atom', 'point1', 'point2', 'point3'.  See
	///  Stub docs for more details.
	inline
	void right_stub_atoms( AtomNameList const & atoms ) {
		assert( atoms.size() == 4 );
		right_stub_atoms_ = atoms;
	}


	/// @brief the rotation-translation representing an explicit transform from
	///  the stub on pose_left to the stub on pose_right, i.e. the representation
	///  of the "jump"
	inline
	void rt( RT const & transform ) {
		rt_ = transform;
	}


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
	/// @return empty set if no fixed positions
	virtual
	Positions original_fixed_positions() const;


	/// @brief return set of any mutable positions necessary with respect to the original
	///  interval and original Pose numbering
	/// @remarks Used for ensuring build regions for instructions do not overlap and
	///  so that jumps may be placed correctly.
	/// @return empty set if no mutable positions
	virtual
	Positions original_mutable_positions() const;


public: // virtual object descriptor


	/// @brief does this object create undefined backbone in the modified region?
	/// @return false
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
	void modify_impl( Pose & pose_left );


protected: // accessors


	/// @brief connect at this position on 'pose_left' passed into modify()
	/// @remarks this position can shift if listening to a Pose/Conformation and the number
	///  of residues changes
	inline
	Size left_position() const {
		return left_position_;
	}


	/// @brief connect at this position on 'pose_right'
	inline
	Size right_position() const {
		return right_position_;
	}


	/// @brief generate the NamedStubID for the left stub atoms on position 'p'
	inline
	NamedStubID left_named_stub_id( Size const p ) const {
		return NamedStubID( left_stub_atoms_, p );
	}


	/// @brief generate the NamedStubID for the left right atoms on position 'p'
	inline
	NamedStubID right_named_stub_id( Size const p ) const {
		return NamedStubID( right_stub_atoms_, p );
	}


protected: // mutators


	/// @brief connect at this position on 'pose_left' passed into modify()
	/// @remarks this position can shift if listening to a Pose/Conformation and the number
	///  of residues changes
	inline
	void left_position( Size const i ) {
		left_position_ = i;
	}


	/// @brief connect at this position on 'pose_right'
	inline
	void right_position( Size const i ) {
		right_position_ = i;
	}


protected: // virtual mutators


	/// @brief do the actual reset of intervals, positions, etc to initial state
	virtual
	void reset_accounting_impl();


private: // initialization


	/// @brief initialization on construction
	void init();


private: // data


	/// @brief connect at this position on 'pose_left' passed into modify()
	/// @remarks this position can shift if listening to a Pose/Conformation and the number
	///  of residues changes
	Size left_position_;


	/// @brief connect at this position on 'pose_right'
	/// @remarks this position can shift if listening to a Pose/Conformation and the number
	///  of residues changes
	Size right_position_;


	/// @brief connect this pose to the right of pose_left when modify( pose_left ) is called
	Pose pose_right_;


	/// @brief apply the transform over the jump connection between pose_left and
	///  pose_right?  default False
	bool use_rt_;


	/// @brief Use these atoms to compute the stub on pose_left.
	///  Default ["CA", "N", "CA", "C"].
	/// @details This is an array with 4 elements in the same order as Stub
	///  initialization: 'center atom', 'point1', 'point2', 'point3'.  See
	///  Stub docs for more details.
	AtomNameList left_stub_atoms_;


	/// @brief Use these atoms to compute the stub on pose_right.
	///  Default ["CA", "N", "CA", "C"].
	/// @details This is an array with 4 elements in the same order as Stub
	///  initialization: 'center atom', 'point1', 'point2', 'point3'.  See
	///  Stub docs for more details.
	AtomNameList right_stub_atoms_;


	/// @brief the rotation-translation representing an explicit transform from
	///  the stub on pose_left to the stub on pose_right, i.e. the representation
	///  of the "jump"
	RT rt_;


	/// @brief tracks the numbering of the new region specifying the connected Pose
	/// @details before modify_impl() holds [1, pose_right.n_residue()]; after
	///  modify holds [1 + pose_left.n_residue(), pose_right.n_residue + pose_left.n_residue()]
	Interval interval_;


};


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_ConnectRight_HH */
