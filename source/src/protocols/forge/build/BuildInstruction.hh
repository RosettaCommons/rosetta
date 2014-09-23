// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/forge/build/BuildInstruction.hh
/// @brief tracks modifications to be made and is capable of residue length changes on a Pose
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_BuildInstruction_hh
#define INCLUDED_protocols_forge_build_BuildInstruction_hh

// unit headers
#include <protocols/forge/build/BuildInstruction.fwd.hh>

// package headers
#include <protocols/forge/build/Interval.hh>

// project headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// utility headers
#include <utility/signals/Link.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

// C++ headers
#include <set>
#include <string>

#include <utility/vector1.hh>



namespace protocols {
namespace forge {
namespace build {


// enclose enum in descriptive namespace to prevent conflicts
namespace BuildInstructionState {
	/// @brief describes the state of the BuildInstruction
	/// @details There are currently three possible BuildInstruction states.
	///  <em>READY</em> indicates the BuildInstruction has been reset and is
	///  ready to modify a Pose.  <em>WAITING_ON_DEPENDENCIES</em> indicates
	///  the BuildInstruction is waiting for its dependencies to be satisfied
	///  before allowing modifications to proceed.  <em>MODIFY_WAS_SUCCESSFUL</em>
	///  indicates the BuildInstruction has finished modifications to the Pose,
	///  and its residue indexing is now consistent with the newly modified Pose.
	enum Enum {
		READY,
		MODIFY_WAS_SUCCESSFUL,
		WAITING_ON_DEPENDENCIES
	};
}


/// @brief tracks modifications to be made and is capable of making residue length changes on a Pose
class BuildInstruction : public utility::pointer::ReferenceCount {


private: // typedefs


	typedef utility::pointer::ReferenceCount Super;


public: // typedefs


	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::chemical::ResidueTypeSet ResidueTypeSet;
	typedef core::chemical::ResidueTypeSetCAP ResidueTypeSetCAP;
	typedef core::chemical::ResidueTypeSetCOP ResidueTypeSetCOP;
	typedef core::conformation::signals::LengthEvent LengthEvent;
	typedef core::kinematics::MoveMap MoveMap;
	typedef core::pose::Pose Pose;

	typedef utility::vector1< BuildInstructionCAP > BuildInstructionCAPs;

	typedef utility::signals::Link Link;

	typedef std::set< Size > Positions;
	typedef std::string String;


public: // construct/destruct


	/// @brief default constructor
	BuildInstruction();


	/// @brief interval constructor
	/// @param[in] i the residue range to operate on
	/// @param[in] rts the residue type set to use, default FA_STANDARD
	BuildInstruction(
		Interval const & i,
		ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
	);


	/// @brief copy constructor
	BuildInstruction( BuildInstruction const & rval );


	/// @brief default destructor
	virtual
	~BuildInstruction();


public: // assignment


	/// @brief copy assignment
	BuildInstruction & operator =( BuildInstruction const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	BuildInstructionOP clone() const = 0;


public: // Pose modification methods


	/// @brief modify this pose
	/// @return the state of the BuildInstruction
	BuildInstructionState::Enum modify( Pose & pose );


public: // query state


	/// @brief return the state of this BuildInstruction
	inline
	BuildInstructionState::Enum state() const {
		return state_;
	}


	/// @brief Is the BuildInstruction's state at <em>READY</em>?
	/// @remarks <em>READY</em> indicates the BuildInstruction has been reset
	///  and is ready to modify a Pose.
	bool ready() const;


	/// @brief Is the BuildInstruction's state at <em>WAITING_ON_DEPENDENCIES</em>?
	/// @remarks <em>WAITING_ON_DEPENDENCIES</em> indicates
	///  the BuildInstruction is waiting for its dependencies to be satisfied
	///  before allowing modifications to proceed.
	bool waiting_on_dependencies() const;


	/// @brief Is the BuildInstruction's state at <em>MODIFY_WAS_SUCCESSFUL</em>?
	/// @remarks <em>MODIFY_WAS_SUCCESSFUL</em> indicates the BuildInstruction has
	///  finished modifications to the Pose, and its residue indexing is now
	///  consistent with the newly modified Pose.
	bool modify_was_successful() const;


public: // virtual Conformation observer interface


	/// @brief update indexing on residue append
	virtual
	void on_residue_append( LengthEvent const & event ) = 0;


	/// @brief update indexing on residue prepend
	virtual
	void on_residue_prepend( LengthEvent const & event ) = 0;


	/// @brief update indexing on residue delete
	virtual
	void on_residue_delete( LengthEvent const & event ) = 0;


public: // common Conformation observer interface


	/// @brief attach to a Pose's conformation
	void attach_to( Pose & pose );


	/// @brief detach from a Pose's conformation
	void detach_from();


	/// @brief update any indexing wrt length change to Pose/Conformation being watched
	void on_length_change( LengthEvent const & event );


public: // accessors


	/// @brief detach after modify()?
	/// @details Ordinarily a modify() call consists of attaching as a
	///  length obs of the Pose, actually doing the modify, and then
	///  detaching as a length obs of the Pose.  This flag governs whether
	///  the last detach occurs.
	inline
	bool detach_after_modify() const {
		return detach_after_modify_;
	}


	/// @brief return original residue range of this instruction
	inline
	Interval const & original_interval() const {
		return original_interval_;
	}


	/// @brief the residue type set being used
	inline
	ResidueTypeSet const & residue_type_set() const {
		ResidueTypeSetCOP rts( rts_ ); // Fix me: returning reference to temporairly locked pointer
		return *rts;
	}


	/// @brief does this BuildInstruction have dependencies?
	inline
	bool has_dependencies() const {
		return !dependencies_.empty();
	}


	/// @brief the number of dependencies this BuildInstruction has
	inline
	Size n_dependencies() const {
		return dependencies_.size();
	}


	/// @brief the list of instructions whose modify() must complete before
	///  the modify() for this instruction may be called successfully
	inline
	BuildInstructionCAPs const & dependencies() const {
		return dependencies_;
	}


public: // virtual accessors


	/// @brief is the original interval storing valid information, or is empty
	///  or being used for something else?
	virtual
	bool original_interval_valid() const = 0;


	/// @brief a copy of the working range of residues specifying the modified region
	/// @details This residue range can change wrt length changes in Pose /Conformation
	///  being watched.
	virtual
	Interval interval() const = 0;


	/// @brief return a copy of the set of positions within the new region
	///  that were pre-existing in the original Pose prior to modify()
	virtual
	Positions preexisting_positions() const = 0;


	/// @brief return a copy of the set of positions that are "new" and did
	///  not exist in the original Pose.
	virtual
	Positions new_positions() const = 0;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has a defined conformation.  E.g. existing or copied residues.
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	Positions defined_positions() const = 0;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has an undefined conformation.  E.g. newly created residues.
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	Positions undefined_positions() const = 0;


	/// @brief return a copy of the MoveMap that defines the moveable/fixed
	///  positions/dofs for this instruction
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	MoveMap movemap() const = 0;


public: // mutators


	/// @brief set detach after modify()
	/// @details Ordinarily a modify() call consists of attaching as a
	///  length obs of the Pose, actually doing the modify, and then
	///  detaching as a length obs of the Pose.  This flag governs whether
	///  the last detach occurs.
	inline
	void detach_after_modify( bool const flag ) {
		detach_after_modify_ = flag;
	}


	/// @brief reset intervals, positions, etc to initial state and drop
	///  observers.  State set to READY.
	void reset_accounting();


public: // virtual mutators


	/// @brief add an instruction to this BuildInstruction's dependency list
	/// @remarks before this instruction's modify() may complete properly, all
	///  instructions in the dependency list must have completed modify() first
	virtual
	void add_dependency_to( BuildInstructionCAP i );


	/// @brief is this instruction dependent upon the given instruction?
	virtual
	bool dependent_on( BuildInstructionCAP i ) const;


	/// @brief clear the list of dependencies
	inline
	void clear_dependencies() {
		dependencies_.clear();
	}


public: // original positions


	/// @brief return the set of positions within the original interval that
	///  will be kept in this BuildInstruction
	virtual
	Positions original_kept_positions() const = 0;


	/// @brief return set of positions within the original interval that will
	///  be deleted in this BuildInstruction
	virtual
	Positions original_deleted_positions() const = 0;


public: // instruction comparison


	/// @brief compares fixed and mutable positions to determine compatibility with
	///  another instruction
	/// @remarks override this to obtain custom compatibility check
	virtual
	bool compatible_with( BuildInstruction const & rval ) const;


	/// @brief return set of any fixed positions necessary with respect to the original
	///  interval and original Pose numbering
	/// @remarks Used for ensuring build regions for instructions do not overlap and
	///  so that jumps may be placed correctly.
	/// @return empty set if no fixed positions necessary
	virtual
	Positions original_fixed_positions() const = 0;


	/// @brief return set of any mutable positions necessary with respect to the original
	///  interval and original Pose numbering
	/// @remarks Used for ensuring build regions for instructions do not overlap and
	///  so that jumps may be placed correctly.
	virtual
	Positions original_mutable_positions() const = 0;


public: // virtual object descriptor


	/// @brief does this object create undefined backbone in the modified region?
	virtual
	bool creates_undefined_backbone() const = 0;


protected: // accessors


	/// @brief access to the Conformation length observer link
	inline
	Link const & length_obs_link() const {
		return length_obs_link_;
	}


protected: // dependency handling


	/// @brief are dependencies satisfied so that modify_impl() can complete
	///  successfully?
	/// @return True if modify_impl() can be called, False otherwise.
	/// @remarks This function will automatically be checked within modify()
	///  prior to running modify_impl().
	virtual
	bool dependencies_satisfied() const;


protected: // virtual Pose modification methods


	/// @brief do the actual work of modifying the Pose
	virtual
	void modify_impl( Pose & pose ) = 0;


protected: // set state


	/// @brief set the BuildInstruction's current state
	inline
	void state( BuildInstructionState::Enum const s ) {
		state_ = s;
	}


protected: // virtual mutators


	/// @brief do the actual reset of intervals, positions, etc to initial state
	virtual
	void reset_accounting_impl() = 0;


private: // data


	/// @brief Describes the state of the BuildInstruction.
	/// @remarks See documentation of BuildInstructionState::Enum and the
	///  query state, set state, and <tt>modify()</tt> functions for more
	///  information.
	BuildInstructionState::Enum state_;


	/// @brief defines the original residue range via closed interval [first, second]
	///  or an original anchor position [anchor, anchor] that's marked as invalid
	/// @remarks no mutable access to this data
	Interval original_interval_;


	/// @brief the residue type set to use
	ResidueTypeSetCAP rts_;


	/// @brief Conformation length observer link
	Link length_obs_link_;


	/// @brief detach after performing modify()
	/// @details Ordinarily a modify() call consists of attaching as a
	///  length obs of the Pose, actually doing the modify, and then
	///  detaching as a length obs of the Pose.  This flag governs whether
	///  the last detach occurs.
	bool detach_after_modify_;


	/// @brief if this BuildInstruction's modify() is dependent on any other
	///  BuildInstructions whose modify() must complete first, we store access
	///  to those instructions here
	/// @warning Due to copy ambiguities, dependencies are currently not copied!
	///  They must be externally rebuilt, e.g. by the BuildManager.
	BuildInstructionCAPs dependencies_;


};


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_BuildInstruction_HH */
