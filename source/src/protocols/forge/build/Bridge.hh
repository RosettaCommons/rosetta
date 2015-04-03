// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/build/Bridge.hh
/// @brief  connect two contiguous but disjoint sections of a Pose into one
///         continuous section
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_Bridge_hh
#define INCLUDED_protocols_forge_build_Bridge_hh

// unit headers
#include <protocols/forge/build/Bridge.fwd.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief connect two contiguous but disjoint sections of a Pose into one
///  continuous section
/// @details Anchor residues [i,j] to bridge must be adjacent to each other in the
///  Pose (i+1 = j) and adjacent to a cutpoint.  Both i and j will be idealized
///  and marked as new moveable positions since the psi @ i and the phi @ j are
///  undefined.  Depending on the fold tree, this will cause a random swing
///  either downstream or upstream of the bridge!
class Bridge : public BuildInstruction {


private: // typedefs


	typedef BuildInstruction Super;


public: // typedefs


	typedef Super::Size Size;
	typedef Super::Real Real;

	typedef Super::ResidueTypeSetCAP ResidueTypeSetCAP;
	typedef Super::LengthEvent LengthEvent;
	typedef Super::MoveMap MoveMap;
	typedef Super::Pose Pose;

	typedef Super::Positions Positions;
	typedef Super::String String;


public: // construct/destruct


	/// @brief default constructor
	Bridge();


	/// @brief sec.struct only constructor (poly-alanine for new region)
	/// @param[in] interval build bridge using these two residues as anchor positions
	/// @param[in] ss the secondary structure desired, also defines length of new bridge,
	///  region between the anchor positions, can be empty
	/// @param[in] rts the residue type set to use, default FA_STANDARD
	/// @remarks length of the *one-letter* aa must equal the length of ss
	Bridge(
		Interval const & i,
		String const & ss,
		ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
	);


	/// @brief full constructor
	/// @param[in] interval build bridge using these two residues as anchor positions
	/// @param[in] ss the secondary structure desired, also defines length of new bridge,
	///  region between the anchor positions, can be empty
	/// @param[in] aa the annotated amino acid sequence desired, default is poly-alanine
	/// @param[in] rts the residue type set to use, default FA_STANDARD
	/// @remarks length of the *one-letter* aa must equal the length of ss
	Bridge(
		Interval const & i,
		String const & ss,
		String const & aa,
		ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )
	);


	/// @brief copy constructor
	Bridge( Bridge const & rval );


	/// @brief default destructor
	virtual
	~Bridge();


public: // assignment


	/// @brief copy assignment
	Bridge & operator =( Bridge const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	BuildInstructionOP clone() const;


public: // virtual accessors


	/// @brief is the original interval storing valid information, or is empty
	///  or being used for something else?
	/// @return true, stores valid interval
	inline
	virtual
	bool original_interval_valid() const {
		return true;
	}


	/// @brief a copy of the working range of residues specifying the bridged region
	///  including the anchors
	/// @details This residue range can change wrt length changes in Pose /Conformation
	///  being watched.
	inline
	virtual
	Interval interval() const {
		return interval_;
	}


	/// @brief return a copy of the set of positions within the new region
	///  that were pre-existing in the original Pose prior to modify()
	/// @return A set containing two positions -- interval.left and interval.right.
	virtual
	Positions preexisting_positions() const;


	/// @brief return a copy of the set of positions that are "new" and did
	///  not exist in the original Pose.
	/// @return A set containing positions spanning [interval.left+1, interval.right-1].
	virtual
	Positions new_positions() const;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has a defined conformation.  E.g. existing or copied residues.
	/// @return A set containing two positions -- interval.left and interval.right.
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	Positions defined_positions() const;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has an undefined conformation.  E.g. newly created residues.
	/// @return A set containing positions spanning [interval.left+1, interval.right-1].
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	Positions undefined_positions() const;


	/// @brief return a copy of the MoveMap that defines the moveable/fixed
	///  positions/dofs for this instruction
	/// @return a MoveMap with [interval.left, interval.right] bb & chi set to true
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
	/// @return A set containing the endpoints of the original interval.
	virtual
	Positions original_kept_positions() const;


	/// @brief return set of positions within the original interval that will
	///  be deleted in this BuildInstruction
	/// @return A set containing the positions in [original.left+1, original.right-1].
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
		return true;
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


private: // data


	/// @brief the anchor positions of the bridge
	/// @remarks this can shift if listening to a Pose/Conformation and the number
	///  of residues changes
	Interval interval_;


	/// @brief secondary structure string, also defines length of the bridge
	String ss_;


	/// @brief annotated amino acid string, length of the one-letter version
	///  must be equal to length of ss
	String aa_;


};


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_Bridge_HH */
