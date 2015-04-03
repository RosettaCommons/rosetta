// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/forge/build/SegmentRebuild.hh
/// @brief instruction to rebuild a segment
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_SegmentRebuild_hh
#define INCLUDED_protocols_forge_build_SegmentRebuild_hh

// unit headers
#include <protocols/forge/build/SegmentRebuild.fwd.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>

// project headers
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace forge {
namespace build {


/// @brief instruction to rebuild a segment
/// @remarks Handles both cut based rebuilding (i.e. loops) and continuous
///  rebuilding adjacent to the boundary of a chain, such as n-term/c-term
///  extensions.
class SegmentRebuild : public BuildInstruction {


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
	SegmentRebuild();


	/// @brief sec.struct only constructor (poly-alanine for new region)
	/// @param[in] interval rebuild this range of residues
	/// @param[in] ss the secondary structure desired, also defines length of new build region
	/// @param[in] rts the residue type set to use, default FA_STANDARD
	/// @param[in] keep_known_bb_torsions_at_junctions Attempt to keep the omega
	///  at original_interval().left-1, the phi at original_interval().left, and
	///  the psi+omega at original_interval().right present from the original Pose
	///  in the modified Pose.
	/// @remarks length of the *one-letter* aa must equal the length of ss
	SegmentRebuild(
		Interval const & i,
		String const & ss,
		ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ),
		bool const keep_known_bb_torsions_at_junctions = false
	);


	/// @brief full constructor
	/// @param[in] interval rebuild this range of residues
	/// @param[in] ss the secondary structure desired, also defines length of new build region
	/// @param[in] aa the annotated amino acid sequence desired, default is poly-alanine
	/// @param[in] rts the residue type set to use, default FA_STANDARD
	/// @param[in] keep_known_bb_torsions_at_junctions Attempt to keep the omega
	///  at original_interval().left-1, the phi at original_interval().left, and
	///  the psi+omega at original_interval().right present from the original Pose
	///  in the modified Pose.
	/// @remarks length of the *one-letter* aa must equal the length of ss
	SegmentRebuild(
		Interval const & i,
		String const & ss,
		String const & aa,
		ResidueTypeSetCAP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ),
		bool const keep_known_bb_torsions_at_junctions = false
	);


	/// @brief copy constructor
	SegmentRebuild( SegmentRebuild const & rval );


	/// @brief default destructor
	virtual
	~SegmentRebuild();


public: // assignment


	/// @brief copy assignment
	SegmentRebuild & operator =( SegmentRebuild const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	BuildInstructionOP clone() const;


public: // accessors


	/// @brief get secondary structure string
	inline
	String const & ss() const {
		return ss_;
	}


	/// @brief get annotated amino acid string
	inline
	String const & aa() const {
		return aa_;
	}


	/// @brief Attempt to keep the omega at original_interval().left-1, the phi
	///  at original_interval().left and the psi+omega at original_interval().right
	///  present from the original Pose in the modified Pose?  Default False
	inline
	bool keep_known_bb_torsions_at_junctions() const {
		return keep_known_bb_torsions_at_junctions_;
	}


public: // virtual accessors


	/// @brief is the original interval storing valid information, or is empty
	///  or being used for something else?
	/// @return true, stores valid interval
	inline
	virtual
	bool original_interval_valid() const {
		return true;
	}


	/// @brief a copy of the working range of residues specifying the modified region
	/// @remarks this can change if listening to Conformation LengthEvents
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
	/// @return A set of positions spanning the entire modified interval -- all
	///  positions are undefined.
	virtual
	Positions new_positions() const;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has a defined conformation.  E.g. existing or copied residues.
	/// @return An empty set -- no positions are defined.
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	Positions defined_positions() const;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has an undefined conformation.  E.g. newly created residues.
	/// @return A set of positions spanning the entire modified interval -- all
	///  positions are undefined.
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
	/// @return An empty set -- no positions are kept.
	virtual
	Positions original_kept_positions() const;


	/// @brief return set of positions within the original interval that will
	///  be deleted in this BuildInstruction
	/// @return A set containing all positions in the original interval.
	virtual
	Positions original_deleted_positions() const;


public: // instruction comparison


	/// @brief return set of any fixed positions necessary with respect to the original
	///  interval and original Pose numbering
	/// @remarks Used for ensuring build regions for instructions do not overlap and
	///  so that jumps may be placed correctly.  If SegmentRebuild has dependencies,
	///  the set of positions shrinks from the endpoints of [original_left - 1, original_right + 1]
	///  down to an empty set in the assumption that the user is requesting an
	///  advanced feature and knows what they're doing around the endpoints,
	///  e.g. rebuilding directly adjacent to a swapped out section.
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
	/// @return true
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


	/// @brief the working range of residues
	/// @remarks this range can shift if listening to a Pose/Conformation and the number
	///  of residues changes
	Interval interval_;


	/// @brief secondary structure string, also defines length of extension
	String ss_;


	/// @brief annotated amino acid string, length of the one-letter version
	///  must be equal to length of ss
	String aa_;


	/// @brief Attempt to keep the omega at original_interval().left-1, the phi
	///  at original_interval().left and the psi+omega at original_interval().right
	///  present from the original Pose in the modified Pose?
	/// @details If True, during modify(), will (1) set the omega of interval_.left-1
	///  in the newly modified Pose equal to the omega of the original Pose at
	///  original_interval().left-1, (2) set the phi of interval_.left in
	///  the newly modified Pose equal to the phi of the original Pose in
	///  original_interval().right, and (3) set the psi+omega of interval_.right
	///  in the newly modified Pose equal to the original psi+omega of
	///  original_interval().right.  Default False.
	bool keep_known_bb_torsions_at_junctions_;


};


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_SegmentRebuild_HH */
