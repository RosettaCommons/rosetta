// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/forge/build/SegmentInsert.hh
/// @brief  insert an external segment flanked by new regions
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_build_SegmentInsert_hh
#define INCLUDED_protocols_forge_build_SegmentInsert_hh

// unit headers
#include <protocols/forge/build/SegmentInsert.fwd.hh>

// package headers
#include <protocols/forge/build/BuildInstruction.hh>

// project headers
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>

#include <protocols/forge/build/SegmentInsertConnectionScheme/SegmentInsertConnectionScheme_Enum.hh>

namespace protocols {
namespace forge {
namespace build {


/// @brief insert an external segment flanked by new regions
/// @remarks Similar to SegmentRebuild, but with the addition of a defined
///  segment located in the modified region.
class SegmentInsert : public BuildInstruction {


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
	SegmentInsert();


	/// @brief sec.struct only constructor (poly-ala for flanking regions)
	/// @param[in] interval The interval between which the insert will span.
	///  To perform a pure insertion without replacing any residues
	///  within a region, use an interval with a zero as the left endpoint, e.g.
	///  [0, insert_after_this_residue].  If inserting before the first residue
	///  of the Pose then interval = [0,0].  If inserting after the last residue
	///  of the Pose then interval = [0, last_residue].
	/// @param[in] ss The secondary structure specifying the flanking regions,
	///  with a character '^' specifying where the insert is to be placed.
	/// @param[in] insert The Pose to insert.
	/// @param[in] keep_known_bb_torsions_at_junctions Attempt to keep the omega
	///  at original_interval().left-1, the phi at original_interval().left, and
	///  the psi+omega at original_interval().right present from the original Pose
	///  in the modified Pose.  This should be false for pure insertions.
	/// @param[in] connection_scheme Connect insertion on its N-side, C-side,
	///  or decide randomly between the two (default RANDOM).
	SegmentInsert(
		Interval const & i,
		String const & ss,
		Pose const & insert,
		bool const keep_known_bb_torsions_at_junctions = false,
		SegmentInsertConnectionScheme::Enum connection_scheme = protocols::forge::build::SegmentInsertConnectionScheme::RANDOM_SIDE
	);


	/// @brief sec.struct + aa constructor
	/// @param[in] interval The interval between which the insert will span.
	///  To perform a pure insertion without replacing any residues
	///  within a region, use an interval with a zero as the left endpoint, e.g.
	///  [0, insert_after_this_residue].  If inserting before the first residue
	///  of the Pose then interval = [0,0].  If inserting after the last residue
	///  of the Pose then interval = [0, last_residue].
	/// @param[in] ss The secondary structure specifying the flanking regions,
	///  with a character '^' specifying where the insert is to be placed.
	/// @param[in] aa The annotated amino acid string specifying the flanking
	///  regions, with a character '^' specifying where the insert is to be
	///  placed.
	/// @param[in] insert The Pose to insert.
	/// @param[in] keep_known_bb_torsions_at_junctions Attempt to keep the omega
	///  at original_interval().left-1, the phi at original_interval().left, and
	///  the psi+omega at original_interval().right present from the original Pose
	///  in the modified Pose.  This should be false for pure insertions.
	/// @param[in] connection_scheme Connect insertion on its N-side, C-side,
	///  or decide randomly between the two (default RANDOM).
	/// @remarks length of the *one-letter* aa must equal the length of ss
	SegmentInsert(
		Interval const & i,
		String const & ss,
		String const & aa,
		Pose const & insert,
		bool const keep_known_bb_torsions_at_junctions = false,
		SegmentInsertConnectionScheme::Enum connection_scheme = protocols::forge::build::SegmentInsertConnectionScheme::RANDOM_SIDE
	);


	/// @brief copy constructor
	SegmentInsert( SegmentInsert const & rval );


	/// @brief default destructor
	virtual
	~SegmentInsert();


public: // assignment


	/// @brief copy assignment
	SegmentInsert & operator =( SegmentInsert const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	BuildInstructionOP clone() const;


public: // accessors


	/// @brief the character representing the insertion point
	/// @return '^'
	inline
	static
	char insertion_char() {
		return '^';
	}


	/// @brief performing a pure insertion/extension?
	/// @return True if SegmentInsert was initialized with an interval whose
	///  left endpoint == 0, False otherwise.
	inline
	bool performing_pure_insertion() const {
		return original_interval().left == 0;
	}


	/// @brief Attempt to keep the phi at original_interval().left and the psi+omega
	///  at original_interval().right present from the original Pose in the modified
	///  Pose?  Not applicable for pure insertions.  Default False.
	inline
	bool keep_known_bb_torsions_at_junctions() const {
		return keep_known_bb_torsions_at_junctions_;
	}


	/// @brief the pose to insert
	Pose const & insert_pose() const;


	/// @brief get secondary structure string (includes insertion point as '^')
	inline
	String const & ss() const {
		return ss_;
	}


	/// @brief get annotated amino acid string (includes insertion point as '^')
	inline
	String const & aa() const {
		return aa_;
	}


	/// @brief get the insert connection scheme -- whether the insertion will
	///  be glued on the N-side, C-side, or if SegmentInsert will pick randomly
	///  between the two
	inline
	SegmentInsertConnectionScheme::Enum insert_connection_scheme() const {
		return insert_connection_scheme_;
	}


	/// @brief get the absolute index of the insertion point with respect to the
	///  flanking regions (i.e. the index inside the ss string)
	Size insertion_point_absolute_index() const;


	/// @brief get the residue at the start of the insertion relative to the
	///  modified interval (flanking positions are not part of the insertion!)
	Size insertion_start_residue() const;


	/// @brief get the residue at the end of the insertion relative to the
	///  modified interval (flanking positions are not part of the insertion!)
	Size insertion_end_residue() const;


	/// @brief get the number of flanking residues to the left of the insertion
	///  point
	Size flanking_left_nres() const;


	/// @brief get the number of flanking residues to the right of the insertion
	///  point
	Size flanking_right_nres() const;


	/// @brief get the ss string of the flanking residues to the left of the
	///  insertion point
	String flanking_left_ss() const;


	/// @brief get the ss string of the flanking residues to the right of the
	///  insertion point
	String flanking_right_ss() const;


	/// @brief get the annotated aa string of the flanking residues to the left
	///  of the insertion point
	String flanking_left_aa() const;


	/// @brief get the annotated aa string of the flanking residues to the right
	///  of the insertion point
	String flanking_right_aa() const;


	/// @brief a torsion (bb/chi) specific override movemap indexed wrt the insert Pose
	///  (residue indices may only be within the range [1, insert_pose.n_residue()]
	/// @remarks When generating the movemap(), this torsion movemap will be enforced.
	///  Only *explicit* settings of TorsionType, MoveMapTorsionID, and TorsionID will
	///  be honored.  Implicit false settings are ignored.
	inline
	MoveMap const & insert_pose_torsion_override_movemap() const {
		return insert_pose_torsion_override_movemap_;
	}


public: // virtual accessors


	/// @brief is the original interval storing valid information, or is empty
	///  or being used for something else?
	/// @return always false, original interval invalid
	inline
	virtual
	bool original_interval_valid() const {
		return false;
	}


	/// @brief a copy of the working range of residues specifying the modified region
	/// @details This residue range can change wrt length changes in Pose /Conformation
	///  being watched.
	virtual
	Interval interval() const;


	/// @brief return a copy of the set of positions within the new region
	///  that were pre-existing in the original Pose prior to modify()
	/// @return An empty set, there are no pre-existing positions.
	virtual
	Positions preexisting_positions() const;


	/// @brief return a copy of the set of positions that are "new" and did
	///  not exist in the original Pose.
	virtual
	Positions new_positions() const;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has a defined conformation.  E.g. existing or copied residues.
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	Positions defined_positions() const;


	/// @brief return a copy of the set of positions within the newly modified
	///  region that has an undefined conformation.  E.g. newly created residues.
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	Positions undefined_positions() const;


	/// @brief return a copy of the MoveMap that defines the moveable/fixed
	///  positions/dofs for this instruction
	/// @return TO BE FILLED IN
	/// @details This set can change wrt length changes in Pose/Conformation being
	///  watched.
	virtual
	MoveMap movemap() const;


public: // mutators


	/// @brief set a torsion (bb/chi) specific override movemap indexed wrt the insert
	///  Pose (residue indices may only be within the range [1, insert_pose.n_residue()]
	/// @remarks When generating the movemap(), this torsion movemap will be enforced.
	///  Only *explicit* settings of TorsionType, MoveMapTorsionID, and TorsionID will
	///  be honored.  Implicit false settings are ignored.
	void insert_pose_torsion_override_movemap( MoveMap const & mm );


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


private: // initialization


	/// @brief init to be called during non-default constructors
	void init();


private: // data


	/// @brief the working range of residues
	/// @remarks this range can shift if listening to a Pose/Conformation and the number
	///  of residues changes
	Interval interval_;


	/// @brief secondary structure string defining the flanking regions and the
	///  insertion point
	/// @details The insertion point is marked as '^'.  For instance, "L^LLL"
	///  means one residue to the left of the insertion point and three residues
	///  to the right of the insertion point.
	String ss_;


	/// @brief annotated amino acid string including the insertion point
	/// @details Length of the one-letter version must be equal to length of ss.
	/// The insertion point must exist in the same spot as the insertion point in
	///  secondary structure string.
	String aa_;


	/// @brief Attempt to keep the omega at original_interval().left-1, the phi
	///  at original_interval().left and the psi+omega at original_interval().right
	///  present from the original Pose in the modified Pose?  This should be
	///  false for pure insertions.
	/// @details If True, during modify(), will (1) set the omega of interval_.left-1
	///  in the newly modified Pose equal to the omega of the original Pose at
	///  original_interval().left-1, (2) set the phi of interval_.left in
	///  the newly modified Pose equal to the phi of the original Pose in
	///  original_interval().right, and (3) set the psi+omega of interval_.right
	///  in the newly modified Pose equal to the original psi+omega of
	///  original_interval().right.  Default False.
	bool keep_known_bb_torsions_at_junctions_;


	/// @brief connect insertion on its N-side, C-side, or decide randomly
	///  between the two (default RANDOM_SIDE)
	SegmentInsertConnectionScheme::Enum insert_connection_scheme_;


	/// @brief insert this Pose
	Pose insert_pose_;


	/// @brief a torsion (bb/chi) specific override movemap indexed wrt the insert Pose
	///  (residue indices may only be within the range [1, insert_pose.n_residue()]
	/// @remarks When generating the movemap(), this torsion movemap will be enforced.
	///  Only *explicit* settings of TorsionType, MoveMapTorsionID, and TorsionID will
	///  be honored.  Implicit false settings are ignored.
	MoveMap insert_pose_torsion_override_movemap_;


};


} // namespace build
} // namespace forge
} // namespace protocols


#endif /* INCLUDED_protocols_forge_build_SegmentInsert_HH */
