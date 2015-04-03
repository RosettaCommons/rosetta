// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/match_enzdes_util/LigandConformer.hh
/// @brief  Declaration of a class to hold a ligand conformation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_toolbox_match_enzdes_util_LigandConformer_hh
#define INCLUDED_protocols_toolbox_match_enzdes_util_LigandConformer_hh

// Unit headers
#include <protocols/toolbox/match_enzdes_util/LigandConformer.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// ObjexxFCL headers
#include <ObjexxFCL/FArray2D.fwd.hh>

// Numeric headers
#include <numeric/HomogeneousTransform.hh>
#include <numeric/geometry/hashing/SixDHasher.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>

#include <core/id/AtomID.fwd.hh>
//#include <protocols/match/Hit.fwd.hh>
#include <utility/fixedsizearray1.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace match_enzdes_util {

class LigandConformer : public utility::pointer::ReferenceCount {
public:
	typedef utility::pointer::ReferenceCount      parent;
	typedef core::Size                              Size;
	typedef core::Real                              Real;
	typedef core::Vector                          Vector;
	typedef numeric::geometry::hashing::Real6      Real6;
	typedef numeric::HomogeneousTransform< Real > HTReal;

public:
	LigandConformer();

	LigandConformer( LigandConformer const & );

	virtual ~LigandConformer();

	Real
	atom1_atom2_distance() const;

	Real
	atom2_atom3_distance() const;

	/// @brief Returns an angle in degrees between the three downstream atoms.
	Real
	atom1_atom2_atom3_angle() const;

	/// @brief returns the distance between orientation atom 1 and orientation atom 2
	Real
	oatom1_oatom2_distance() const;

	/// @brief returns the distance between orientation atom 2 and orientation atom 3
	Real
	oatom2_oatom3_distance() const;

	/// @brief Returns an angle in degrees between the three orientation atoms.
	Real
	oatom1_oatom2_oatom3_angle() const;


	void
	coordinates_from_orientation(
		Real6 const & orientation,
		utility::vector1< core::id::AtomID > const & atom_indices,
		utility::vector1< Vector > & atom_coords
	) const;

	/// @brief Specify the residue, with coordinates, that's being used as the downstream
	/// partner.  This class is meant to be used in conjuction with the ClassicMatchAglrotihm,
	/// and therefore the initialization routines are specific for that algorithm.  In this
	/// initialization function, one must list atoms "D1, D2 and D3" in the convention of
	/// describing the rigid-body orientation between three atoms of the upstream partner
	/// (atoms U3, U2 & U1) and three atoms of the downstream partner (atoms D1, D2 & D3) in terms
	/// of 2 angles, 1 distance, and 3 dihedrals.  The user must also list the 3 atoms used to
	/// define the orientation frame of the downstream ligand.  It is essential to the
	/// matching algorithm that the same three orientation atoms are used for all LigandConformers.
	void
	initialize_from_residue(
		Size D1,
		Size D2,
		Size D3,
		Size orientation_atom1,
		Size orientation_atom2,
		Size orientation_atom3,
		core::conformation::Residue const & residue
	);

	void ignore_h_collisions( bool setting );

	/// @brief The orientaton frame at orientation atom 3 given
	/// the coordinate frame at D3 (this frame is called frame3)
	Real6
	global_orientation_from_frame3(
		HTReal const & frame3
	) const;

	/// @brief The orientation frame at orientation atom 3 given
	/// orientation atom 3's xyz coordinates and the euler angles
	/// describing the frame
	HTReal
	frame_from_global_orientation(
		Real6 orientation
	) const;

	void
	move_atoms_to_collcheck_begin( utility::vector1< Size > const & restype_atnos_to_move_early );

	inline
	Size
	n_collision_check_atoms() const {
		return collision_check_id_2_restype_id_.size();
	}

	inline
	Size
	restype_id_2_collision_check_id( Size restype_atomno ) const {
		return restype_id_2_collision_check_id_[ restype_atomno ];
	}

	inline
	Size
	collision_check_id_2_restype_id( Size coll_check_id ) const {
		return collision_check_id_2_restype_id_[ coll_check_id ];
	}

	inline
	Vector
	coordinate_in_D3_frame( Size restype_atomno, HTReal const & frame3 ) const {
		return frame3 * points_in_D3_frame_[ restype_atomno ];
	}

	inline
	Vector
	coordinate_in_global_frame( Size restype_atomno, HTReal const & orientation_frame ) const {
		return orientation_frame * points_in_global_orintation_frame_[ restype_atomno ];
	}

	/// @ brief helper function to get the coordinates in 2D FArray format
	void
	get_global_coords_as_FArray2D(
		ObjexxFCL::FArray2D< numeric::Real > & coords,
		HTReal const & orientation_frame,
		utility::vector1< core::Size > const & restype_atomnos
	) const;

	core::chemical::ResidueType
	get_lig_restype() const;

private:

	void
	create_collcheck_ordering( utility::vector1< bool > selected, Size count_from );

private:

	core::chemical::ResidueTypeCOP ligand_restype_;

	/// The indices of the three atoms defining the orientation of the
	/// ligand in the global coordinate frame
	/// These indices are in the restype indexing of atoms.
	utility::fixedsizearray1< Size, 3 >         orientation_atoms_;
	utility::fixedsizearray1< Vector, 3 >       oats_in_D3_frame_;
	HTReal                                      oframe_in_D3frame_;

	/// The coordinates of all the ligand atoms in the global orientation frame.
	utility::vector1< Vector >                  points_in_global_orintation_frame_;

	/// The indices for the three atoms defining the location of the downstream partner
	/// from the upstream partner.  D1 D2 and D3.  These indices are in the restype indexing of atoms.
	utility::fixedsizearray1< Size, 3 >         atoms_123_;
	Real d12_; // distance from D1 to D2
	Real d23_; // distance drom D2 to D3
	Real ang123_; /// angle between D1, D2 and D3

	/// The coordinates of the other ligand atoms in the coordinate frame from atom D3.
	utility::vector1< Vector >       points_in_D3_frame_;
	bool ignore_h_collisions_;
	utility::vector1< Size >         collision_check_id_2_restype_id_;
	utility::vector1< Size >         restype_id_2_collision_check_id_;

};

}
}
}

#endif
