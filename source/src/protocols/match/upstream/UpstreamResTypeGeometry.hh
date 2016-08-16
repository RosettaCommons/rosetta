// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/upstream/UpstreamResTypeGeometry.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_upstream_UpstreamResTypeGeometry_hh
#define INCLUDED_protocols_match_upstream_UpstreamResTypeGeometry_hh

// Unit headers
#include <protocols/match/upstream/UpstreamResTypeGeometry.fwd.hh>

// Project headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/types.hh>

// Numeric headers
#include <numeric/HomogeneousTransform.hh>
#include <numeric/xyzVector.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace upstream {

/// @brief A simple class that describes the geometry for a particular
/// residue type.  It describes the coordinate frame geometry for the
/// fourth atom defining each chi dihedral. The fourth atom is called
/// the "chi tip" atom, as it's at the tip of the growing kinematic chain
/// when building chi i. This class also describes the location
/// of the atoms controlled by each chi which are not the chi-tip atoms;
/// it measures their location in the coordinate frame of the chi-tip atom.
///
/// @details To generate the coordinate of the chi-tip atom, the
/// stored coordinate frame is multiplied
/// by the coordinate frame at the third atom after that coordinate frame
/// has been multipled by the chi-angle-z-axis rotation HT.
class UpstreamResTypeGeometry : public utility::pointer::ReferenceCount
{
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~UpstreamResTypeGeometry();
	typedef core::Size                            Size;
	typedef core::Real                            Real;
	typedef core::Vector                          Vector;
	typedef numeric::HomogeneousTransform< Real > HTReal;

public:
	UpstreamResTypeGeometry();
	UpstreamResTypeGeometry( core::chemical::ResidueType const & );

	void initialize_from_residue_type( core::chemical::ResidueType const & );

public:

	/// @brief the name of the residue type used to generate this geometry
	std::string const & name() const {
		return restype_name_;
	}

	/// @brief the number of atoms in this residue type
	Size natoms() const {
		return controlling_chi_for_atom_.size();
	}

	Size nchi() const {
		return chitip_atoms_.size();
	}

	bool atom_controlled_by_any_chi( Size atomno ) const {
		return controlling_chi_for_atom_[ atomno ] != 0;
	}

	bool atom_is_chitip( Size atomno ) const {
		return controlling_chi_for_atom_[ atomno ] != 0 && which_point_for_atom_[ atomno ] == 0;
	}

	utility::vector1< Size > const &
	controlling_chi_for_atom() const { return controlling_chi_for_atom_; }

	utility::vector1< Size > const &
	which_point_for_atom() const { return which_point_for_atom_; }

	utility::vector1< Size > const &
	chitip_atoms() const { return chitip_atoms_; }

	Size
	chitip_atom( Size chi ) const {
		return chitip_atoms_[ chi ];
	}

	HTReal const &
	pre_chitip_transform( Size chi ) const {
		return pre_chitip_transforms_[ chi ];
	}

	utility::vector1< HTReal > const &
	ht_for_chitip_atoms() const { return ht_for_chitip_atoms_; }

	HTReal const &
	ht_for_chitip_atom( Size chi ) const {
		return ht_for_chitip_atoms_[ chi ];
	}

	Size
	n_nonchitip_atoms_for_chi( Size chi ) const {
		return nonchitip_atoms_[ chi ].size();
	}

	utility::vector1< utility::vector1< Size > > const &
	nonchitip_atoms() const { return nonchitip_atoms_; }

	Size
	nonchitip_atom( Size chi, Size which_nonchitip_atom_for_chi ) const {
		return nonchitip_atoms_[ chi ][ which_nonchitip_atom_for_chi ];
	}

	utility::vector1< utility::vector1< Vector > > const &
	points_for_nonchitip_atoms() const { return points_for_nonchitip_atoms_; }

	utility::vector1< Vector > const &
	points_for_nonchitip_atoms( Size chi ) const {
		return points_for_nonchitip_atoms_[ chi ];
	}

	/// @brief Convenience function: get the coordinate in the chitip frame
	/// for a particular atom.  The atom must be a non-chitip atom that is
	/// not part of the backbone (it must be controlled by a chi angle).
	Vector const &
	point_for_nonchitip_atom( Size atom ) {
		assert( atom_controlled_by_any_chi( atom ) && !atom_is_chitip( atom ) );
		return points_for_nonchitip_atoms_[ controlling_chi_for_atom_[ atom ] ]
			[ which_point_for_atom_[ atom ] ];
	}


	Size N_atom_id() const { return N_atom_id_; }
	Size CA_atom_id() const { return CA_atom_id_; }
	Size C_atom_id() const { return C_atom_id_; }
	Size O_atom_id() const { return O_atom_id_; }
	Size CB_atom_id() const { return CB_atom_id_; }
	Size H_atom_id() const  { return H_atom_id_; }
	Size HA_atom_id() const  { return HA_atom_id_; }


	bool has_N_atom() const { return N_atom_id_ != 0; }
	bool has_CA_atom() const { return CA_atom_id_ != 0; }
	bool has_C_atom() const { return C_atom_id_ != 0; }
	bool has_O_atom() const { return O_atom_id_ != 0; }
	bool has_CB_atom() const { return CB_atom_id_ != 0; }
	bool has_H_atom() const  { return H_atom_id_ != 0; }
	bool has_HA_atom() const  { return HA_atom_id_ != 0; }

	bool
	atom_has_nonchi_coordinate( Size restype_atomid ) const;

	Vector const &
	coordinate_for_nonchi_atom_in_ideal_frame( Size restype_atomid ) const;


private:
	/// Data

	std::string restype_name_;

	utility::vector1< Size >   controlling_chi_for_atom_;
	utility::vector1< Size >   which_point_for_atom_;

	utility::vector1< Size >   chitip_atoms_;
	utility::vector1< HTReal > pre_chitip_transforms_;
	utility::vector1< HTReal > ht_for_chitip_atoms_;


	utility::vector1< utility::vector1< Size > >   nonchitip_atoms_;
	utility::vector1< utility::vector1< Vector > > points_for_nonchitip_atoms_;

	/// The ideal frame is defined at Calpha with the half point between N and C
	/// in the N plane-halfpoint-Calpha plane.
	/// i.g. HTReal( N, halfpoint, CAlpha);
	/// Non-chi dependendent atoms are measured from the ideal coordinate in this frame.
	/// This includes Cbeta and the Halphas.  It also includes O and H, but since their
	/// geometry depends on phi and psi, this data would be inappropriate for them.
	utility::vector1< Size   > restype_atom_id_2_nonchi_atom_id_;
	utility::vector1< Size   > nonchi_atom_id_2_restype_atom_id_;
	utility::vector1< Vector > nonchi_atoms_in_ideal_frame_;

	Size N_atom_id_;
	Size CA_atom_id_;
	Size C_atom_id_;
	Size O_atom_id_;
	Size CB_atom_id_;
	Size H_atom_id_;
	Size HA_atom_id_;

};


}
}
}

#endif

