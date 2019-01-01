// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/util.hh
/// @brief   utility functions for core/scoring/nmr classes
/// @details last Modified: 05/26/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_util_HH
#define INCLUDED_core_scoring_nmr_util_HH

// Package headers
#include <core/io/nmr/AtomSelection.fwd.hh>
#include <core/scoring/nmr/types.hh>

// Project headers
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace core {
namespace scoring {
namespace nmr {

NMR_VALUE_AVERAGING_TYPE
convert_string_to_averaging_type(std::string const & averaging_type);

SINGLE_NMR_VALUE_WEIGHTING
convert_string_to_weighting_scheme(std::string const & weighting_scheme);

/// @brief determine the RDC type from the spins' atom names
RDC_TYPE
rdc_type_from_atom_names(std::pair< core::io::nmr::AtomSelection,core::io::nmr::AtomSelection > const & spinsAB);

RDC_NORM_TYPE
convert_string_to_normalization_type(std::string const & normalization_type);

PRE_SPIN_TYPE
pre_spin_type_from_atom_name(core::io::nmr::AtomSelection const & atom);

PRE_RATE_TYPE
convert_string_to_rate_type(std::string const & rate_type);

std::string
convert_rdc_type_to_string(RDC_TYPE const & type);

std::string
convert_pre_spin_type_to_string(PRE_SPIN_TYPE const & type);

/// @brief function that handles conversion from pseudoatom identifier for degenerate protons to fullatom name
utility::vector1<id::AtomID>
lookup_pseudoprotons(
	Size const rsd,
	std::string const & atomname,
	pose::Pose const & pose
);

/// @brief creates a rotation matrix given three Euler angles
///        and a convention as determined at runtime
Matrix
rotation_matrix_from_euler_angles(
	Vector const & angles,
	EULER_CONVENTION convention
);

/// @brief determines the three Euler angles given a rotation matrix
///        and the respective Euler convention as determined at runtime
Vector
euler_angles_from_rotation_matrix(
	Matrix const & rotM,
	EULER_CONVENTION convention
);

/// @brief orders the three eigenvalues of NMRTensor and brings NMRTensor
///        in unique tensor representation
void
order_tensor_parameters(
	utility::vector1<Real> & params,
	EULER_CONVENTION convention
);

/// @brief Utility function that transforms an NMR tensor from the frame of residue_from
///        into the coordinate frame of the symmetric residue_to given a symmetric pose.
///        Returns the original matrix if the pose is not symmetric or if residue_from and
///        residue_to belong to the same subunit residue_to.
Matrix
apply_tensor_transformation(
	pose::Pose & pose,
	Matrix const & Min,
	Size resid_from,
	Size resid_to
);

/// @brief Utility function that transforms a vector of metal (or spinlabel) coordinates from
///        the frame of residue_from into the coordinate frame of the symmetric residue_to given a symmetric pose.
///        Returns the original vector if the pose is not symmetric or if residue_from and
///        residue_to belong to the same subunit residue_to.
Vector
apply_vector_transformation(
	pose::Pose & pose,
	Vector const & Vin,
	Size resid_from,
	Size resid_to
);

/// @brief Utility function that rotates a vector from the frame of residue_from
///        into the coordinate frame of the symmetric residue_to given a symmetric pose.
///        Returns the original vector if the pose is not symmetric or if residue_from and
///        residue_to belong to the same subunit residue_to.
Vector
apply_vector_rotation(
	pose::Pose & pose,
	Vector const & Vin,
	Size resid_from,
	Size resid_to
);

/// @brief Create a local coordinate frame from backbone coordinates.
///        The local frame is defined by a homogeneous transform object
///        which can be created from three axes and one point.
///        The center is located at CA. The z'-axis is along CA - CB.
///        The x'-axis is within the CA-CB - CA-C plane and the
///        y'-axis is perpendicular to the CB - CA - C plane.
HT
define_local_frame(
	Vector const & CA,
	Vector const & CB,
	Vector const & C
);

/// @brief Auxiliary PCS function
/// @params par: Tensor values [xM, yM, zM, Xax, Xrh]
/// rotM: Rotation matrix to transform spin coordinates in tensor frame
/// spin_coord: Spin xyz coordinates
/// scal: optional scaling factor
Real
pcs_func(
	Vec5 const & par,
	Matrix const & rotM,
	Vector const & spin_coord,
	Real const & scal=1.0
);

} // nmr
} // scoring
} // core

#endif // INCLUDED_core_scoring_nmr_util_HH
