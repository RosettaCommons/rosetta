// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/symmetry/SymmetryTransform.cc
/// @brief  Headers for the SymmetryTransform class, which stores information about a
/// transform between two symmetry subunits.
/// @details Currently stores both a HomogeneousTransform for the translation and rotation, a boolean for whether
/// this transform involves a mirror operation, and two vectors defining the mirror transformation plane.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_core_conformation_symmetry_SymmetryTransform_hh
#define INCLUDED_core_conformation_symmetry_SymmetryTransform_hh

#include <core/conformation/symmetry/SymmetryTransform.fwd.hh>

// Unit Headers
#include <core/conformation/symmetry/SymmData.fwd.hh>

//core
#include <core/types.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>

//numeric
#include <numeric/HomogeneousTransform.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>
#include <iosfwd>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace conformation {
namespace symmetry {


//Symm_info

class SymmetryTransform : public utility::pointer::ReferenceCount {

public: //Constructor, destructor, copy, clone:

	/// @brief Default constructor.
	///
	SymmetryTransform();

	/// @brief Data constructor.
	/// @param[in] p1,p2,p3 xyzVectors for initializing the HomogeneousTransform.
	/// @param[in] mirror_z Should the copy be mirrored about the z-axis?
	SymmetryTransform(
		numeric::xyzVector< core::Real > const & p1,
		numeric::xyzVector< core::Real > const & p2,
		numeric::xyzVector< core::Real > const & p3,
		bool const mirror_z=false
	);

	/// @brief Default destructor.
	///
	virtual ~SymmetryTransform();

	/// @brief Clone operator.
	/// @details Copy this object and return an owning pointer to the copy.
	SymmetryTransformOP clone() const;

public: //Public methods:

	/// @brief Get the HomogeneousTransform matrix.
	///
	inline numeric::HomogeneousTransform< core::Real > const &
	ht() const {
		return htransform_;
	}

	/// @brief Set whether this HomogeneousTransform mirrors about the z-axis or not.
	///
	void set_mirror_z( bool const val );

	/// @brief Get whether this HomogeneousTransform mirrors aout the z-axis or not.
	///
	bool mirror_z() const { return mirror_z_; }


private: //Private member variables:

	/// @brief The HomogeneousTransform object that describes the translational and rotational part
	/// of the symmetry transform.
	numeric::HomogeneousTransform< core::Real > htransform_;

	/// @brief Should this HomogeneousTransform mirror the object about the z-axis?
	/// @details False by default.
	bool mirror_z_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // SymmetryTransform class

} // symmetry
} // conformation
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_symmetry_SymmetryTransform )
#endif // SERIALIZATION

#endif
