// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/symmetry/SymmetryTransform.cc
/// @brief  Stores information about a transform between two symmetry subunits.
/// @details Currently stores both a HomogeneousTransform for the translation and rotation, a boolean for whether
/// this transform involves a mirror operation, and two vectors defining the mirror transformation plane.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit Headers
#include <core/conformation/symmetry/SymmetryTransform.hh>

// Package headers
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/id/types.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/RT.hh>

#include <core/conformation/Conformation.hh>
#include <numeric/xyzMatrix.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.functions.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

// C++ headers
#include <iostream>

// core utilities
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/HomogeneousTransform.srlz.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "core.conformation.SymmetryTransform" );

namespace core {
namespace conformation {
namespace symmetry {

/// @brief Default constructor.
///
SymmetryTransform::SymmetryTransform():
	htransform_(),
	mirror_z_(false)
	//TODO -- initialize private member variables here
{}

/// @brief Data constructor.
/// @param[in] p1,p2,p3 xyzVectors for initializing the HomogeneousTransform.
SymmetryTransform::SymmetryTransform(
	numeric::xyzVector< core::Real > const & p1,
	numeric::xyzVector< core::Real > const & p2,
	numeric::xyzVector< core::Real > const & p3,
	bool const mirror_z
) :
	htransform_( p1, p2, p3 ),
	mirror_z_(false)
{
	set_mirror_z(mirror_z);
}

/// @brief Destructor
SymmetryTransform::~SymmetryTransform()
{}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the copy.
SymmetryTransformOP
SymmetryTransform::clone() const {
	return SymmetryTransformOP( new SymmetryTransform(*this) );
}

/// @brief Set whether this HomogeneousTransform mirrors about the z-axis or not.
///
void
SymmetryTransform::set_mirror_z( bool const val ) {
	if ( mirror_z_ == val ) return; //Do nothing if we've already set the value this way.
	else {
		mirror_z_ = val;
		htransform_ = numeric::HomogeneousTransform<core::Real>(
			htransform_.rotation_matrix()*kinematics::Jump::mirror_z_transform,
			htransform_.point() );
	}
	return;
}


} // symmetry
} // conformation
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::symmetry::SymmetryTransform::save( Archive & arc ) const {
	arc( CEREAL_NVP( htransform_ ) ); // numeric::HomogeneousTransform< core::Real >
	arc( CEREAL_NVP( mirror_z_ ) ); // _Bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::symmetry::SymmetryTransform::load( Archive & arc ) {
	arc( htransform_ ); // numeric::HomogeneousTransform< core::Real >
	arc( mirror_z_ ); // _Bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::symmetry::SymmetryTransform );
CEREAL_REGISTER_TYPE( core::conformation::symmetry::SymmetryTransform )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_symmetry_SymmetryTransform )
#endif // SERIALIZATION
