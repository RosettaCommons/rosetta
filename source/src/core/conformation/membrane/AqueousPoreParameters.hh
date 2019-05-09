// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/membrane/AqueousPoreParameters.hh
/// @brief A class for defining an aqueous pore
/// @author Rebecca Alford (rfalford12@gmail.com)


#ifndef INCLUDED_core_conformation_membrane_AqueousPoreParameters_hh
#define INCLUDED_core_conformation_membrane_AqueousPoreParameters_hh

// Unit Headers
#include <core/conformation/membrane/AqueousPoreParameters.fwd.hh>

// Numeric Headers
#include <numeric/MathMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>
#include <numeric/interpolation/spline/SimpleInterpolator.hh>
#include <numeric/cubic_polynomial.hh>
#include <core/types.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <cstdlib>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace conformation {
namespace membrane {

/// @brief A class for defining an aqueous pore
class AqueousPoreParameters : public utility::pointer::ReferenceCount {

	typedef utility::vector1< numeric::CubicPolynomial > piecewise_poly;

public:

	// Default construction
	AqueousPoreParameters();

	/// @brief Construct an AqueousPoreParameters object from scratch
	AqueousPoreParameters(
		core::Real const min_center_x,
		core::Real const max_center_x,
		core::Real const min_center_y,
		core::Real const max_center_y,
		core::Real const min_major_radius,
		core::Real const max_major_radius,
		core::Real const min_minor_radius,
		core::Real const max_minor_radius,
		core::Real const min_rotation_angle,
		core::Real const max_rotation_angle,
		utility::vector1< core::Real > boundaries,
		piecewise_poly pore_center_x,
		piecewise_poly pore_center_y,
		piecewise_poly pore_major_radius,
		piecewise_poly pore_minor_radius,
		piecewise_poly pore_rotation_angle
	);

	AqueousPoreParameters(AqueousPoreParameters const & src);
	virtual ~AqueousPoreParameters();

	AqueousPoreParametersOP
	clone() const;

public: // public accessor methods

	core::Real boundaries( core::Size index ) const;

	// plain functions
	core::Real pore_center_x( core::Real const zcoord ) const;
	core::Real pore_center_y( core::Real const zcoord ) const;
	core::Real pore_major_radius( core::Real const zcoord ) const;
	core::Real pore_minor_radius( core::Real const zcoord ) const;
	numeric::MathMatrix< core::Real > pore_rotation( core::Real const zcoord ) const;

	// derivatives - temporary rotation angle func until I derive the full form
	core::Real pore_center_x_deriv( core::Real const zcoord ) const;
	core::Real pore_center_y_deriv( core::Real const zcoord ) const;
	core::Real pore_major_radius_deriv( core::Real const zcoord ) const;
	core::Real pore_minor_radius_deriv( core::Real const zcoord ) const;
	core::Real pore_rotation_deriv( core::Real const zcood ) const;

private: // helper methods

	core::Real
	eval_piecewise_cubic_polynomial(
		utility::vector1< numeric::CubicPolynomial > piecewise_poly,
		core::Real const min_score,
		core::Real const max_score,
		core::Real const zcoord
	) const;

	core::Real
	eval_piecewise_cubic_polynomial_deriv(
		utility::vector1< numeric::CubicPolynomial > piecewise_poly,
		core::Real const zcoord
	) const;

private:

	core::Real min_center_x_;
	core::Real max_center_x_;
	core::Real min_center_y_;
	core::Real max_center_y_;

	core::Real min_major_radius_;
	core::Real max_major_radius_;
	core::Real min_minor_radius_;
	core::Real max_minor_radius_;

	core::Real min_rotation_angle_;
	core::Real max_rotation_angle_;

	utility::vector1< core::Real > boundaries_;

	piecewise_poly pore_center_x_;
	piecewise_poly pore_center_y_;
	piecewise_poly pore_major_radius_;
	piecewise_poly pore_minor_radius_;
	piecewise_poly pore_rotation_angle_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //core
} //conformation
} //membrane

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_conformation_membrane_AqueousPoreParameters )
#endif // SERIALIZATION


#endif //INCLUDED_core_conformation_membrane_AqueousPoreParameters_hh





