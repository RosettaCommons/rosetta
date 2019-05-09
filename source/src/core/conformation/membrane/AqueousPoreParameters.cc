// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/conformation/membrane/AqueousPoreParameters.cc
/// @brief A class for defining an aqueous pore
/// @author Rebecca Alford (rfalford12@gmail.com)

// Utility Headers
#include <core/conformation/membrane/AqueousPoreParameters.hh>

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

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "core.conformation.membrane.AqueousPoreParameters" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace conformation {
namespace membrane {

typedef utility::vector1< numeric::CubicPolynomial > piecewise_poly;

AqueousPoreParameters::AqueousPoreParameters() :
	utility::pointer::ReferenceCount(),
	min_center_x_(0),
	max_center_x_(0),
	min_center_y_(0),
	max_center_y_(0),
	min_major_radius_(0),
	max_major_radius_(0),
	min_minor_radius_(0),
	max_minor_radius_(0),
	min_rotation_angle_(0),
	max_rotation_angle_(0),
	boundaries_(),
	pore_center_x_(0),
	pore_center_y_(0),
	pore_major_radius_(),
	pore_minor_radius_(),
	pore_rotation_angle_()
{}

AqueousPoreParameters::AqueousPoreParameters(
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
) : utility::pointer::ReferenceCount(),
	min_center_x_( min_center_x ),
	max_center_x_( max_center_x ),
	min_center_y_( min_center_y ),
	max_center_y_( max_center_y ),
	min_major_radius_( min_major_radius ),
	max_major_radius_( max_major_radius ),
	min_minor_radius_( min_minor_radius ),
	max_minor_radius_( max_minor_radius ),
	min_rotation_angle_( min_rotation_angle ),
	max_rotation_angle_( max_rotation_angle ),
	boundaries_( boundaries ),
	pore_center_x_( pore_center_x ),
	pore_center_y_( pore_center_y ),
	pore_major_radius_( pore_major_radius ),
	pore_minor_radius_( pore_minor_radius ),
	pore_rotation_angle_( pore_rotation_angle )
{}

AqueousPoreParameters::~AqueousPoreParameters() {}

AqueousPoreParameters::AqueousPoreParameters( AqueousPoreParameters const & src )
: utility::pointer::ReferenceCount( src ),
	min_center_x_( src.min_center_x_ ),
	max_center_x_( src.max_center_x_ ),
	min_center_y_( src.min_center_y_ ),
	max_center_y_( src.max_center_y_ ),
	min_major_radius_( src.min_major_radius_ ),
	max_major_radius_( src.max_major_radius_ ),
	min_minor_radius_( src.min_minor_radius_ ),
	max_minor_radius_( src.max_minor_radius_ ),
	min_rotation_angle_( src.min_rotation_angle_ ),
	max_rotation_angle_( src.max_rotation_angle_ ),
	boundaries_( src.boundaries_ ),
	pore_center_x_( src.pore_center_x_ ),
	pore_center_y_( src.pore_center_y_ ),
	pore_major_radius_( src.pore_major_radius_ ),
	pore_minor_radius_( src.pore_minor_radius_ ),
	pore_rotation_angle_( src.pore_rotation_angle_ )
{}

AqueousPoreParametersOP
AqueousPoreParameters::clone() const {
	return AqueousPoreParametersOP( new AqueousPoreParameters( *this ) );
}

core::Real
AqueousPoreParameters::boundaries( core::Size index ) const {
	return boundaries_[ index ];
}

core::Real
AqueousPoreParameters::pore_center_x( core::Real const zcoord ) const {
	return eval_piecewise_cubic_polynomial( pore_center_x_, min_center_x_, max_center_x_, zcoord );
}

core::Real
AqueousPoreParameters::pore_center_y( core::Real const zcoord ) const {
	return eval_piecewise_cubic_polynomial( pore_center_y_, min_center_y_, max_center_y_, zcoord );
}

core::Real
AqueousPoreParameters::pore_major_radius( core::Real const zcoord ) const {
	return eval_piecewise_cubic_polynomial( pore_major_radius_, min_major_radius_, max_major_radius_, zcoord );
}

core::Real
AqueousPoreParameters::pore_minor_radius( core::Real const zcoord ) const {
	return eval_piecewise_cubic_polynomial( pore_minor_radius_, min_minor_radius_, max_minor_radius_, zcoord );
}

numeric::MathMatrix< core::Real >
AqueousPoreParameters::pore_rotation( core::Real const zcoord ) const {

	using namespace numeric;

	core::Real theta( eval_piecewise_cubic_polynomial( pore_rotation_angle_, min_rotation_angle_, max_rotation_angle_, zcoord ) );
	core::Real cos_theta( std::cos( theta ) );
	core::Real sin_theta( std::sin( theta ) );

	MathMatrix< core::Real > rot_matrix( 2, 2 );
	rot_matrix(0,0) = cos_theta;
	rot_matrix(0,1) = -sin_theta;
	rot_matrix(1,0) = sin_theta;
	rot_matrix(1,1) = cos_theta;
	return rot_matrix;
}

core::Real
AqueousPoreParameters::pore_center_x_deriv( core::Real const zcoord ) const {
	return eval_piecewise_cubic_polynomial_deriv( pore_center_x_, zcoord );
}

core::Real
AqueousPoreParameters::pore_center_y_deriv( core::Real const zcoord ) const {
	return eval_piecewise_cubic_polynomial_deriv( pore_center_y_, zcoord );
}

core::Real
AqueousPoreParameters::pore_major_radius_deriv( core::Real const zcoord ) const {
	return eval_piecewise_cubic_polynomial_deriv( pore_major_radius_, zcoord );
}

core::Real
AqueousPoreParameters::pore_minor_radius_deriv( core::Real const zcoord ) const {
	return eval_piecewise_cubic_polynomial_deriv( pore_minor_radius_, zcoord );
}

core::Real
AqueousPoreParameters::pore_rotation_deriv( core::Real const zcoord ) const {
	return eval_piecewise_cubic_polynomial_deriv( pore_rotation_angle_, zcoord );
}

core::Real
AqueousPoreParameters::eval_piecewise_cubic_polynomial(
	utility::vector1< numeric::CubicPolynomial > piecewise_poly_func,
	core::Real const min_score,
	core::Real const max_score,
	core::Real const zcoord
) const {

	using namespace numeric;
	using namespace numeric::interpolation::spline;

	if ( zcoord < boundaries_.front() ) {
		return min_score;
	} else if ( zcoord >= boundaries_.back() ) {
		return max_score;
	} else if ( zcoord >= boundaries_.front() && zcoord <= boundaries_.back() ) {
		for ( core::Size ii = 1; ii < boundaries_.size(); ++ii ) {
			if ( boundaries_[ii] <= zcoord && boundaries_[ii+1] > zcoord ) {
				core::Real result = numeric::eval_cubic_polynomial( zcoord, piecewise_poly_func[ii] );
				return result;
			}
		}
	} else {
		utility_exit_with_message( "Polynomial exceeded defined boundaries" );
	}

	// Should never reach this poing
	return numeric::eval_cubic_polynomial( 0, piecewise_poly_func[1] );

}

core::Real
AqueousPoreParameters::eval_piecewise_cubic_polynomial_deriv(
	utility::vector1< numeric::CubicPolynomial > piecewise_poly,
	core::Real const zcoord
) const {

	using namespace numeric;
	using namespace numeric::interpolation::spline;

	if ( zcoord < boundaries_.front() || zcoord >= boundaries_.back() ) {
		return 0.0;
	} else if ( zcoord >= boundaries_.front() && zcoord <= boundaries_.back() ) {
		for ( core::Size ii = 1; ii < boundaries_.size(); ++ii ) {
			if ( boundaries_[ii] <= zcoord && boundaries_[ii+1] > zcoord ) {
				return numeric::cubic_polynomial_deriv( zcoord, piecewise_poly[ii] );
			}
		}
	} else {
		utility_exit_with_message( "Polynomial exceeded defined boundaries" );
	}

	// Should never reach this poing
	return numeric::cubic_polynomial_deriv( 0, piecewise_poly[1] );

}



} //core
} //conformation
} //membrane







#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::membrane::AqueousPoreParameters::save( Archive & arc ) const {
	arc( CEREAL_NVP( min_center_x_ ) ); // core::Real
	arc( CEREAL_NVP( max_center_x_ ) ); // core::Real
	arc( CEREAL_NVP( min_center_y_ ) ); // core::Real
	arc( CEREAL_NVP( max_center_y_ ) ); // core::Real
	arc( CEREAL_NVP( min_major_radius_ ) ); // core::Real
	arc( CEREAL_NVP( max_major_radius_ ) ); // core::Real
	arc( CEREAL_NVP( min_minor_radius_ ) ); // core::Real
	arc( CEREAL_NVP( max_minor_radius_ ) ); // core::Real
	arc( CEREAL_NVP( min_rotation_angle_ ) ); // core::Real
	arc( CEREAL_NVP( max_rotation_angle_ ) ); // core::Real
	arc( CEREAL_NVP( boundaries_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( pore_center_x_ ) ); // piecewise_poly
	arc( CEREAL_NVP( pore_center_y_ ) ); // piecewise_poly
	arc( CEREAL_NVP( pore_major_radius_ ) ); // piecewise_poly
	arc( CEREAL_NVP( pore_minor_radius_ ) ); // piecewise_poly
	arc( CEREAL_NVP( pore_rotation_angle_ ) ); // piecewise_poly
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::membrane::AqueousPoreParameters::load( Archive & arc ) {
	arc( min_center_x_ ); // core::Real
	arc( max_center_x_ ); // core::Real
	arc( min_center_y_ ); // core::Real
	arc( max_center_y_ ); // core::Real
	arc( min_major_radius_ ); // core::Real
	arc( max_major_radius_ ); // core::Real
	arc( min_minor_radius_ ); // core::Real
	arc( max_minor_radius_ ); // core::Real
	arc( min_rotation_angle_ ); // core::Real
	arc( max_rotation_angle_ ); // core::Real
	arc( boundaries_ ); // utility::vector1<core::Real>
	arc( pore_center_x_ ); // piecewise_poly
	arc( pore_center_y_ ); // piecewise_poly
	arc( pore_major_radius_ ); // piecewise_poly
	arc( pore_minor_radius_ ); // piecewise_poly
	arc( pore_rotation_angle_ ); // piecewise_poly
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::membrane::AqueousPoreParameters );
CEREAL_REGISTER_TYPE( core::conformation::membrane::AqueousPoreParameters )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_membrane_AqueousPoreParameters )
#endif // SERIALIZATION
