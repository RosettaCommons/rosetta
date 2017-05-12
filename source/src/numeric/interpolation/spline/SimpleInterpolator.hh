// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/numeric/interpolation/spline/SimpleInterpolator.hh
/// @brief  Interpolation with cubic splines
/// @author Will Sheffler


#ifndef INCLUDED_numeric_interpolation_spline_SimpleInterpolator_hh
#define INCLUDED_numeric_interpolation_spline_SimpleInterpolator_hh

#include <numeric/interpolation/spline/Interpolator.hh>
#include <utility/vector1.hh>
#include <numeric/interpolation/spline/spline_functions.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace numeric {
namespace interpolation {
namespace spline {

class SimpleInterpolator;
typedef utility::pointer::shared_ptr< SimpleInterpolator > SimpleInterpolatorOP;
typedef utility::pointer::shared_ptr< SimpleInterpolator const > SimpleInterpolatorCOP;

class SimpleInterpolator : public Interpolator {

public:

	SimpleInterpolator(
		utility::vector1<Real> const & x,
		utility::vector1<Real> const & y,
		Real lbdy,
		Real ubdy
	);

	SimpleInterpolator();

	InterpolatorOP clone() const override;

	void interpolate( Real x, Real & y, Real & dy ) const override;

	/// @brief serialize the Interpolator to a json_spirit object
	utility::json_spirit::Value serialize() const override;
	/// @brief deserialize a json_spirit object to a Interpolator
	void deserialize(utility::json_spirit::mObject data) override;

	bool operator == ( Interpolator const & other ) const override;
	bool same_type_as_me( Interpolator const & other ) const override;

	/// @brief accessor for the x-values defining the interpolation knots
	utility::vector1< Real > const & x()   const { return x_; }

	/// @brief accessor for the y-values specified for each knot
	utility::vector1< Real > const & y()   const { return y_; }

	/// @brief accessors for the second derivative values specified for each knot
	utility::vector1< Real > const & ddy() const { return ddy_; }

private:

	utility::vector1<Real> x_, y_, ddy_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

typedef utility::pointer::shared_ptr< Interpolator > InterpolatorOP;

} // end namespace spline
} // end namespace interpolation
} // end namespace numeric

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( numeric_interpolation_spline_SimpleInterpolator )
#endif // SERIALIZATION


#endif
