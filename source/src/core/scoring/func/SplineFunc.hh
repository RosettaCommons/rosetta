// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (C) 199x-2009 University of Washington
// (C) 199x-2009 University of California Santa Cruz
// (C) 199x-2009 University of California San Francisco
// (C) 199x-2009 Johns Hopkins University
// (C) 199x-2009 University of North Carolina, Chapel Hill
// (C) 199x-2009 Vanderbilt University

/// @file   core/scoring/func/SplineFunc.hh
/// @brief  Constraint function for looking up data in a histogram over which a spline is created
/// @author Stephanie Hirst (stephanie.j.hirst@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_func_SplineFunc_hh
#define INCLUDED_core_scoring_func_SplineFunc_hh

// Unit Headers
#include <core/scoring/func/SplineFunc.fwd.hh>

// Package Headers
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/Func.hh>

// Project Headers

// Utility and Numeric Headers
#include <utility/vector1.hh>
#include <numeric/interpolation/spline/Interpolator.hh>

// C++ Headers
#include <iosfwd>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

// TAG, cutoff, slope, intercept
typedef utility::vector1<std::tuple<std::string, platform::Real, platform::Real, platform::Real>> boundary_fn_type;

class SplineFunc : public Func {
public:

	// @brief SplineFunc construction and destruction
	SplineFunc();

	// @brief SplineFunc construction from code
	SplineFunc(
		std::string const & KB_description,
		core::Real const weight,
		core::Real const exp_val,
		core::Real const bin_size,
		utility::vector1<core::Real> const & bins_vect,
		utility::vector1<core::Real> const & potential_vect,
		boundary_fn_type const & boundary_functions = boundary_fn_type());

	~SplineFunc() override;

	/// @brief returns a clone of this SplineFunc
	FuncOP clone() const override { return utility::pointer::make_shared< SplineFunc >( *this ); }
	bool operator == ( Func const & other ) const override;
	bool same_type_as_me( Func const & other ) const override;

	/// @brief Compare this SplineFunc to another one, and determine whether
	/// the two are equal within some threshold.
	/// @details The float threshold is used for comparing floating-point values.
	bool is_approximately_equal( SplineFunc const & other, core::Real const float_threshold ) const;

	/// @brief return SplineFunc member variables
	core::Real get_exp_val() const;

	std::string const & get_filename() const;

	std::string const & get_KB_description() const;

	core::Real get_weight() const;

	core::Real get_bin_size() const;

	core::Size get_bins_vect_size() const;

	core::Size get_potential_vect_size() const;

	core::Real get_lower_bound_x() const;

	core::Real get_upper_bound_x() const;

	core::Real get_lower_bound_y() const;

	core::Real get_upper_bound_y() const;

	core::Real get_lower_bound_dy() const;

	core::Real get_upper_bound_dy() const;

	utility::vector1<core::Real> const & get_bins_vect() const;

	utility::vector1<core::Real> const & get_potential_vect() const;

	/// @brief Initialize this SplineFunc from the given izstream.
	/// @details Triggers read from disk UNLESS the first entry (filename) is "NONE" (upper or lower case or any mixture).
	/// In that case, the spline is read from the rest of the line.
	void read_data( std::istream &in ) override;

	/// @brief Returns the value of this SplineFunc evaluated at distance x.
	core::Real func( core::Real const x ) const override;

	/// @brief Returns the value of the first derivative of this SplineFunc at distance x.
	core::Real dfunc( core::Real const x ) const override;

	/// @brief show the definition of this SplineFunc to the specified output stream.
	void show_definition( std::ostream &out ) const override;

	/// @brief show some sort of stringified representation of the violations for this constraint.
	core::Size show_violations( std::ostream &out, core::Real x, core::Size verbose_level, core::Real threshold = 1 ) const override;

private:

	/// @brief Are two values equal within some threshold?
	bool equal_within_thresh( core::Real const val1, core::Real const val2, core::Real const threshold ) const;

private:

	core::Real exp_val_;
	std::string filename_;
	std::string KB_description_;
	core::Real weight_;
	core::Real bin_size_;
	core::Real lower_bound_x_;
	core::Real lower_bound_y_;
	core::Real upper_bound_x_;
	core::Real upper_bound_y_;
	core::Real lower_bound_dy_;
	core::Real upper_bound_dy_;
	utility::vector1<core::Real> bins_vect_;
	utility::vector1<core::Real>::size_type bins_vect_size_;
	utility::vector1<core::Real> potential_vect_;
	utility::vector1<core::Real>::size_type potential_vect_size_;
	numeric::interpolation::spline::InterpolatorOP interpolator_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // SplineFunc class

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_SplineFunc )
#endif // SERIALIZATION


#endif
