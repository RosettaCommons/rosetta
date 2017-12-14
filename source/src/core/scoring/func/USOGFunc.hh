// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/func/USOGFunc.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef CORE_SCORING_CONSTRAINTS_USOGFUNC_HH_
#define CORE_SCORING_CONSTRAINTS_USOGFUNC_HH_

// Package headers
#include <core/scoring/func/USOGFunc.fwd.hh>
#include <core/scoring/func/Func.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <iostream>
#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

/// @brief Unnormalized, unbounded sum of Gaussians constraint
class USOGFunc : public Func {
public:
	static Real background_prob;

	/// @brief Used in conjunction with read_data() to initialize a new instance
	USOGFunc() {};

	/// @brief Constructs a new instance with a single gaussian
	USOGFunc(core::Real mean, core::Real std_dev, core::Real weight=1);

	/// @brief Constructs a new instance from the specified lists of means,
	/// standard deviations, and weights. Assumes that all lists have equal
	/// length and weights sum to 1.
	USOGFunc( utility::vector1<core::Real> const & means,
		utility::vector1<core::Real> const & std_devs,
		utility::vector1<core::Real> const & weights);

	/// @brief No-op virtual destructor
	~USOGFunc() {}

	FuncOP clone() const;
	virtual bool operator == ( Func const & other ) const;
	virtual bool same_type_as_me( Func const & other ) const;

	/// @brief Returns a value representing this function evaluated at a given point
	core::Real func(const core::Real x) const;

	/// @brief Returns a value representing the derivative of this function evaluated at a given point
	core::Real dfunc(const core::Real x) const;

	/// @brief Initializes this function from the given input stream
	void read_data(std::istream& in);

	/// @brief Writes the definition of this function to the specific output stream
	void show_definition(std::ostream& out) const;

	/// @brief Returns the number of Gaussian components
	core::Size numGaussians() const;

private:
	/// @brief Resets all information associated with this instance
	void resetInstance();

	utility::vector1<core::Real> means_;
	utility::vector1<core::Real> std_devs_;
	utility::vector1<core::Real> weights_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


/// @brief Reads and returns a single floating point value from the specified
/// input stream, aborting if either the failbit or badbit is flipped.
core::Real readValueOrDie(std::istream& in);

}  // namespace constraints
}  // namespace scoring
}  // namespace core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_func_USOGFunc )
#endif // SERIALIZATION


#endif  // CORE_SCORING_CONSTRAINTS_USOGFUNC_HH_
