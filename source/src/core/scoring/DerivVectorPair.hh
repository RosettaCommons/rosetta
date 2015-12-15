// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/DerivVectorPair.hh
/// @brief  Class for storing a pair of derivative vectors, f1 and f2, used in our
///         internal-geometry minimization algorithm.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_DerivVectorPair_hh
#define INCLUDED_core_scoring_DerivVectorPair_hh

// Unit headers
#include <core/scoring/DerivVectorPair.fwd.hh>

// Package headers
#include <core/types.hh>

// Numeric headers
#include <numeric/xyzVector.hh>


namespace core {
namespace scoring {

/// @brief A glorified struct for holding f1/f2 vectors used to represent the
/// derivatives in our internal-geometry based minimization algorithms.
class DerivVectorPair
{
public:
	DerivVectorPair() : f1_( 0.0 ), f2_( 0.0 ) {}
	Vector & f1() { return f1_; }
	Vector & f2() { return f2_; }
	Vector const & f1() const { return f1_; }
	Vector const & f2() const { return f2_; }

private:
	Vector f1_;
	Vector f2_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // scoring
} // core

#endif
