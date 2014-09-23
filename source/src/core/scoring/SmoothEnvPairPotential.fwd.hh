// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SmoothEnvPairPotential.cc
/// @brief  Smooth, differentiable version of centroid env and pair terms
/// @author Frank DiMaio

#ifndef INCLUDED_core_scoring_SmoothEnvPairPotential_fwd_hh
#define INCLUDED_core_scoring_SmoothEnvPairPotential_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <core/types.hh>

//Auto Headers
namespace core {
namespace scoring {

template <class T>
class SigmoidWeightedCenList;

// templated typedefs are not allowed
// do this instead
typedef utility::pointer::shared_ptr< SigmoidWeightedCenList<Real> > SigmoidWeightedCenListRealOP;
typedef utility::pointer::shared_ptr< SigmoidWeightedCenList<numeric::xyzVector<Real> > > SigmoidWeightedCenListVectorOP;

class SmoothEnvPairPotential;

typedef utility::pointer::shared_ptr< SmoothEnvPairPotential > SmoothEnvPairPotentialOP;
typedef utility::pointer::shared_ptr< SmoothEnvPairPotential const > SmoothEnvPairPotentialCOP;

} // ns scoring
} // ns core

#endif
