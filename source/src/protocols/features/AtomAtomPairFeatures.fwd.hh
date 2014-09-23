// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/AtomAtomPairFeatures.fwd.hh
/// @brief  report atom-atom pair geometry and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_AtomAtomPairFeatures_fwd_hh
#define INCLUDED_protocols_features_AtomAtomPairFeatures_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols{
namespace features{

class AtomAtomPairFeatures;
typedef utility::pointer::shared_ptr< AtomAtomPairFeatures > AtomAtomPairFeaturesOP;
typedef utility::pointer::shared_ptr< AtomAtomPairFeatures const > AtomAtomPairFeaturesCOP;

}// namespace
}// namespace

#endif // include guard
