// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/matdes/MatDesPointMutationCalculator.fwd.hh
/// @brief	this is a modified version of chris king's PointMutationCalculator with additional functionality that is currently not compatible with all of the ParetoOpt functionality. please note that this has been checked into master in its current state in response to requests from others to use this modified version of chris king's GreedyOptMutationMover. although this is still a somewhat developmental piece of code, it has currently been left in src/protocols/matdes/ to avoid issues with intra-library level dependencies.  
/// @author jacob bale (balej@uw.edu)

#ifndef INCLUDED_protocols_matdes_MatDesPointMutationCalculator_fwd_hh
#define INCLUDED_protocols_matdes_MatDesPointMutationCalculator_fwd_hh


// Utility headers
// AUTO-REMOVED #include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace matdes {

class MatDesPointMutationCalculator;
typedef utility::pointer::owning_ptr< MatDesPointMutationCalculator >  MatDesPointMutationCalculatorOP;
typedef utility::pointer::owning_ptr< MatDesPointMutationCalculator const >  MatDesPointMutationCalculatorCOP;


} // namespace matdes
} // namespace protocols

#endif
