// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/GaussianChainTripleFunc.fwd.hh
/// @brief forward declaration for Mixture functions
/// @author Rhiju Das


#ifndef INCLUDED_core_scoring_constraints_GaussianChainTripleFunc_fwd_hh
#define INCLUDED_core_scoring_constraints_GaussianChainTripleFunc_fwd_hh

#include <utility/pointer/owning_ptr.fwd.hh>

// C++ Headers

namespace core {
namespace scoring {
namespace constraints {

class GaussianChainTripleFunc;
typedef utility::pointer::owning_ptr< GaussianChainTripleFunc > GaussianChainTripleFuncOP;
typedef utility::pointer::owning_ptr< GaussianChainTripleFunc const > GaussianChainTripleFuncCOP;

}
}
}

#endif
