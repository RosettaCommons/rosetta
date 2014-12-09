// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/optimization/AtomTreeMinimizer.hh
/// @brief  High-level atom tree minimizer class for symmetrical minimization
/// @author Ingemar Andre

#ifndef INCLUDED_core_optimization_symmetry_SymAtomTreeMinimizer_fwd_hh
#define INCLUDED_core_optimization_symmetry_SymAtomTreeMinimizer_fwd_hh


namespace core {
namespace optimization {
namespace symmetry {


// Forward
class SymAtomTreeMinimizer;
typedef utility::pointer::owning_ptr< AtomTreeMinimizer       > SymAtomTreeMinimizerOP;
typedef utility::pointer::owning_ptr< AtomTreeMinimizer const > SymAtomTreeMinimizerCOP;


} // symmetry
} // namespace optimization
} // namespace core


#endif // INCLUDED_core_optimization_symmetry_SymAtomTreeMinimizer_FWD_HH
