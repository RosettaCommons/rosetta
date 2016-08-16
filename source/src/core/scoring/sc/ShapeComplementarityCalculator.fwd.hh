// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/sc/ShapeComplementarityCalculator.fwd.hh
/// @brief  Headers for the Shape Complementarity Calculator
/// @author Luki Goldschmidt (luki@mbi.ucla.edu)


#ifndef INCLUDED_core_scoring_sc_ShapeComplementarityCalculator_fwd_hh
#define INCLUDED_core_scoring_sc_ShapeComplementarityCalculator_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace sc {

// Forward
class ShapeComplementarityCalculator;

// Types
typedef utility::pointer::shared_ptr< ShapeComplementarityCalculator > ShapeComplementarityCalculatorOP;
typedef utility::pointer::shared_ptr< ShapeComplementarityCalculator const > ShapeComplementarityCalculatorCOP;


} // namespace sc
} // namespace scoring
} // namespace core

#endif
