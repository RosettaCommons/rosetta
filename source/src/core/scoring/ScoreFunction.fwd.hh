// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.fwd.hh
/// @brief  core::scoring::ScoreFunction forward declarations
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_core_scoring_ScoreFunction_fwd_hh
#define INCLUDED_core_scoring_ScoreFunction_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <boost/shared_ptr.hpp>


namespace core {
namespace scoring {


// Forward
class ScoreFunction;

typedef utility::pointer::shared_ptr< ScoreFunction > ScoreFunctionOP;
typedef utility::pointer::shared_ptr< ScoreFunction const > ScoreFunctionCOP;
typedef boost::shared_ptr < ScoreFunction > ScoreFunctionSP;

} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_ScoreFunction_FWD_HH
