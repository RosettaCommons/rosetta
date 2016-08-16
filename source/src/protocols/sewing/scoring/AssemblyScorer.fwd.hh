// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AssemblyScorer.fwd.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_sewing_scoring_AssemblyScorer_fwd_hh
#define INCLUDED_protocols_sewing_scoring_AssemblyScorer_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace sewing  {
namespace scoring {

class AssemblyScorer;
typedef utility::pointer::shared_ptr< AssemblyScorer > AssemblyScorerOP;
typedef utility::pointer::shared_ptr< AssemblyScorer const > AssemblyScorerCOP;

class AssemblyScoreFunction;
typedef utility::pointer::shared_ptr< AssemblyScoreFunction > AssemblyScoreFunctionOP;
typedef utility::pointer::shared_ptr< AssemblyScoreFunction const > AssemblyScoreFunctionCOP;

} //scoring namespace
} //sewing namespace
} //protocols namespace

#endif


