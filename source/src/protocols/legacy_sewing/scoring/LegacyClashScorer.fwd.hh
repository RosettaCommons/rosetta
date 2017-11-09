// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyClashScorer.fwd.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_legacy_sewing_scoring_LegacyClashScorer_fwd_hh
#define INCLUDED_protocols_legacy_sewing_scoring_LegacyClashScorer_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace legacy_sewing  {
namespace scoring {

class LegacyClashScorer;
typedef utility::pointer::shared_ptr< LegacyClashScorer > LegacyClashScorerOP;
typedef utility::pointer::shared_ptr< LegacyClashScorer const > LegacyClashScorerCOP;

} //scoring namespace
} //legacy_sewing namespace
} //protocols namespace

#endif


