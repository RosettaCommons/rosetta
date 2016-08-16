// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/mm/MMTorsionScore.fwd.hh
/// @brief  core::scoring::mm::MMTorsionScore forward declarations
/// @author P. Douglas Renfrew (renfrew@nyu.edu)


#ifndef INCLUDED_core_scoring_mm_MMTorsionScore_fwd_hh
#define INCLUDED_core_scoring_mm_MMTorsionScore_fwd_hh

#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

namespace core {
namespace scoring {
namespace mm {


// Forward
class MMTorsionScore;

typedef  utility::pointer::weak_ptr< MMTorsionScore > MMTorsionScoreAP;
typedef  utility::pointer::weak_ptr< MMTorsionScore const > MMTorsionScoreCAP;
typedef  utility::pointer::shared_ptr< MMTorsionScore > MMTorsionScoreOP;
typedef  utility::pointer::shared_ptr< MMTorsionScore const > MMTorsionScoreCOP;

} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_mm_torsion_score_FWD_HH
