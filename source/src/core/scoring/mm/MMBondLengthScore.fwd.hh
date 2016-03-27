// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/mm/MMBondLengthScore.fwd.hh
/// @brief  Molecular mechanics bond length score class
/// @author Frank DiMaio (based on Colin Smith's MMBondAngle potential)


#ifndef INCLUDED_core_scoring_mm_MMBondLengthScore_fwd_hh
#define INCLUDED_core_scoring_mm_MMBondLengthScore_fwd_hh

#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

namespace core {
namespace scoring {
namespace mm {


// Forward
class MMBondLengthScore;

typedef  utility::pointer::weak_ptr< MMBondLengthScore > MMBondLengthScoreAP;
typedef  utility::pointer::weak_ptr< MMBondLengthScore const > MMBondLengthScoreCAP;
typedef  utility::pointer::shared_ptr< MMBondLengthScore > MMBondLengthScoreOP;
typedef  utility::pointer::shared_ptr< MMBondLengthScore const > MMBondLengthScoreCOP;

} // namespace mm
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_mm_mm_BondLength_score_FWD_HH
