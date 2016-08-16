// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/AdaptiveScoreHistogram.fwd.hh
/// @brief  AdaptiveScoreHistogram forward declaration
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_AdaptiveScoreHistogram_fwd_hh
#define INCLUDED_protocols_frag_picker_scores_AdaptiveScoreHistogram_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace frag_picker {
namespace scores {

/// @brief forward declaration for score histogram
class AdaptiveScoreHistogram;

typedef utility::pointer::shared_ptr<AdaptiveScoreHistogram>
	AdaptiveScoreHistogramOP;
typedef utility::pointer::shared_ptr<AdaptiveScoreHistogram const>
	AdaptiveScoreHistogramCOP;

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_AtomBasedConstraintsScore_FWD_HH */
