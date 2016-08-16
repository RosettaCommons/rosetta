// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/CommonFragmentComparators.fwd.hh
/// @brief provides a few ways to compare fragment candidates (forward definitions)
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_CommonFragmentComparators_fwd_hh
#define INCLUDED_protocols_frag_picker_CommonFragmentComparators_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace frag_picker {

class CompareQueryPosition;
class CompareTotalScore;
class CompareScoreComponent;

typedef utility::pointer::shared_ptr<CompareQueryPosition> CompareQueryPositionOP;
typedef utility::pointer::shared_ptr<CompareQueryPosition const> CompareQueryPositionCOP;

typedef utility::pointer::shared_ptr<CompareTotalScore> CompareTotalScoreOP;
typedef utility::pointer::shared_ptr<CompareTotalScore const> CompareTotalScoreCOP;

typedef utility::pointer::shared_ptr<CompareScoreComponent> CompareScoreComponentOP;
typedef utility::pointer::shared_ptr<CompareScoreComponent const> CompareScoreComponentCOP;

} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_CommonFragmentComparators_FWD_HH */
