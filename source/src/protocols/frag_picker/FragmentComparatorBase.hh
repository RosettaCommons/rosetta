// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/FragmentComparatorBase.hh
/// @brief defines a base class for all fragment comparators
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_FragmentComparatorBase_hh
#define INCLUDED_protocols_frag_picker_FragmentComparatorBase_hh

// package headers
#include <protocols/frag_picker/FragmentComparatorBase.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <core/types.hh>
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>


namespace protocols {
namespace frag_picker {

class FragmentComparatorBase : public utility::pointer::ReferenceCount {
public:
	virtual bool operator() (
		ScoredCandidate first_candidate,
		ScoredCandidate second_candidate
	) = 0;
};

} // frag_picker
} // protocols

#endif
