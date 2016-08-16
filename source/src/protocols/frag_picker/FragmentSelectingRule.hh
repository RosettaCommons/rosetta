// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/FragmentSelectingRule.hh
/// @brief provides a selector that picks best fragments based on their total score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_FragmentSelectingRule_hh
#define INCLUDED_protocols_frag_picker_FragmentSelectingRule_hh

#include <protocols/frag_picker/FragmentSelectingRule.fwd.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>


// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace frag_picker {

/// @brief selects a given number of fragments just by selecting the best ones
class FragmentSelectingRule: public utility::pointer::ReferenceCount {
public:

	/// @brief  Constructor sets the desired number of fragments.
	FragmentSelectingRule(Size frags_per_pos) {
		frags_per_pos_ = frags_per_pos;
	}

	/// @brief  Selects desired number of fragments from a given candidates
	virtual void select_fragments( ScoredCandidatesVector1 const&, ScoredCandidatesVector1& ) = 0;

	/// @brief Says how many fragments will be selected for each position in a query sequence
	inline Size frags_per_pos() {

		return frags_per_pos_;
	}

	virtual ~FragmentSelectingRule() {
	}

private:
	Size frags_per_pos_;
};

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_FragmentSelectingRule_HH */
