// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/CompositeFragmentSelector.cc
/// @brief provides a selector that picks best fragments based on their total score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/CompositeFragmentSelector.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace frag_picker {

/// @brief  Selects desired number of fragments from a given set of candidates
void CompositeFragmentSelector::select_fragments(
	ScoredCandidatesVector1 const& in,
	ScoredCandidatesVector1& out )
{
	if ( selectors_.size()==0 ) return;
	if ( selectors_.size()==1 ) {
		selectors_[1]->select_fragments(in,out);
		return;
	}
	Size n = selectors_.size();
	ScoredCandidatesVector1 tmp, tmp_tmp;
	selectors_[1]->select_fragments(in,tmp);
	for ( Size i=2; i<=n-1; i++ ) {
		selectors_[i]->select_fragments(tmp,tmp_tmp);
		tmp.clear();
		for ( Size j=1; j<=tmp_tmp.size(); j++ ) {
			tmp.push_back(tmp_tmp[j]);
		}
	}
	selectors_[n]->select_fragments(tmp,out);
}

} // frag_picker
} // protocols
