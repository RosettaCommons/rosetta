// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/DiversifyDihedralsSelector.hh
/// @brief provides a selector that picks best fragments based on their total score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_DiversifyDihedralsSelector_hh
#define INCLUDED_protocols_frag_picker_DiversifyDihedralsSelector_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/FragmentSelectingRule.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <core/types.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace frag_picker {

/// @brief selects fragments by running several selectors
class DiversifyDihedralsSelector: public FragmentSelectingRule {
public:

	/// @brief  Constructor sets the desired number of fragments and crmsd cutoff.
	DiversifyDihedralsSelector(core::Size frags_per_pos, core::Real cutoff) : FragmentSelectingRule(frags_per_pos) {
		cutoff_ = cutoff;
	}

	/// @brief  Selects desired number of fragments from a given set of candidates
	virtual void select_fragments( ScoredCandidatesVector1 const&, ScoredCandidatesVector1& );

	virtual ~DiversifyDihedralsSelector() {
	}

	core::Real dihedral_rmsd(FragmentCandidateOP,FragmentCandidateOP);

private:
	core::Real cutoff_;
	utility::vector1<core::Real> phi_;
	utility::vector1<core::Real> psi_;
};

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_DiversifyDihedralsSelector_HH */
