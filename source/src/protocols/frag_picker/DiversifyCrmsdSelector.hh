// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/DiversifyCrmsdSelector.hh
/// @brief provides a selector that picks best fragments based on their total score
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_DiversifyCrmsdSelector_hh
#define INCLUDED_protocols_frag_picker_DiversifyCrmsdSelector_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/FragmentSelectingRule.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>


#include <core/types.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <utility/vector1.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {

/// @brief selects fragments by running several selectors
class DiversifyCrmsdSelector: public FragmentSelectingRule {
public:

	/// @brief  Constructor sets the desired number of fragments and crmsd cutoff.
	DiversifyCrmsdSelector(core::Size frags_per_pos, core::Real cutoff) : FragmentSelectingRule(frags_per_pos) {
		cutoff_ = cutoff;
	}

	/// @brief  Selects desired number of fragments from a given set of candidates
	virtual void select_fragments( ScoredCandidatesVector1 const&, ScoredCandidatesVector1& );

	virtual ~DiversifyCrmsdSelector() { }

private:
	void copy_coordinates(FragmentCandidateOP src, ObjexxFCL::FArray2D_double & dst);
	core::Real cutoff_;
	ObjexxFCL::FArray2D_double fi_;
	ObjexxFCL::FArray2D_double fj_;
	ObjexxFCL::FArray1D_double weights_;
};

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_DiversifyCrmsdSelector_HH */
