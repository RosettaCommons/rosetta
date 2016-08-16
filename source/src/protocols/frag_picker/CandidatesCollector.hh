// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/CandidatesCollector.hh
/// @brief  Pure virtual base class for a container holding fragment candidates
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_CandidatesCollector_hh
#define INCLUDED_protocols_frag_picker_CandidatesCollector_hh

// type headers
#include <utility/pointer/ReferenceCount.hh>

// package headers
#include <protocols/frag_picker/CandidatesCollector.fwd.hh>
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <protocols/frag_picker/scores/FragmentScoreManager.fwd.hh>

#ifdef PYROSETTA
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>
#endif

namespace protocols {
namespace frag_picker {

/// TODO write collector for mixed fragment lengths
/// @brief A base class for collecting fragments.
/// @details The purpose of a collector is to keep fragment candidates to the end
/// of vall processing. Then a selector will go through all the candidates stored
/// in a collector and select the final fragments
/// @see GrabAll collector for a possible implementation
class CandidatesCollector: public utility::pointer::ReferenceCount {
public:
	/// @brief  Insert a fragment candidate to the container
	virtual bool add( ScoredCandidate ) = 0;

	/// @brief removes all candidates from the container
	virtual void clear() = 0;

	/// @brief inserts candidates from another collector
	/// Candidates may or may not get inserted depending on the candidate and type of storage
	virtual void insert( Size, CandidatesCollectorOP ) = 0;

	/// @brief  Check how many candidates have been already collected for a given position
	virtual Size count_candidates(Size) const = 0;

	/// @brief  Check how many candidates have been already collected for all positions
	virtual Size count_candidates() const = 0;

	/// @brief  Check the size of query sequence that this object knows.
	/// This is mainly to be ale to check if it is the same as in the other parts of
	/// fragment picking machinery.
	virtual Size query_length() const =0;

	virtual ScoredCandidatesVector1 & get_candidates( Size position_in_query) = 0;

	/// @brief Describes what has been collected
	virtual void print_report(
		std::ostream & output,
		scores::FragmentScoreManagerOP scoring
	) const = 0;
};

} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_CandidatesCollector_HH */
