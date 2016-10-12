// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/quota/QuotaCollector.hh
/// @brief  Pure virtual base class for a container holding fragment candidates
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_quota_QuotaCollector_hh
#define INCLUDED_protocols_frag_picker_quota_QuotaCollector_hh

// type headers
#include <utility/pointer/ReferenceCount.hh>

// package headers
#include <protocols/frag_picker/quota/QuotaCollector.fwd.hh>
#include <protocols/frag_picker/quota/QuotaPool.hh>
#include <protocols/frag_picker/CandidatesCollector.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.fwd.hh>

#include <core/fragment/SecondaryStructure.fwd.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace frag_picker {
namespace quota {

/// @brief A base class for collecting fragments.
/// @details The purpose of a collector is to keep fragment candidates to the end
/// of vall processing. Then a selector will go through all the candidates stored
/// in a collector and select the final fragments
/// @see GrabAll collector for a possible implementation
class QuotaCollector : public CandidatesCollector {
public:

	QuotaCollector(core::Size query_size, core::Size frag_size)  { storage_.resize(query_size-frag_size+1); frag_size_ = frag_size; }

	/// @brief  Insert a fragment candidate to the container
	bool add(std::pair<FragmentCandidateOP, scores::FragmentScoreMapOP>);

	/// @brief removes all candidates from the container
	void clear();

	/// @brief  Check how many candidates have been already collected for a given position
	core::Size count_candidates(core::Size) const;

	/// @brief  Check how many candidates have been already collected for all positions
	core::Size count_candidates() const;

	/// @brief  Check the size of query sequence that this object knows.
	/// This is mainly to be ale to check if it is the same as in the other parts of
	/// fragment picking machinery.
	core::Size query_length() const { return storage_.size(); }

	//Undefined, commenting out to fix PyRosetta build  ScoredCandidatesVector1 const & get_candidates(core::Size position_in_query) const;
	ScoredCandidatesVector1 & get_candidates(core::Size position_in_query);

	/// @brief Describes what has been collected
	void print_report(std::ostream & output,
		scores::FragmentScoreManagerOP scoring
	) const;

	/// @brief Inserts candidates from another QuotaCollector for a give position in the query
	/// Candidates may or may not get inserted depending on the candidate
	void insert(core::Size pos, CandidatesCollectorOP collector) {
		QuotaCollectorOP c = utility::pointer::dynamic_pointer_cast< protocols::frag_picker::quota::QuotaCollector > ( collector );
		if ( c == 0 ) {
			utility_exit_with_message("Cant' cast candidates' collector to QuotaCollector. Is quota set up correctly?");
		}
		for ( core::Size j=1; j<=storage_[pos].size(); ++j ) {
			ScoredCandidatesVector1 & content = c->get_pool(pos, j)->get_candidates(0);
			for ( core::Size l=1; l<=content.size(); l++ ) storage_[pos][j]->push( content[l] );
		}
	}

	/// @brief list all registered pools with their capacity
	void list_pools(std::ostream & where) const;

	/// @brief prints the number of quota pools for a given positio in query
	core::Size count_pools(core::Size position) const {
		return storage_[position].size();
	}

	QuotaPoolCOP get_pool( core::Size position, core::Size pool_id ) const {
		return storage_[position][pool_id];
	}

	QuotaPoolOP get_pool( core::Size position, core::Size pool_id ) {
		return storage_[position][pool_id];
	}

	void add_pool(core::Size position,QuotaPoolOP the_pool) {
		if ( position <= storage_.size() ) {
			storage_[position].push_back(the_pool);
		}
	}

	void renormalize_quota_pools();

	void attach_secondary_structure_pools(core::Real,core::fragment::SecondaryStructureOP,
		std::string, core::Size, utility::vector1<core::Size>,utility::vector1<core::Real>,core::Size);
private:
	core::Size frag_size_;
	ScoredCandidatesVector1 frags_for_pos_;
	utility::vector1< utility::vector1<QuotaPoolOP> > storage_;
};

} // quota
} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_quota_QuotaCollector_HH */
