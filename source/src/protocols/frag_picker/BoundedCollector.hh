// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/BoundedCollector.hh
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_BoundedCollector_hh
#define INCLUDED_protocols_frag_picker_BoundedCollector_hh

// package headers
#include <protocols/frag_picker/LazySortedVector1.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/BoundedCollector.fwd.hh>
#include <protocols/frag_picker/CandidatesCollector.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>

#include <protocols/frag_picker/CommonFragmentComparators.hh>

// utility headers
#include <core/types.hh>


/// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

#include <iostream>

#include <utility/vector1.hh>

//#include <algorithm>

namespace protocols {
namespace frag_picker {


/// @brief Keeps the N best fragments candidates for the final selection
/// @details The purpose of a collector is to keep the best fragment candidates to the end
/// of vall processing. In particular, the capacity of this collector may be larger than
/// the number of fragments one wants to get
template< class StrictWeakOrdering >
class BoundedCollector: public CandidatesCollector {
public:
	typedef utility::pointer::shared_ptr< BoundedCollector< StrictWeakOrdering > > BoundedCollectorOP;
	typedef utility::pointer::shared_ptr< BoundedCollector< StrictWeakOrdering > const> BoundedCollectorCOP;

public:

	/// @brief create a collector for a given size of a query sequence
	BoundedCollector(
		Size query_size,
		Size max_frags_per_pos,
		StrictWeakOrdering fragment_comparator,
		Size n_score_terms,
		Size buffer_factor = 5
	) {

		FragmentCandidateOP worst_f( new FragmentCandidate(1,1,0,1) );
		scores::FragmentScoreMapOP worst_s( new scores::FragmentScoreMap(n_score_terms) );
		for(Size i=1;i<=n_score_terms;i++)
		    worst_s->set_score_component(99999.9999,i);

		for (Size i = 1; i <= query_size; i++) {
			LazySortedVector1< std::pair< FragmentCandidateOP,
					scores::FragmentScoreMapOP >, StrictWeakOrdering > queue(
					fragment_comparator, max_frags_per_pos,max_frags_per_pos*buffer_factor);
			queue.set_worst( std::pair<FragmentCandidateOP,
			                scores::FragmentScoreMapOP>(worst_f,worst_s) );
			storage_.push_back(queue);
		}
	}

	/// @brief  Insert a fragment candidate to the container
	inline bool add( ScoredCandidate new_canditate) {

		return storage_[new_canditate.first->get_first_index_in_query()].push(
				new_canditate);
	}

	/// @brief  Check how many candidates have been already collected for a given position
	/// APL Note: you cannot have inlined virtual functions
	inline Size count_candidates(Size seq_pos) const {
		return storage_[seq_pos].size();
	}

	/// @brief  Check how many candidates have been already collected for all positions
	/// APL Note: you cannot have inlined virtual functions
	inline Size count_candidates() const {

		Size response = 0;
		for(Size i=1;i<=storage_.size();++i)
			response += storage_[i].size();
		return response;
	}

	/// @brief  Check the size of query sequence that this object knows.
	/// This is mainly to be able to check if it is the same as in the other parts of
	/// fragment picking machinery.
	inline Size query_length() const {
		return storage_.size();
	}

	/// @brief Inserts candidates from another Collector for a give position in the query
	/// Candidates may or may not get inserted depending on the candidate
	void insert(Size pos, CandidatesCollectorOP collector) {
		BoundedCollectorOP c = utility::pointer::dynamic_pointer_cast< protocols::frag_picker::BoundedCollector<class protocols::frag_picker::CompareTotalScore> > ( collector );
		if (c == 0) {
			utility_exit_with_message("Cant' cast candidates' collector to BoundedCollector.");
		}
		ScoredCandidatesVector1 & content = c->get_candidates(pos);
		for(Size l=1;l<=content.size();l++) storage_.at(pos).push( content[l] );
	}

	/// @brief returns all stored fragment candidates that begins at a given position in a query
	inline 	ScoredCandidatesVector1 & get_candidates( Size position_in_query ) {
		return storage_.at(position_in_query).expose_data();
	}

	inline void clear() {

	    for (Size i_pos = 1; i_pos <= storage_.size(); ++i_pos)
		storage_[i_pos].clear();
	}

	/// @brief prints how many candidates have been collected for each position and how good they are
	void print_report(std::ostream& out, scores::FragmentScoreManagerOP scoring) const {
		using namespace ObjexxFCL::format;
		out
				<< "\n pos  count   best     worst  | pos  count   best    worst   | pos  count    best    worst  |\n";
		Size cnt = 0;
		for (Size i_pos = 1; i_pos <= storage_.size(); ++i_pos) {

			if (storage_[i_pos].size() <= 1) {
				out << I(4, i_pos) << "      0                   |";
			} else {
				out << I(4, i_pos) << " " << I(6, storage_[i_pos].size())
						<< " ";
				out << F(8, 2, scoring->total_score(
						storage_[i_pos].peek_front().second)) << " ";
				out << F(8, 2, scoring->total_score(
						storage_[i_pos].peek_back().second)) << " |";
			}
			++cnt;
			if (cnt % 3 == 0)
				out << '\n';
		}
		out << std::endl;
	}

private:
	utility::vector1< LazySortedVector1< std::pair< FragmentCandidateOP, scores::FragmentScoreMapOP> , StrictWeakOrdering> > storage_;
};

// Concrete version for PyRosetta
class BoundedCollector_CompareTotalScore : public BoundedCollector<CompareTotalScore>
{
public:
	BoundedCollector_CompareTotalScore(Size query_size, Size max_frags_per_pos, CompareTotalScore fragment_comparator,Size n_score_terms,Size buffer_factor = 5)
	 : BoundedCollector<CompareTotalScore>(query_size, max_frags_per_pos, fragment_comparator, n_score_terms, buffer_factor) {};
};

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_BoundedCollector_HH */
