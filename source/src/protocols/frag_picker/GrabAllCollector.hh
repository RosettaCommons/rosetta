// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/GrabAllCollector.hh
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_GrabAllCollector_hh
#define INCLUDED_protocols_frag_picker_GrabAllCollector_hh

// package headers
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/CandidatesCollector.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.fwd.hh>

// utility headers
#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {

// Always always declare your OP typedefs if you're defining a polymorphic class.
// I shouldn't have to do this for you.
class GrabAllCollector;
typedef utility::pointer::shared_ptr< GrabAllCollector > GrabAllCollectorOP;
typedef utility::pointer::shared_ptr< GrabAllCollector const > GrabAllCollectorCOP;


/// @brief Keeps all fragments candidates for the final selection
/// @details The purpose of a collector is to keep fragment candidates to the end
/// of vall processing. This simple implementation keeps all candidates which might
/// use huge memory
class GrabAllCollector: public CandidatesCollector {
public:

	/// @brief create a collector for a given size of a query sequence
	GrabAllCollector(Size query_size)
	{
		storage_.resize(query_size);
		for ( Size i = 1; i <= query_size; i++ ) {
			utility::vector1<std::pair<FragmentCandidateOP,
				scores::FragmentScoreMapOP> > vec;
			storage_[i] = vec;
		}
	}

	/// @brief  Insert a fragment candidate to the container
	inline bool add( ScoredCandidate new_canditate ) override
	{
		storage_[new_canditate.first->get_first_index_in_query()].push_back(
			new_canditate);
		return true;
	}

	inline void clear() override
	{
		for ( Size i=1; i<storage_.size(); ++i ) {
			storage_[i].clear();
		}
	}

	/// @brief prints how many fragments have been collected for each position
	void print_report(
		std::ostream & output,
		scores::FragmentScoreManagerOP
	) const override {
		output << "\n pos  | count |  pos  | count | pos  | count |\n";
		for ( Size i = 1; i <= storage_.size(); ++i ) {
			output << I(5, i) << " |" << I(6, storage_[i].size()) << " |";
			if ( i % 3 == 0 ) {
				output << '\n';
			}
		}
	}

	/// @brief  Check how many candidates have been already collected for a given position
	inline Size count_candidates(Size seq_pos) const override {
		return storage_[seq_pos].size();
	}

	/// @brief  Check how many candidates have been already collected for all positions
	inline Size count_candidates() const override {

		Size response = 0;
		for ( Size i=1; i<=storage_.size(); ++i ) {
			response += storage_[i].size();
		}
		return response;
	}

	/// @brief  Check the size of query sequence that this object knows.
	/// This is mainly to be able to check if it is the same as in the other parts of
	/// fragment picking machinery.
	inline Size query_length() const override {
		return storage_.size();
	}

	/// @brief Inserts candidates from another Collector for a give position in the query
	inline void insert(Size pos, CandidatesCollectorOP collector) override {
		GrabAllCollectorOP c = utility::pointer::dynamic_pointer_cast< protocols::frag_picker::GrabAllCollector > ( collector );
		if ( c == nullptr ) {
			utility_exit_with_message("Cant' cast candidates' collector to GrabAllCollector.");
		}
		for ( Size j=1; j<=storage_[pos].size(); ++j ) {
			ScoredCandidatesVector1 & content = c->get_candidates(pos);
			for ( Size l=1; l<=content.size(); l++ ) storage_[pos].push_back( content[l] );
		}
	}

	/// @brief returns all stored fragment candidates that begins at a given position in a query
	inline ScoredCandidatesVector1 & get_candidates(Size position_in_query) override {
		return storage_.at(position_in_query);
	}

private:
	utility::vector1< ScoredCandidatesVector1 > storage_;
};

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_GrabAllCollector_HH */
