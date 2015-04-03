// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/quota/ABEGO_SS_Pool.hh
/// @brief class for a quota pool based on secondary structure prediction and ABEGO torsion bin classification
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_quota_ABEGO_SS_Pool_hh
#define INCLUDED_protocols_frag_picker_quota_ABEGO_SS_Pool_hh

#include <protocols/frag_picker/quota/ABEGO_SS_Pool.fwd.hh>
#include <protocols/frag_picker/quota/QuotaPool.hh>
#include <protocols/frag_picker/quota/ABEGO_SS_Map.hh>

// package headers
#include <protocols/frag_picker/LazySortedVector1.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/CommonFragmentComparators.hh>

// utility headers
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>

// C++ headers
#include <string>

//Auto Headers
namespace protocols {
namespace frag_picker {
namespace quota {

typedef LazySortedVector1<std::pair<FragmentCandidateOP,
                scores::FragmentScoreMapOP>, CompareByScoreCombination> BoundedQuotaContainer;
typedef utility::pointer::shared_ptr<BoundedQuotaContainer>
		BoundedQuotaContainerOP;

/// @brief represents a single pool used by quota selector
class ABEGO_SS_Pool : public QuotaPool {
public:
	/// @brief Creates a pool of a given size and name
	/// @param size - total number of fragments from all the pools at a given position
	/// @param name - name assigned to this pool. This in general may be any string that
	///	later allows one control pool's behavior from a flag file
	/// @param abego_ss_pairs - what types of secondary structure and ABEGO this pool is accepting
	/// @param score_components_id - which scores will be used to sort this pool
	/// @param weights - weights for the scores that in general may be different than these used for fragment picking
	/// @param fraction - fraction of this pool in the entire population
	ABEGO_SS_Pool(Size,std::string,utility::vector1<std::pair<Size,Size> >,
	    utility::vector1<Size>,utility::vector1<Real>,Real,Size,Size);

	bool accepts(char ss_type,char abego_type) { return ss_abego_types_->check_status(ss_type,abego_type); }

	virtual ~ABEGO_SS_Pool();

	/// @brief Says how many fragments (in total) may fit into this pool
	virtual Size total_size() const { return storage_->size(); }

	/// @brief Says how many fragments are currently in this pool
	virtual Size current_size() const { return storage_->count_inserted();}

	/// @brief Says how many fragments can still be inserted into this pool
	virtual Size size_left() const { return total_size() - current_size(); }

	virtual bool could_be_accepted(ScoredCandidate candidate) const;

	void resize(Size new_size) {
	    storage_->resize(new_size,new_size*buffer_factor_);
	}

	virtual void set_fraction(Real new_fraction) {

		QuotaPool::set_fraction(new_fraction);
		this_size_ = (Size)(total_size_*new_fraction);
		if(this_size_<20) {
			this_size_ = 20;
		}
		storage_->resize(this_size_,this_size_*buffer_factor_);
	}

	/// @brief Push a fragment candidate into the container
	virtual void push(ScoredCandidate candidate) {
		storage_->push(candidate);
	}

	// Stuff inherited from CandidatesCollector base

	/// @brief  Insert a fragment candidate to the container
	virtual bool add(ScoredCandidate);

	/// @brief removes all candidates from the container
	void clear() { storage_->clear(); }

	/// @brief  Check how many candidates have been already collected for a given position
	/// @details This is a very special case - collector will be used only for a given position.
	/// Thus it returns the total number of inserted candidates, as count_candidates() does
	Size count_candidates(Size) const { return current_size(); }

	/// @brief  Check how many candidates have been already collected for all positions
	Size count_candidates() const { return current_size(); }

	/// @brief  Check the size of query sequence that this object knows.
	/// @details This is a very special case - collector will be used only for a given position and it does NOT
	/// know the tolal size. Thus it returns always 0
	Size query_length() const { return 0; }

	/// @brief Inserts candidates from another Collector for a give position in the query
	/// Candidates may or may not get inserted depending on the candidate
	void insert(Size, CandidatesCollectorOP collector) {
		ABEGO_SS_PoolOP c = utility::pointer::dynamic_pointer_cast< protocols::frag_picker::quota::ABEGO_SS_Pool > ( collector );
		if (c == 0)
			utility_exit_with_message("Cant' cast candidates' collector to ABEGO_SS_Pool.");
		ScoredCandidatesVector1 & content = c->get_candidates(0);
		for(Size l=1;l<=content.size();l++) storage_->push( content[l] );
	}

	/// @brief  Returns all the candidate in this pool
	ScoredCandidatesVector1 & get_candidates( Size /*position_in_query*/ ) {
		return storage_->expose_data();
	}

	/// @brief Describes what has been collected
	void print_report(std::ostream &, scores::FragmentScoreManagerOP) const;

	virtual Real quota_score(ScoredCandidate candidate) const {

	    Real t2(0);
	    for(Size i=1;i<=components_.size();i++)
		t2 += candidate.second->at( components_[i] ) * weights_[i];

	    return t2;
	}

private:
	static char all_abego_[];
	static char all_ss_[];
	Size total_size_;
	Size this_size_;
	Size buffer_factor_;
	ABEGO_SS_MapOP ss_abego_types_;
	utility::vector1<Size> components_;
	utility::vector1<Real> weights_;
	BoundedQuotaContainerOP storage_;
};

} // quota
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_quota_ABEGO_SS_Pool_HH */
