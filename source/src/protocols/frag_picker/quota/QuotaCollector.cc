// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/quota/QuotaCollector.hh
/// @brief  Pure virtual base class for a container holding fragment candidates
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


// package headers
#include <protocols/frag_picker/quota/QuotaCollector.hh>
#include <protocols/frag_picker/quota/QuotaPool.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.fwd.hh>
#include <protocols/frag_picker/quota/SecondaryStructurePool.hh>
#include <core/fragment/SecondaryStructure.hh>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>

#include <iostream>
#include <iomanip>

namespace protocols {
namespace frag_picker {
namespace quota {

static THREAD_LOCAL basic::Tracer trQuotaCollector(
	"protocols.frag_picker.quota.QuotaCollector");

bool QuotaCollector::add(std::pair<FragmentCandidateOP, scores::FragmentScoreMapOP> candidate) {

	Size pos = candidate.first->get_first_index_in_query();
	bool if_inserted = false;
	for ( Size j=1; j<=storage_[pos].size(); ++j ) {
		//     trQuotaCollector.Debug<< "Trying: "<<storage_[pos][j]->get_pool_name()<<" pos: "<<pos<<" pool_id: "<<j<<" of "<<storage_[pos].size()<<std::endl;
		bool tmp = storage_[pos][j]->add(candidate);
		if_inserted = if_inserted || tmp;
	}

	return if_inserted;
}

void QuotaCollector::list_pools(std::ostream & where) const {

	for ( Size i=1; i<storage_.size(); ++i ) {
		where << std::setw(3) << i;
		for ( Size j=1; j<storage_[i].size(); ++j ) {
			where << std::setw(10) << storage_[i][j]->get_pool_name() << "(" << std::setw(3) << storage_[i][j]->total_size() << ") ";
		}
		where << std::endl;
	}
}


void QuotaCollector::clear() {

	for ( Size i=1; i<=storage_.size(); ++i ) {
		for ( Size j=1; j<storage_[i].size(); ++j ) {
			storage_[i][j]->clear();
		}
	}
}

Size QuotaCollector::count_candidates(Size pos) const {

	Size cnt = 0;
	for ( Size j=1; j<storage_[pos].size(); ++j ) {
		cnt += storage_[pos][j]->count_candidates();
	}

	return cnt;
}

Size QuotaCollector::count_candidates() const {
	Size cnt = 0;
	for ( Size i=1; i<=storage_.size(); ++i ) {
		for ( Size j=1; j<storage_[i].size(); ++j ) {
			cnt += storage_[i][j]->count_candidates();
		}
	}

	return cnt;
}

/// @brief Describes what has been collected
void QuotaCollector::print_report(
	std::ostream & output,
	scores::FragmentScoreManagerOP //scoring
) const {

	using namespace ObjexxFCL::format;

	output<<"QuotaCollector contains the following number of fragments at each position:\n";
	for ( Size i=1; i<=storage_.size(); ++i ) {
		output<< I(4, i);
		for ( Size j=1; j<=storage_[i].size(); ++j ) {
			output<<" "<<RJ(10,storage_[i][j]->get_pool_name())
				<<" "<<I(3,storage_[i][j]->current_size());
		}
		output<<"\n";
	}
	output<<std::endl;
}

ScoredCandidatesVector1 & QuotaCollector::get_candidates(Size position_in_query) {

	frags_for_pos_.clear();
	for ( Size j=1; j<=storage_[position_in_query].size(); ++j ) {
		for ( Size k=1; k<=storage_[position_in_query][j]->count_candidates(); k++ ) {
			ScoredCandidatesVector1 & content = storage_[position_in_query][j]->get_candidates(0);
			for ( Size l=1; l<=content.size(); l++ ) {
				frags_for_pos_.push_back( content[l] );
			}
		}
	}
	trQuotaCollector << frags_for_pos_.size()<<" candidates collected for pos. "<<position_in_query<<std::endl;

	return frags_for_pos_;
}


void QuotaCollector::renormalize_quota_pools() {

	for ( Size i=1; i<=storage_.size(); ++i ) {
		Real sum = 0.0;
		for ( Size j=1; j<=i; ++j ) {
			sum += storage_[i][j]->get_fraction();
		}
		for ( Size j=1; j<=i; ++j ) {
			storage_[i][j]->set_fraction( storage_[i][j]->get_fraction()/sum );
		}
	}
}


void QuotaCollector::attach_secondary_structure_pools(
	Real prediction_fraction,
	core::fragment::SecondaryStructureOP prediction,
	std::string name,
	Size n_candidates,
	utility::vector1<Size> components,
	utility::vector1<Real> weights,
	Size n_scores
) {

	Size n;
	Real f;
	Size middle = frag_size_ / 2 + 1;
	for ( Size i=1; i<=storage_.size(); ++i ) {

		trQuotaCollector.Trace << "Quota pools at pos: "<<i<<std::endl;
		f = prediction->helix_fraction(i+middle-1)*prediction_fraction;
		n = (Size) (n_candidates * f );
		Size buffer_factor = 3;
		if ( n == 0 ) {
			trQuotaCollector.Debug << "Pool >"<<name<<"< for H would have size 0, not created"<<std::endl;
		} else {
			SecondaryStructurePoolOP pool( new SecondaryStructurePool(n_candidates,name+":H",'H',components,weights,f,n_scores,buffer_factor) );
			storage_[i].push_back( pool );
			if ( trQuotaCollector.Trace.visible() ) {
				trQuotaCollector.Trace << "Pool >"<<name<<":E< added at query pos: "<<i<<" as number: "<<storage_[i].size()<<std::endl;
				trQuotaCollector.Trace<<"Scoring scheme for the quota pool sorting is:";
				for ( Size l=1; l<=weights.size(); l++ ) {
					trQuotaCollector.Trace<<"\n\t"<<components[l]<<"\t"<<weights[l];
				}
				trQuotaCollector.Trace<<std::endl;
			}
		}

		f = prediction->sheet_fraction(i+middle-1)*prediction_fraction;
		n = (Size) (n_candidates * f );
		if ( n == 0 ) {
			trQuotaCollector.Debug << "Pool >"<<name<<"< for E would have size 0, not created"<<std::endl;
		} else {
			SecondaryStructurePoolOP pool( new SecondaryStructurePool(n_candidates,name+":E",'E',components,weights,f,n_scores,buffer_factor) );
			storage_[i].push_back( pool );
			if ( trQuotaCollector.Trace.visible() ) {
				trQuotaCollector.Trace << "Pool >"<<name<<":E< added at query pos: "<<i<<" as number: "<<storage_[i].size()<<std::endl;
				trQuotaCollector.Trace<<"Scoring scheme for the quota pool sorting is:";
				for ( Size l=1; l<=weights.size(); l++ ) {
					trQuotaCollector.Trace<<"\n\t"<<components[l]<<"\t"<<weights[l];
				}
				trQuotaCollector.Trace<<std::endl;
			}
		}

		f = prediction->loop_fraction(i+middle-1)*prediction_fraction;
		n = (Size) (n_candidates * f );
		if ( n == 0 ) {
			trQuotaCollector.Debug << "Pool >"<<name<<"< for L would have size 0, not created"<<std::endl;
		} else {
			SecondaryStructurePoolOP pool( new SecondaryStructurePool(n_candidates,name+":L",'L',components,weights,f,n_scores,buffer_factor) );
			storage_[i].push_back( pool );
			if ( trQuotaCollector.Trace.visible() ) {
				trQuotaCollector.Trace << "Pool >"<<name<<":E< added at query pos: "<<i<<" as number: "<<storage_[i].size()<<std::endl;
				trQuotaCollector.Trace<<"Scoring scheme for the quota pool sorting is:";
				for ( Size l=1; l<=weights.size(); l++ ) {
					trQuotaCollector.Trace<<"\n\t"<<components[l]<<"\t"<<weights[l];
				}
				trQuotaCollector.Trace<<std::endl;
			}
		}
	}
}

} // quota
} // frag_picker
} // protocols

