// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/quota/QuotaSelector.hh
/// @brief provides a quota selector based on difefrent secondary structure predicitons
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/quota/QuotaSelector.hh>
#include <protocols/frag_picker/quota/QuotaPool.hh>
#include <protocols/frag_picker/FragmentSelectingRule.hh>
#include <protocols/frag_picker/FragmentPicker.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>

// utility headers
#include <core/types.hh>
#include <basic/Tracer.hh>

// C++
#include <algorithm>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace quota {

static THREAD_LOCAL basic::Tracer trQuotaSelector(
	"protocols.frag_picker.quota.QuotaSelector");

QuotaSelector::QuotaSelector(Size frags_per_pos,Size pos_in_query,QuotaCollectorOP collector) :
	FragmentSelectingRule(frags_per_pos),collector_(collector) {

	trQuotaSelector.Trace << "Selector set up for pos. "<<pos_in_query<<std::endl;
	q_pos_ = pos_in_query;
	collector_ = collector;
}


Size QuotaSelector::next_from_pool( ScoredCandidatesVector1 const& pool,Size recently_taken,std::set<Size> & in_use) {
	if ( recently_taken >= pool.size() ) {
		return 0;
	}
	recently_taken++;

	for ( Size i=recently_taken; i<=pool.size(); i++ ) {
		if ( pool[i].first == 0 ) {
			trQuotaSelector.Warning << "Fragment candidate at pos. " << i
				<< " from a pool is null - nothing to select\n"
				<< " maybe your vall database is to small?" << std::endl;
			return 0;
		}
		if ( in_use.find( pool[i].first->key() ) == in_use.end() ) {
			in_use.insert( pool[i].first->key() );
			return i;
		}
	}

	return 0;
}


struct quota_limits_sorter_biggest_first {
	bool operator() (std::pair<Size,Size> i,std::pair<Size,Size> j) { return (i.second>j.second); }
} mysorter_biggest_first;

struct quota_limits_sorter_smallest_first {
	bool operator() (std::pair<Size,Size> i,std::pair<Size,Size> j) { return (i.second<j.second); }
} mysorter_smallest_first;

struct missing_fraction_sorter_biggest_first {
	bool operator() (std::pair<Size,Real> i,std::pair<Size,Real> j) { return (i.second>j.second); }
} mysorter_missing_fraction_biggest_first;


void QuotaSelector::push_the_limits(utility::vector1<Size> & q_limits,Size target_total) {

	Size total = 0;
	utility::vector1< std::pair<Size,Size> > tmp;
	for ( Size i=1; i<=q_limits.size(); i++ ) {
		total += q_limits[i];
		tmp.push_back( std::pair<Size,Size>(i,q_limits[i]) );
	}
	Size diff = target_total - total;
	std::sort(tmp.begin(),tmp.end(),mysorter_biggest_first);
	trQuotaSelector.Trace << "Pushing the quota limits from total "<<total<<" to "<<target_total
		<<", the first pool has size "<<tmp[1].second<<std::endl;
	// AMW: cppcheck notes that this is just a while true, so let's make it a real one (from 2==2)
	while ( true ) {
		for ( Size i=1; i<=q_limits.size(); i++ ) {
			if ( diff==0 ) {
				return;
			}
			q_limits[ tmp[i].first ] ++;
			diff--;
		}
	}
}


void QuotaSelector::push_the_limits_to_the_winner(utility::vector1<Size> & q_limits,Size target_total) {

	Size total = 0;
	utility::vector1< std::pair<Size,Size> > tmp;
	for ( Size i=1; i<=q_limits.size(); i++ ) {
		total += q_limits[i];
		tmp.push_back( std::pair<Size,Size>(i,q_limits[i]) );
	}
	Size diff = target_total - total;
	std::sort(tmp.begin(),tmp.end(),mysorter_biggest_first);
	trQuotaSelector.Trace << "Pushing the quota limits from total "<<total<<" to "<<target_total
		<<", the first pool has size "<<tmp[1].second<<" and takes extra "<<diff<<std::endl;
	q_limits[ tmp[1].first ] += diff;
}

void QuotaSelector::push_the_limits(utility::vector1<Size> & q_limits,Size target_total,
	utility::vector1<Real> & q_fractions) {


	Size total = 0;
	for ( Size i=1; i<=q_limits.size(); i++ ) {
		total += q_limits[i];
	}
	Size diff = target_total - total;

	trQuotaSelector.Trace << "Pushing the quota limits from total "<<total<<" to "<<target_total
		<<std::endl;

	utility::vector1< std::pair<Size,Real> > tmp;
	while ( diff > 0 ) {
		tmp.clear();
		for ( Size i=1; i<=q_limits.size(); i++ ) {
			Real missing = q_fractions[i] - ((Real) q_limits[i]) / target_total;
			tmp.push_back( std::pair<Size,Real>(i,missing) );
		}
		std::sort(tmp.begin(),tmp.end(),mysorter_missing_fraction_biggest_first);
		q_limits[ tmp[1].first ] ++;
		// for(Size i=1;i<=q_limits.size();i++)
		//     trQuotaSelector.Trace << "pool id,old fraction, new fraction : "<<tmp[i].first<<" "<<tmp[i].second<<" "<<
		//  ((Real) q_limits[ tmp[i].first ]) / target_total<<std::endl;
		diff--;
	}
}


void QuotaSelector::select_fragments_200(
	ScoredCandidatesVector1 const&,
	ScoredCandidatesVector1& output_selection
) {

	utility::vector1<Size> last_selected(collector_->count_pools(q_pos_),0);
	Size n = 0;
	std::set<Size> nice_collection;
	bool something_taken = false;
	utility::vector1<Size> q_limits;
	utility::vector1<Size> q_cnt;
	Size total = 0;
	for ( Size i=1; i<=collector_->count_pools(q_pos_); ++i ) {
		q_limits.push_back( collector_->get_pool(q_pos_,i)->total_size() );
		total += collector_->get_pool(q_pos_,i)->total_size();
		q_cnt.push_back(0);
	}
	for ( Size i=1; i<=collector_->count_pools(q_pos_); ++i ) {
		if ( total == 0 ) {
			q_limits[i] = 0;
		} else {
			q_limits[i] = (q_limits[i] * frags_per_pos() / total);
		}
	}

	push_the_limits( q_limits, frags_per_pos() );
	do {
		something_taken = false;
		for ( Size i=1; i<=collector_->count_pools(q_pos_); ++i ) {
			ScoredCandidatesVector1 const& fs = collector_->get_pool(q_pos_,i)->get_candidates(0);
			if ( q_cnt[i] == q_limits[i] ) {
				continue;
			}
			Size ifrag = next_from_pool( fs, last_selected[i], nice_collection );
			if ( ifrag > 0 ) {
				output_selection.push_back( fs[ifrag] );
				n++;
				something_taken = true;
				q_cnt[i]++;
			}
		}
	} while((n < frags_per_pos()) && (something_taken==true));
	trQuotaSelector.Trace<<"Position "<<q_pos_<<" done, "<< output_selection.size()<<" fragments selected."<<std::endl;
}


void QuotaSelector::select_fragments_25_200(
	ScoredCandidatesVector1 const&,
	ScoredCandidatesVector1& output_selection
) {
	if ( collector_->count_candidates() == 0 ) {
		return;
	}
	utility::vector1<Size> last_selected(collector_->count_pools(q_pos_),0);
	Size frag_size = collector_->get_pool(q_pos_,1)->get_candidates(0)[1].first->get_length();
	Size n = 0;
	std::set<Size> nice_collection;
	bool something_taken = false;
	utility::vector1<Size> q_limits;
	utility::vector1<Size> q_cnt;
	utility::vector1<Real> q_frac;
	Real total = 0.0;
	for ( Size i=1; i<=collector_->count_pools(q_pos_); ++i ) {
		total += collector_->get_pool(q_pos_,i)->get_fraction();
	}
	total = 25.0 / total;
	for ( Size i=1; i<=collector_->count_pools(q_pos_); ++i ) {
		q_limits.push_back( (Size)(collector_->get_pool(q_pos_,i)->get_fraction() * total) );
		q_cnt.push_back(0);
		q_frac.push_back( collector_->get_pool(q_pos_,i)->get_fraction() );
	}

	push_the_limits( q_limits, 25, q_frac );
	if ( trQuotaSelector.Trace.visible() ) {
		trQuotaSelector.Trace<<"Position: "<<q_pos_<<" frag_size:"<<frag_size<<" pool allowances for the first 25 fragments are:\n";
		for ( Size i=1; i<=collector_->count_pools(q_pos_); ++i ) {
			trQuotaSelector.Trace<<collector_->get_pool(q_pos_,i)->get_pool_name()<<" "<<q_limits[i]
				<<" ("<<collector_->get_pool(q_pos_,i)->get_fraction()*25<<")\n";
		}
		trQuotaSelector.Trace<<std::endl;
	}
	do {
		something_taken = false;
		for ( Size i=1; i<=collector_->count_pools(q_pos_); ++i ) {
			if ( q_cnt[i] == q_limits[i] ) {
				continue;
			}
			QuotaPoolOP the_pool = collector_->get_pool(q_pos_,i);
			ScoredCandidatesVector1 const& fs = the_pool->get_candidates(0);
			Size ifrag = next_from_pool( fs, last_selected[i], nice_collection );
			if ( ifrag > 0 ) {
				fs[ifrag].second->set_quota_score( the_pool->quota_score( fs[ifrag] ) );
				fs[ifrag].first->set_pool_name( the_pool->get_pool_name() );

				output_selection.push_back( fs[ifrag] );
				n++;
				something_taken = true;
				q_cnt[i]++;
			}
		}
	} while((n < 25) && (something_taken==true));

	utility::vector1<Real> tmp(FragmentPicker::log_25_.max_pools(),1000000000000.0);
	for ( Size i=1; i<=collector_->count_pools(q_pos_); ++i ) {
		QuotaPoolOP p = collector_->get_pool(q_pos_,i);
		tmp[FragmentPicker::log_25_.tag_map_[p->get_pool_name()]] = (q_limits[i] - p->get_fraction()*25);
	}
	FragmentPicker::log_25_.log(frag_size,q_pos_,tmp);

	Size fr_per_pos = frags_per_pos();
	total = 0.0;
	for ( Size i=1; i<=collector_->count_pools(q_pos_); ++i ) {
		total += collector_->get_pool(q_pos_,i)->get_fraction();
	}
	total = ((Real) fr_per_pos) / total;
	for ( Size i=1; i<=collector_->count_pools(q_pos_); ++i ) {
		q_limits[i] =  (Size)(collector_->get_pool(q_pos_,i)->get_fraction() * total);
	}
	push_the_limits( q_limits, fr_per_pos , q_frac);
	if ( trQuotaSelector.Trace.visible() ) {
		trQuotaSelector.Trace<<"Position: "<<q_pos_<<" frag_size:"<<frag_size<<" pool allowances are:\n";
		for ( Size i=1; i<=collector_->count_pools(q_pos_); ++i ) {
			trQuotaSelector.Trace<<collector_->get_pool(q_pos_,i)->get_pool_name()<<" "<<q_limits[i]
				<<" ("<<collector_->get_pool(q_pos_,i)->get_fraction()*fr_per_pos<<")\n";
		}
		trQuotaSelector.Trace<<std::endl;
	}
	do {
		something_taken = false;
		for ( Size i=1; i<=collector_->count_pools(q_pos_); ++i ) {
			if ( q_cnt[i] == q_limits[i] ) {
				continue;
			}
			QuotaPoolOP the_pool = collector_->get_pool(q_pos_,i);
			ScoredCandidatesVector1 const& fs = the_pool->get_candidates(0);
			Size ifrag = next_from_pool( fs, last_selected[i], nice_collection );
			last_selected[i] = ifrag;
			if ( ifrag > 0 ) {
				fs[ifrag].second->set_quota_score( the_pool->quota_score( fs[ifrag] ) );
				fs[ifrag].first->set_pool_name( the_pool->get_pool_name() );
				output_selection.push_back( fs[ifrag] );
				n++;
				something_taken = true;
				q_cnt[i]++;
			}
		}
	} while((n < fr_per_pos) && (something_taken==true));

	utility::vector1<Real> tmp2(FragmentPicker::log_200_.max_pools(),1000000000000.0);
	for ( Size i=1; i<=collector_->count_pools(q_pos_); ++i ) {
		QuotaPoolOP p = collector_->get_pool(q_pos_,i);
		tmp2[FragmentPicker::log_200_.tag_map_[p->get_pool_name()]] = (q_limits[i] - p->get_fraction()*fr_per_pos);
	}
	FragmentPicker::log_200_.log(frag_size,q_pos_,tmp2);

	// trQuotaSelector.Debug<<"Position "<<q_pos_<<" done, "<< output_selection.size()<<" fragments selected."<<std::endl;
}


} // quota
} // frag_picker
} // protocols
