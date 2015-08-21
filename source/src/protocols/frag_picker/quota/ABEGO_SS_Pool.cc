// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/quota/ABEGO_SS_Pool.hh
/// @brief provides a quota selector based on difefrent secondary structure predicitons
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/quota/ABEGO_SS_Pool.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>
#include <protocols/frag_picker/CommonFragmentComparators.hh>

// utility headers
#include <utility/vector1.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

//C++
#include <sstream>

namespace protocols {
namespace frag_picker {
namespace quota {

static thread_local basic::Tracer trABEGO_SS_Pool(
	"protocols.frag_picker.quota.ABEGO_SS_Pool");

/// @brief Creates a pool of a given size and name
/// @param total_size - how many fragments will be selected (in total in all pools)
/// @param name - name assigned to this pool. This in general may be any string that
/// later allows one control pool's behavior from a flag file
/// @param ss_type - what is the type of secondary structur this pool is accepting
/// @param fraction - fraction of the total number of fragments that goes into this pool
ABEGO_SS_Pool::ABEGO_SS_Pool(Size total_size,std::string pool_name,
	utility::vector1< std::pair<Size,Size> > ss_abego_types,
	utility::vector1<Size> which_components,utility::vector1<Real> weights,Real fraction,Size n_scores,Size buffer_factor = 5) :
	QuotaPool(pool_name,fraction) {

	ss_abego_types_ = ABEGO_SS_MapOP( new ABEGO_SS_Map(ss_abego_types) );

	buffer_factor_ = buffer_factor;
	debug_assert ( which_components.size() == weights.size() );
	for ( Size i=1; i<=which_components.size(); i++ ) {
		components_.push_back( which_components[i] );
		weights_.push_back( weights[i] );
	}
	total_size_ = total_size;
	this_size_ = (Size)(fraction * total_size);
	if ( this_size_<20 ) {
		this_size_ = 20;
	}

	CompareByScoreCombination ordering(components_, weights_);
	storage_ = BoundedQuotaContainerOP( new BoundedQuotaContainer(ordering,this_size_,this_size_*buffer_factor) );
	FragmentCandidateOP worst_f( new FragmentCandidate(1,1,0,1) );
	scores::FragmentScoreMapOP worst_s( new scores::FragmentScoreMap(n_scores) );
	storage_->set_worst(std::pair<FragmentCandidateOP, scores::FragmentScoreMapOP>(worst_f,worst_s));
	for ( Size i=1; i<=n_scores; i++ ) {
		worst_s->set_score_component(99999.9999,i);
	}

	if ( trABEGO_SS_Pool.Debug.visible() ) {
		trABEGO_SS_Pool.Debug<< "Creating a pool >"<<pool_name<<"< for: \n";
		//     for(Size i=1;i<=ss_abego_types_.size();i++) {
		//  trABEGO_SS_Pool.Debug<<ss_abego_types_[i].first<<":"<<ss_abego_types_[i].second<<"\n";
		//     }
		trABEGO_SS_Pool.Debug<< ss_abego_types_->show_valid();
		trABEGO_SS_Pool.Debug<<"\nholding "<<this_size_<<" candidates, quota fraction is: "<<fraction<<std::endl;
	}
}


ABEGO_SS_Pool::~ABEGO_SS_Pool() {}

bool ABEGO_SS_Pool::could_be_accepted(ScoredCandidate candidate) const {

	VallResidueOP r = candidate.first->get_middle_residue();
	Size abego_bin = torsion2big_bin_id(r->phi(),r->psi(),r->omega());

	return ss_abego_types_->check_status(candidate.first->get_middle_ss(),abego_bin);
}


bool ABEGO_SS_Pool::add(ScoredCandidate candidate) {

	if ( could_be_accepted(candidate) ) {
		return storage_->push( candidate );
	}
	return false;
}

void ABEGO_SS_Pool::print_report(
	std::ostream & out,
	scores::FragmentScoreManagerOP //manager
) const
{
	out << get_pool_name() << " collected "<<current_size()<<std::endl;
}

} // quota
} // frag_picker
} // protocols

