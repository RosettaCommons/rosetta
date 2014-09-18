// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/quota/TorsionBinPool.hh
/// @brief provides a quota selector based on torsion bin prediciton
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/TorsionBinPool.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>

// utility headers
#include <utility/vector1.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>

namespace protocols {
namespace frag_picker {
namespace quota {


static thread_local basic::Tracer trTorsionBinPool(
		"protocols.frag_picker.quota.TorsionBinPool");

SecondaryStructurePool::~TorsionBinPool() {}

/// @brief prints information on which fragments can be accepted by this pool and how many of them
void TorsionBinPool::show_availability(std::ostream & where) const {
	 where << std::setw(10) << get_pool_name()<<" : nH = "<< max_h_ - nh_<<" nE = "<<max_e_ - ne_<<" nL = "<<max_l_ - nl_<<"\n";
}

bool TorsionBinPool::try_fragment(ScoredCandidate & candidate) {

  char ss = candidate.first->get_middle_ss();
  switch(ss) {
    case 'H' :
      if (nh_ < max_h_) {
        nh_++;
	++inserted_;
	return true;
      }
      return false;
    case 'E' :
      if (ne_ < max_e_) {
        ne_++;
	++inserted_;
	return true;
      }
      return false;
    default:
      if (nl_ < max_l_) {
        nl_++;
	++inserted_;
	return true;
      }
      return false;
  }
}

void TorsionBinPool::restart(Size query_seq_pos,Size fragment_size) {

    frag_size_ = fragment_size;
    Size query_center_frag_pos = query_seq_pos + fragment_size / 2;
    max_h_ = round( total_size_ * prediction_->helix_fraction( query_center_frag_pos ) );
    max_e_ = round( total_size_ * prediction_->strand_fraction( query_center_frag_pos ) );
    max_l_ = round( total_size_ * prediction_->loop_fraction( query_center_frag_pos ) );
    nh_ = ne_ = nl_ = 0;
}

} // quota
} // frag_picker
} // protocols

