// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/CompositeFragmentSelector.cc
/// @brief provides a selector that picks fragments diferent enough from the fragments that have been already picked
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/DiversifyCrmsdSelector.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <core/pose/Pose.hh>

#include <core/id/NamedAtomID.hh>
//#include <core/io/pdb/pose_io.hh>
#include <numeric/model_quality/rms.hh>
#include <basic/Tracer.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace frag_picker {

using namespace core;

static basic::Tracer trDiversifyCrmsdSelector(
                "protocols.frag_picker.DiversifyCrmsdSelector");

void DiversifyCrmsdSelector::copy_coordinates(FragmentCandidateOP src, ObjexxFCL::FArray2D_double & dst) {

    pose::PoseOP pose = src->get_chunk()->get_pose();
    Size len = src->get_length();
    Size offset = src->get_first_index_in_vall() - 1;

    for (core::Size i = 1; i <= len; i++) {
	id::NamedAtomID idCA("CA", i+offset);
	PointPosition const& xyz = pose->xyz(idCA);
	for (core::Size d = 1; d <= 3; ++d) {
	    dst(d, i) = xyz[d - 1];
	}
    }
}


/// @brief  Selects desired number of fragments from a given set of candidates
void DiversifyCrmsdSelector::select_fragments(
   ScoredCandidatesVector1 const& in,
	 ScoredCandidatesVector1& out )
{

	if(in.size()==0) return;

	Size len = in[1].first->get_length();

	if ((Size) fi_.size2() < len) {
		fj_.redimension(3, len, 0.0);
		fi_.redimension(3, len, 0.0);
	}

	out.push_back( in[1] );
	for(Size i=2;i<=in.size();i++) {
		if(out.size() >= frags_per_pos() ) break;
		bool is_ok = true;
		copy_coordinates(in[i].first,fi_);
		for(Size j=1;j<=out.size();j++) {
		    copy_coordinates(out[j].first,fj_);
		    Real rms = numeric::model_quality::rms_wrapper(len,fi_,fj_);
		    if(rms<cutoff_) {
			is_ok = false;
			trDiversifyCrmsdSelector.Trace<<"Crmsd is "<<rms<<" fragment "<< *in[i].first<<" denied"<<std::endl;;
			break;
		    }
		}
		if(is_ok) {
		    out.push_back( in[i] );
		    trDiversifyCrmsdSelector.Trace<<"Fragment "<< *in[i].first<<" passed"<<std::endl;;
		}
	}
	trDiversifyCrmsdSelector<<out.size()<<" fragments passed through DiversifyCrmsdSelector at query position "
	    << in[1].first->get_first_index_in_query()<<std::endl;
}

} // frag_picker
} // protocols
