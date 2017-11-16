// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/CompositeFragmentSelector.cc
/// @brief provides a selector that picks fragments diferent enough from the fragments that have been already picked
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/DiversifyCrmsdByClustering.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>
#include <core/pose/Pose.hh>

#include <core/id/NamedAtomID.hh>
#include <numeric/model_quality/rms.hh>
#include <basic/Tracer.hh>

// Clustering stuff
#include <numeric/ClusteringTreeNode.hh>
#include <numeric/agglomerative_hierarchical_clustering.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace frag_picker {

using namespace core;

static basic::Tracer trDiversifyCrmsdByClustering(
	"protocols.frag_picker.DiversifyCrmsdByClustering");

void DiversifyCrmsdByClustering::copy_coordinates(FragmentCandidateOP src, ObjexxFCL::FArray2D_double & dst) {

	pose::PoseOP pose = src->get_chunk()->get_pose();
	core::Size len = src->get_length();
	core::Size offset = src->get_first_index_in_vall() - 1;
	for ( core::Size i = 1; i <= len; i++ ) {
		id::NamedAtomID idCA("CA", i+offset);
		PointPosition const& xyz = pose->xyz(idCA);
		for ( core::Size d = 1; d <= 3; ++d ) {
			dst(d, i) = xyz[d - 1];
		}
	}
}


/// @brief  Selects desired number of fragments from a given set of candidates
void DiversifyCrmsdByClustering::select_fragments(
	ScoredCandidatesVector1 const& in,
	ScoredCandidatesVector1& out )
{
	if ( in.size()==0 ) return;

	//-------------- core::Size of fragments
	core::Size len = in[1].first->get_length();
	trDiversifyCrmsdByClustering.Debug << "Diversifying fragments of size "<<len<<" #(in) = "<<in.size()<<std::endl;
	//-------------- Resize container for xyz data
	if ( in.size() > xyz_.size() ) {
		trDiversifyCrmsdByClustering.Trace << "core::Reallocated: to "<<xyz_.size()<<std::endl;
		xyz_.resize(in.size());
	}
	//-------------- Resize each xyz vector to match fragments' size
	for ( core::Size i=1; i<=in.size(); i++ ) {
		if ( (core::Size)xyz_[i].size2() != len ) {
			xyz_[i].redimension(3,len,0.0);
		}
	}
	//-------------- Copy xyz coordinates, setup ids
	utility::vector1<core::Size> ids;
	for ( core::Size i=1; i<=in.size(); i++ ) {
		copy_coordinates(in[i].first,xyz_[i]);
		ids.push_back(i);
	}

	//-------------- Prepare distance matrix
	utility::vector1< utility::vector1<core::Real> > distances(in.size());
	for ( core::Size i=1; i<=in.size(); i++ ) {
		distances[i].resize( in.size() );
	}

	for ( core::Size i=2; i<=in.size(); i++ ) {
		for ( core::Size j=2; j<=in.size(); j++ ) {
			core::Real val = numeric::model_quality::rms_wrapper(len,xyz_[i],xyz_[j]);
			distances[ i ][ j ] = val;
			distances[ j ][ i ] = val;
		}
		distances[i][i] = 0.0;
	}
	numeric::AverageLinkClusterer alc;
	utility::vector1<numeric::ClusteringTreeNodeOP> roots = alc.cluster(distances,frags_per_pos());

	//----------- Retrieve clusters
	for ( core::Size i=1; i<=frags_per_pos(); i++ ) {
		utility::vector1<core::Size> ids_out;
		numeric::get_cluster_data(ids,roots[i],ids_out);
		out.push_back( in[ ids_out[1] ] );
		trDiversifyCrmsdByClustering.Debug << "Fragment "<<i<<" cluster stats: dist="<<roots[i]->distance()
			<<" size="<<roots[i]->size()<<std::endl;
	}

	trDiversifyCrmsdByClustering.Trace<<out.size()<<" fragments passed through DiversifyCrmsdByClustering at query position "
		<< in[1].first->get_first_index_in_query()<<std::endl;
}

} // frag_picker
} // protocols
