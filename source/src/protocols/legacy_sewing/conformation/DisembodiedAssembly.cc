// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/legacy_sewing/conformation/Assembly.cc
///
/// @author Tim Jacobs

//Unit headers
#include <protocols/legacy_sewing/conformation/Assembly.hh>
#include <protocols/legacy_sewing/conformation/DisembodiedAssembly.hh>

//Package headers
#include <protocols/legacy_sewing/util/io.hh>

//Protocol headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_set/AddResiduesRotamerSetOperation.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>
#include <core/scoring/hbonds/HBondSet.hh>

#include <core/conformation/Residue.functions.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <protocols/loops/Loop.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/legacy_sewing.OptionKeys.gen.hh>

#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/xyzTransform.hh>

#include <utility/LexicographicalIterator.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/vector1.hh>

//Numeric headers


namespace protocols {
namespace legacy_sewing  {

static THREAD_LOCAL basic::Tracer TR("protocols.legacy_sewing.DisembodiedAssembly");

DisembodiedAssembly::DisembodiedAssembly():
	Assembly()
{}


AssemblyOP
DisembodiedAssembly::clone(){
	return AssemblyOP(new DisembodiedAssembly(*this));
}


void
DisembodiedAssembly::append_model(
	Model const & model,
	ScoreResult const & edge_score
){

	utility::vector1<core::Size> matched_indices;

	std::map<SewSegment,SewSegment> matching_segments =
		get_matching_model_segments(model, edge_score);

	std::set<SewSegment> mobile_match_segments;
	std::map<SewSegment,SewSegment>::const_iterator it = matching_segments.begin();
	std::map<SewSegment,SewSegment>::const_iterator it_end = matching_segments.end();
	for ( ; it != it_end; ++it ) {
		SewSegment ref_seg = it->first;
		SewSegment mobile_seg = it->second;

		mobile_match_segments.insert(it->second);

		//  matched_segments_.insert(ref_seg);
		//  matched_segments_.insert(mobile_seg);
		for ( core::Size i=1; i<=all_segments_.size(); ++i ) {
			for ( core::Size j=1; j<=all_segments_[i].size(); ++j ) {
				if ( all_segments_[i][j] == ref_seg ) {
					all_segments_[i].push_back(mobile_seg);
					matched_indices.push_back(i);
				}
			}
		}
	}

	// utility::vector1<SewSegment> chimeras = get_chimera_segments(matching_segments, edge_score.second.segment_matches);
	// runtime_assert(matched_indices.size() == chimeras.size());
	// for (core::Size i = 1; i <= chimeras.size(); ++i ) {
	//  all_segments_[matched_indices[i]].push_back(chimeras[i]);
	// }
	//
	// //Replace the assembly segments that were matched on with the new chimera segments
	// core::Size j=1;
	// it = matching_segments.begin();
	// for(; it != it_end; ++it) {
	//  for(core::Size i=1; i<=segments_.size(); ++i) {
	//   if(segments_[i] == it->first) {
	//    segments_[i] = chimeras[j];
	//   }
	//  }
	//  ++j;
	// }

	//now go through all the new model segments and add any that
	//weren't matched and aren't loops connected to any of the matched segments
	//It's possible in the case of directly adjacent segments that a connected segment is
	//not necessarily a loop, so check the 'hash_' field to ensure we aren't removing valid
	//segments
	utility::vector1< std::pair<core::Size, core::Size> > model_to_assembly_indices;
	for ( core::Size i=1; i<=model.segments_.size(); ++i ) {
		if ( mobile_match_segments.find(model.segments_[i]) == mobile_match_segments.end() ) {
			if ( !model.segments_[i].hash_ && model.segments_.has_next(i) && mobile_match_segments.find(model.segments_[model.segments_.next(i)]) != mobile_match_segments.end() ) {
				if ( TR.Debug.visible() ) {
					TR.Debug << "Discarding loop segment " << model.segments_[i].segment_id_ << " as it precedes the matched segment" << std::endl;
				}
				continue;
			} else if ( !model.segments_[i].hash_ && model.segments_.has_previous(i) && mobile_match_segments.find(model.segments_[model.segments_.previous(i)]) != mobile_match_segments.end() ) {
				if ( TR.Debug.visible() ) {
					TR.Debug << "Discarding loop segment " << model.segments_[i].segment_id_ << " as it follows the matched segment" << std::endl;
				}
				continue;
			} else {
				//We're keeping the segment, so add it to the segments_ vector. Since this
				//is the first segment at this 'position', create a vector of segments
				//(currently containing only this segment) and add it to the all_segments_
				//vector
				segments_.push_back(model.segments_[i]);
				utility::vector1<SewSegment> model_seg_vec;
				model_seg_vec.push_back(model.segments_[i]);
				all_segments_.push_back(model_seg_vec);
				model_to_assembly_indices.push_back(std::make_pair(i, segments_.size()));
				if ( TR.Debug.visible() ) {
					TR.Debug << "Appending model segment " << model.segments_[i].segment_id_ << " to Assembly" << std::endl;
				}
			}
		}
	}

	//Finally. do up the connections
	for ( core::Size i=1; i<=model_to_assembly_indices.size(); ++i ) {
		for ( core::Size j=i+1; j<=model_to_assembly_indices.size(); ++j ) {
			if ( model.segments_.has_next(model_to_assembly_indices[i].first) &&
					model.segments_.next(model_to_assembly_indices[j].first) ) {
				segments_.add_connection(model_to_assembly_indices[i].second, model_to_assembly_indices[j].second);
				if ( TR.Debug.visible() ) {
					TR.Debug << "Copying model connection " << model_to_assembly_indices[i].first << "-" << model_to_assembly_indices[j].first
						<< ". New connection between assembly segments " << model_to_assembly_indices[i].second << "-" << model_to_assembly_indices[j].second << std::endl;
				}
			}
		}
	}
}//DisembodiedAssembly::append_model


} //legacy_sewing namespace
} //protocols namespace
