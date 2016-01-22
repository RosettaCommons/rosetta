// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/sewing/conformation/ContinuousAssembly.cc
///
/// @author Tim Jacobs
/// @author Doonam Kim (append_model which is independent of number of secondary structure)

//Unit headers
//#include <protocols/sewing/conformation/Assembly.hh>
#include <protocols/sewing/conformation/ContinuousAssembly.hh>

//Package headers
#include <protocols/sewing/util/io.hh>

//Protocol headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/rotamer_set/AddResiduesRotamerSetOperation.hh>
#include <core/scoring/constraints/ResidueTypeConstraint.hh>

#include <core/conformation/Residue.functions.hh>
#include <core/chemical/VariantType.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

#include <protocols/loops/Loop.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzTransform.hh>

#include <utility/LexicographicalIterator.hh>
#include <utility/vector1.hh>

//Numeric headers


namespace protocols {
namespace sewing  {

static basic::Tracer TR("protocols.sewing.ContinuousAssembly");

ContinuousAssembly::ContinuousAssembly():
	Assembly()
{}

AssemblyOP
ContinuousAssembly::clone(){
	return AssemblyOP( new ContinuousAssembly(*this) );
}


/// @details add a model to assembly
void
ContinuousAssembly::append_model(
	Model const & model,
	ScoreResult const & edge_score
){

	utility::vector1<core::Size> matched_indices;

	std::map<SewSegment,SewSegment> matching_segments =
		get_matching_model_segments(model, edge_score); //mobile_model, edge_score
	std::set<SewSegment> mobile_match_segments;
	std::map<SewSegment,SewSegment>::const_iterator it = matching_segments.begin();
	std::map<SewSegment,SewSegment>::const_iterator it_end = matching_segments.end();
	for ( ; it != it_end; ++it ) {
		SewSegment ref_seg = it->first; // ref_seg from reference (existing) model
		SewSegment mobile_seg = it->second; // mobile_seg from new model

		mobile_match_segments.insert(it->second);

		for ( core::Size i=1; i<=all_segments_.size(); ++i ) {
			for ( core::Size j=1; j<=all_segments_[i].size(); ++j ) {
				if ( all_segments_[i][j] == ref_seg ) {
					all_segments_[i].push_back(mobile_seg);
					matched_indices.push_back(i);
				}
			}
		}
	}//for(; it != it_end; ++it) {

	//Generate the chimera segments
	utility::vector1<SewSegment> chimeras = get_chimera_segments(matching_segments, edge_score.second.segment_matches, model); //matching_segments, segment_matches,                   mobile_model
	runtime_assert(matched_indices.size() == chimeras.size());
	// if(matched_indices.size() != chimeras.size()) {
	//  TR.Warning << "Matched indices isn't the same size as chimeras!" << std::endl;
	//  TR.Warning << "Matched indices (" << matched_indices.size() << "): " << std::endl;
	//  for(core::Size i=1; i<=matched_indices.size(); ++i){
	//   TR.Warning << "\t" << i << "-" << matched_indices[i] << std::endl;
	//  }
	//
	//  TR.Warning << "Chimeras (" << chimeras.size() << "): " << std::endl;
	//  for(core::Size i=1; i<=chimeras.size(); ++i) {
	//   TR.Warning << "\t" << i << " - " << chimeras[i].model_id_ << " " << chimeras[i].segment_id_ << std::endl;
	//  }
	//
	//  TR.Warning << "Matching segments (" << matching_segments.size() << "): " << std::endl;
	//  std::map<SewSegment, SewSegment>::const_iterator it = matching_segments.begin();
	//  std::map<SewSegment, SewSegment>::const_iterator it_end = matching_segments.end();
	//  for(; it != it_end; ++it) {
	//   TR.Warning << it->first.model_id_ << " " << it->first.segment_id_
	//     << " -> " << it->second.model_id_ << " " << it->second.segment_id_ << std::endl;
	//  }
	//
	//  int ref_model_id = matching_segments.begin()->first.model_id_;
	//  Model test_model = regenerate_model(ref_model_id);
	//  std::map<int, Model> test_model_map;
	//  test_model_map.insert(std::make_pair(ref_model_id, test_model));
	//  write_model_file(test_model_map, "regenerated.models");
	//
	//  utility_exit_with_message("Failed to properly generate chimeras");
	// }

	for ( core::Size i = 1; i <= chimeras.size(); ++i ) {
		all_segments_[matched_indices[i]].push_back(chimeras[i]);
	}

	//Replace the assembly segments that were matched on with the new chimera segments
	core::Size j=1;
	it = matching_segments.begin();
	for ( ; it != it_end; ++it ) {
		for ( core::Size i=1; i<=segments_.size(); ++i ) {
			if ( segments_[i] == it->first ) { // 'it->first' could represent chimera segment
				segments_[i] = chimeras[j];
			}
		}
		++j;
	}

	// Doonam changed to deal with any number of segments (not only limited to smotif)
	int segments_in_model_size = model.segments_.size();
	//Check to see if we matched on the first or last segment
	if ( mobile_match_segments.find(model.segments_[1]) == mobile_match_segments.end() ) {
		//Appending to n_terminus
		utility::vector1<SewSegment> add_these_to_front;

		for ( int seg_to_add_i=1; seg_to_add_i < segments_in_model_size; seg_to_add_i++ ) {
			add_these_to_front.push_back(model.segments_[seg_to_add_i]);
		}

		segments_.insert(segments_.begin(), add_these_to_front.begin(), add_these_to_front.end());

		//Update all_segments vector
		for ( int seg_to_add_i = segments_in_model_size-1; seg_to_add_i > 0; seg_to_add_i-- ) {
			utility::vector1<SewSegment> model_seg_vec;
			model_seg_vec.push_back(model.segments_[seg_to_add_i]);
			all_segments_.insert(all_segments_.begin(), model_seg_vec);
		}

		if ( TR.Debug.visible() ) { TR.Debug << "Model added to the N-terminal" << std::endl; }
	} else if ( mobile_match_segments.find(model.segments_.back()) == mobile_match_segments.end() ) {
		//Appending to c_terminus
		for ( int seg_to_add_i=2; seg_to_add_i < segments_in_model_size+1; seg_to_add_i++ ) {
			segments_.push_back(model.segments_[seg_to_add_i]);
		}

		//Update all_segments vector
		for ( int seg_to_add_i=2; seg_to_add_i < segments_in_model_size+1; seg_to_add_i++ ) {
			utility::vector1<SewSegment> model_seg_vec;
			model_seg_vec.push_back(model.segments_[seg_to_add_i]);
			all_segments_.push_back(model_seg_vec);
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "Model added to the C-terminal" << std::endl;
		}
	} else {
		utility_exit_with_message("Didn't match on the first or last segment! OH NO!");
	}

	//erase old connections
	segments_.clear_connections();

	//go through and add new connections from n->c
	for ( core::Size i=1; i<=segments_.size() - 1; ++i ) {
		segments_.add_connection(i, i+1);
	}

	for ( core::Size i=1; i<=segments_.size(); ++i ) {
		TR.Debug << "Model_id: " << segments_[i].model_id_ << ", Segment_ID: " << segments_[i].segment_id_ << std::endl;
		TR.Debug << "Has previous " << segments_.has_previous(i) << ", Has next " << segments_.has_next(i) << std::endl;
	}

} //append_model


} //sewing namespace
} //protocols namespace
