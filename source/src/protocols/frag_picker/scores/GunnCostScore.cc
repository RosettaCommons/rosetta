// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/GunnCostScore.cc
/// @brief  Object that scores a fragment by its crmsd to the native
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/GunnCostScore.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/FragmentPicker.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <protocols/frag_picker/scores/GunnCost.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer trGunnScore(
	"protocols.frag_picker.scores.GunnCostScore");

GunnCostScore::~GunnCostScore() {}

GunnCostScore::GunnCostScore(Size priority, Real lowest_acceptable_value,
	bool use_lowest, core::pose::PoseOP reference_pose,
	utility::vector1<Size> frag_sizes,Size max_chunk_size) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "GunnCostScore") , gunn_cost_( -0.01 ) {
	reference_pose_ = reference_pose;
	n_atoms_ = reference_pose_->size();
	max_chunk_size_ = max_chunk_size;

	for ( Size j=1; j<=frag_sizes.size(); j++ ) {
		frag_sizes_.push_back(frag_sizes[j]);
		trGunnScore << "Preparing Gunn score for fragments size: "<<frag_sizes[j]<<std::endl;

		utility::vector1<GunnTuple> v;
		v.resize( max_chunk_size_-frag_sizes[j]+1 );
		chunk_gunn_data_.push_back( v );
		trGunnScore << "Prepared "<<v.size()<<" Gunn tuples for the vall structure"<<std::endl;

		utility::vector1<GunnTuple> vn;
		vn.resize( reference_pose_->size()-frag_sizes[j]+1 );
		computeGunnTuples(*reference_pose_,frag_sizes[j],vn);
		ref_gunn_data_.push_back( vn );
		trGunnScore << "Prepared "<<vn.size()<<" Gunn tuples for the reference structure"<<std::endl;
	}
}


bool GunnCostScore::score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	core::pose::PoseOP pose = f->get_chunk()->get_pose();

	Size iv = f->get_first_index_in_vall();
	Size iq = f->get_first_index_in_query();
	GunnTuple t1;
	GunnTuple t2;
	gunn_cost_.compute_gunn( *pose, iv, iv + f->get_length() -1, t1);
	gunn_cost_.compute_gunn( *reference_pose_, iq, iq + f->get_length() -1, t2);

	Real score = gunn_cost_.score_tuple(t1,t2);
	empty_map->set_score_component(score, id_);
	if ( (score > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}

	return true;
}


void GunnCostScore::do_caching(VallChunkOP current_chunk) {

	core::pose::PoseOP pose = current_chunk->get_pose();
	for ( Size j=1; j<=frag_sizes_.size(); j++ ) {
		trGunnScore.Debug <<"Caching "<<frag_sizes_[j]<<" buffer size: "<<chunk_gunn_data_[j].size()<<std::endl;
		computeGunnTuples(*pose,frag_sizes_[j],chunk_gunn_data_[j]);
	}
}

void GunnCostScore::computeGunnTuples(core::pose::Pose & pose,Size frag_size,utility::vector1<GunnTuple> & result) {

	trGunnScore.Debug << "Computing Gunn tuples for the vall structure of size: "
		<<pose.size()<<", results buffer size is: "<<result.size()<<std::endl;
	for ( Size i=1; i<=pose.size()-frag_size + 1; i++ ) {
		gunn_cost_.compute_gunn( pose, i, i+frag_size-1, result[i]);
	}
}


bool GunnCostScore::cached_score(FragmentCandidateOP fragment,
	FragmentScoreMapOP scores) {

	Size frag_len = fragment->get_length();
	Size ifr = 1;
	//bool if_found = false;
	for ( ifr=1; ifr<=frag_sizes_.size(); ifr++ ) {
		if ( frag_sizes_[ifr] == frag_len ) {
			//if_found = true;  set but never used ~Labonte
			break;
		}
	}
	// if(! if_found) {
	//     frag_sizes_.push_back( frag_len );
	//     ifr = frag_sizes_.size() + 1;
	//     prepare_chunk_tuples(frag_len);
	// }

	Size iv = fragment->get_first_index_in_vall();
	Size iq = fragment->get_first_index_in_query();

	Real score = gunn_cost_.score_tuple(chunk_gunn_data_[ifr][iv],ref_gunn_data_[ifr][iq]);
	GunnTuple & t = ref_gunn_data_[ifr][iq];
	trGunnScore.Trace<<t.q1<<" ";
	trGunnScore.Trace<<t.q2<<" ";
	trGunnScore.Trace<<t.q3<<" ";
	trGunnScore.Trace<<t.q4<<" ";
	trGunnScore.Trace<<t.q5<<" ";
	trGunnScore.Trace<<t.q6<<std::endl;
	GunnTuple & t2 = chunk_gunn_data_[ifr][iv];
	trGunnScore.Trace<<t2.q1<<" ";
	trGunnScore.Trace<<t2.q2<<" ";
	trGunnScore.Trace<<t2.q3<<" ";
	trGunnScore.Trace<<t2.q4<<" ";
	trGunnScore.Trace<<t2.q5<<" ";
	trGunnScore.Trace<<t2.q6<<std::endl;
	trGunnScore.Trace<<"scoring: "<<iv<<", "<<iq<<" "<<gunn_cost_.score_tuple(t2,t)<<" "<<score<<std::endl;

	scores->set_score_component(score, id_);
	if ( (score > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}

	return true;
}

void GunnCostScore::clean_up() {
}

FragmentScoringMethodOP MakeGunnCostScore::make(Size priority,
	Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker
	, std::string) {


	Size longest_chunk = picker->get_vall()->get_largest_chunk_size();

	utility::vector1<Size> frag_sizes_tmp;
	if ( option[frags::frag_sizes].user() ) {
		frag_sizes_tmp = option[frags::frag_sizes]();
	} else {
		frag_sizes_tmp.push_back(3);
		frag_sizes_tmp.push_back(9);
	}

	if ( option[in::file::s].user() ) {
		trGunnScore
			<< "Reference structure to score fragments by Gunn cost loaded from: "
			<< option[in::file::s]()[1] << std::endl;
		core::pose::PoseOP nativePose( new core::pose::Pose );
		core::import_pose::pose_from_file(*nativePose, option[in::file::s]()[1], core::import_pose::PDB_file);

		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new GunnCostScore(priority,
			lowest_acceptable_value, use_lowest, nativePose,
			frag_sizes_tmp,longest_chunk) );
	}
	utility_exit_with_message(
		"Can't read a reference structure. Provide it with in::file::s flag");

	return NULL;
}

} // scores
} // frag_picker
} // protocols

