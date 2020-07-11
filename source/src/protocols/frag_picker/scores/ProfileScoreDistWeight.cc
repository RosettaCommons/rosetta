// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/ProfileScore.cc
/// @brief  scores a fragment by weighting L1 profile distances by residue type
/// @author Robert Vernon


// type headers
#include <core/types.hh>

#include <protocols/frag_picker/scores/ProfileScoreDistWeight.hh>

// package headers
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

// project headers
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>

// mini headers
#include <core/sequence/SequenceProfile.hh>

// option key includes
#include <basic/options/keys/OptionKeys.hh>

// utils
#include <basic/prof.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer trProfScoreDistWeight(
	"protocols.frag_picker.scores.ProfileScoreDistWeight");

ProfileScoreDistWeight::ProfileScoreDistWeight(
	core::Size priority,
	core::Real lowest_acceptable_value,
	bool use_lowest,
	core::sequence::SequenceProfileOP query_profile,
	core::fragment::SecondaryStructureOP query_ss_prediction,
	std::string query_sequence, core::Size longest_vall_chunk
) :
	CachingScoringMethod(
	priority, lowest_acceptable_value, use_lowest, "ProfileScoreDistWeight" ),
	query_sequence_(query_sequence),
	query_profile_(query_profile),
	query_ss_(query_ss_prediction)
{
	for ( core::Size i = 1; i <= query_profile->length(); ++i ) {
		utility::vector1<core::Real> row(longest_vall_chunk);
		scores_.push_back(row);
	}


	utility::vector1< utility::vector1< utility::vector1 <core::Real> > > temp(
		3, utility::vector1< utility::vector1 <core::Real> > (
		20, utility::vector1<core::Real> (
		20, 0.0
		)
		)
	);
	distance_weights_ = temp;

	utility::io::izstream data("distances.txt");

	if ( !data ) {
		utility_exit_with_message("[ERROR] Unable to open distances.txt");
	}

	std::string line;
	char res1;
	char res2;
	char ss1;
	char ss2;
	core::Real dist;

	std::map<char,core::Size> ss_type_temp;
	ss_type_temp.insert(std::make_pair('H',1));
	ss_type_temp.insert(std::make_pair('E',2));
	ss_type_temp.insert(std::make_pair('L',3));

	ss_type_map_ = ss_type_temp;

	std::map<char,core::Size> aa_order_tmp;
	aa_order_tmp.insert(std::make_pair('A',1));
	aa_order_tmp.insert(std::make_pair('C',2));
	aa_order_tmp.insert(std::make_pair('D',3));
	aa_order_tmp.insert(std::make_pair('E',4));
	aa_order_tmp.insert(std::make_pair('F',5));
	aa_order_tmp.insert(std::make_pair('G',6));
	aa_order_tmp.insert(std::make_pair('H',7));
	aa_order_tmp.insert(std::make_pair('I',8));
	aa_order_tmp.insert(std::make_pair('K',9));
	aa_order_tmp.insert(std::make_pair('L',10));
	aa_order_tmp.insert(std::make_pair('M',11));
	aa_order_tmp.insert(std::make_pair('N',12));
	aa_order_tmp.insert(std::make_pair('P',13));
	aa_order_tmp.insert(std::make_pair('Q',14));
	aa_order_tmp.insert(std::make_pair('R',15));
	aa_order_tmp.insert(std::make_pair('S',16));
	aa_order_tmp.insert(std::make_pair('T',17));
	aa_order_tmp.insert(std::make_pair('V',18));
	aa_order_tmp.insert(std::make_pair('W',19));
	aa_order_tmp.insert(std::make_pair('Y',20));

	aa_order_map_ = aa_order_tmp;

	while ( getline(data, line) ) {
		std::istringstream line_stream(line);
		line_stream >> res1 >> ss1 >> res2 >> ss2 >> dist;


		core::Size res_type( 0 );
		if ( ss1 == ss2 ) {
			res_type = ss_type_map_.find(ss1)->second;
		}

		//std::cout << "AAAA: " << res1 << " " << res2 << " " << ss1 << " " << ss2 << " " << dist << " " << res_type << std::endl;

		if ( res_type != 0 ) {


			core::Size i_res1, i_res2;

			i_res1 = aa_order_map_.find(res1)->second;
			i_res2 = aa_order_map_.find(res2)->second;

			distance_weights_[res_type][i_res1][i_res2] = dist;
			//std::cout << "BBBB: " << res1 << " " << res2 << " " << ss1 << " " << ss2 << " " << res_type << " " << i_res1 << " " << i_res2 << " " << dist << std::endl;
		}
	}
}


void ProfileScoreDistWeight::do_caching(VallChunkOP chunk) {

	std::string tmp = chunk->chunk_key();
	if ( tmp == cached_scores_id_ ) {
		return;
	}
	cached_scores_id_ = tmp;
	//core::Size size_q = query_profile_->length();

	trProfScoreDistWeight.Debug << "caching profile score for " << chunk->get_pdb_id()
		<< " of size " << chunk->size() << std::endl;
	PROF_START( basic::FRAGMENTPICKING_PROFILE_CAHING );
	//for (core::Size i = 1; i <= size_q; ++i) {
	//std::cout << "A " << query_sequence_ << " " << query_sequence_.length() << std::endl;
	for ( core::Size i = 0; i < query_sequence_.length(); ++i ) {
		//std::cout << "B" << chunk->size() << std::endl;
		//utility::vector1<core::Real> query_prof_row = query_profile_->prof_row(i);
		core::Size seqpos_res_id(aa_order_map_.find(query_sequence_[i])->second);


		for ( core::Size j = 1; j <= chunk->size(); ++j ) {
			//std::cout << "C" << chunk->size() << std::endl;

			//utility::vector1<core::Real> tmplt_prof_row = chunk->at(j)->profile();
			core::Real score(0.0);

			core::Size tmplt_res_id (aa_order_map_.find(chunk->at(j)->aa())->second);

			for ( core::Size s = 1; s <= 3; s++ ) {
				core::Real ss_weight(0.0);
				if ( s == 1 ) {
					ss_weight = query_ss_->helix_fraction(i+1);
				}
				if ( s == 2 ) {
					ss_weight = query_ss_->strand_fraction(i+1);
				}
				if ( s == 3 ) {
					ss_weight = query_ss_->loop_fraction(i+1);
				}

				//for (core::Size v = 1; v <= 20; v++) {
				// core::Real distance_weight(0.0);
				//
				// distance_weight = distance_weights_[s][seqpos_res_id][v];
				// score += tmplt_prof_row[v]*ss_weight*distance_weight;
				//}

				core::Real distance_weight(0.0);

				distance_weight = distance_weights_[s][seqpos_res_id][tmplt_res_id];
				score += ss_weight*distance_weight;

				//score += std::abs(tmplt_prof_row[v] - query_prof_row[q])
				// *ss_weight
				// *distance_weight;
				// }

				//std::cout << "SCOREIJ_FINAL " << score << std::endl;
			}

			//std::cout << "SCOREIJ_FINAL " << score << " " <<  i+1 << " " << j << " " << query_sequence_[i] << " " << chunk->at(j)->aa() << " " << query_ss_->helix_fraction(i+1) << " " << query_ss_->strand_fraction(i+1) << " " << query_ss_->loop_fraction(i+1) << " " << distance_weights_[1][seqpos_res_id][tmplt_res_id] << " " << distance_weights_[2][seqpos_res_id][tmplt_res_id] << " " << distance_weights_[3][seqpos_res_id][tmplt_res_id] << std::endl;

			scores_[i+1][j] = score;
		}
	}

	PROF_STOP( basic::FRAGMENTPICKING_PROFILE_CAHING );
	trProfScoreDistWeight.Debug << "precomputed matrix of scores " << scores_.size()
		<< "x" << chunk->size() << std::endl;
}

bool ProfileScoreDistWeight::cached_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map) {

	std::string tmp = f->get_chunk()->chunk_key();
	if ( tmp != cached_scores_id_ ) {
		do_caching(f->get_chunk());
	}

	core::Real totalScore = 0.0;//f->get_length() * 20.0;
	for ( core::Size i = 1; i <= f->get_length(); i++ ) {
		debug_assert(f->get_first_index_in_query() + i - 1 <= scores_.size());
		debug_assert(f->get_first_index_in_vall()
			+ i - 1<= scores_[1].size());
		totalScore += scores_[f->get_first_index_in_query() + i - 1][f->get_first_index_in_vall() + i - 1];
		//std::cout << "TOTALSCORE " << totalScore << " " << scores_[f->get_first_index_in_query() + i - 1][f->get_first_index_in_vall() + i - 1] << " " <<  f->get_first_index_in_query() + i - 1 << " " << f->get_first_index_in_vall() + i - 1 << std::endl;
	}
	totalScore *= (core::Real) 100.0;
	totalScore /= (core::Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore < lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

FragmentScoringMethodOP MakeProfileScoreDistWeight::make(core::Size priority,
	core::Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker, std::string prediction_id) {

	//query_sequence_ = picker->get_query_seq_string();

	//std::cout << "QUERY_SEQUENCE " << query_sequence_ << std::endl;

	core::Size len = picker->get_vall()->get_largest_chunk_size();

	//std::istringstream line_stream(config_line);
	//std::string score_name;
	//core::Size p;
	//core::Real weight;
	//core::Real lowest;
	//std::string prediction_id;
	//line_stream >> score_name >> p >> weight >> lowest >> prediction_id;

	trProfScoreDistWeight << "Profile scoring method is: DistWeight" << std::endl;
	core::fragment::SecondaryStructureOP query_prediction( picker->get_query_ss(prediction_id) );
	if ( ! query_prediction ) {
		utility_exit_with_message("Unable to find secondary structure prediction for " + prediction_id );
	}
	return (FragmentScoringMethodOP) utility::pointer::make_shared< ProfileScoreDistWeight >(priority,
		lowest_acceptable_value, use_lowest, picker->get_query_seq(), query_prediction, picker->get_query_seq_string(),len);
}

} //scores
} // frag_picker
} // protocols
