// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/FragmentCrmsd.cc
/// @brief  Scores fragments by disulfide-linke Calpha distances
/// @author Robert Vernon

#include <protocols/frag_picker/scores/DisulfideDistance.hh>

#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/FragmentPicker.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/io/raw_data/DisulfideFile.hh>

// utils
#include <ObjexxFCL/FArray1D.hh>
#include <basic/prof.hh>


#include <utility>
#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

DisulfideDistance::DisulfideDistance(core::Size priority, core::Real lowest_acceptable_value, bool use_lowest,
	utility::vector1< core::Size > const & disulfide_data, core::Size largest_fragment) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "DisulfideDistance"),
	disulfide_data_(disulfide_data),
	largest_fragment_(largest_fragment)
{
}

void DisulfideDistance::do_caching(VallChunkOP current_chunk) {

	n_res_ = current_chunk->size();

	chunk_ca_distances_.redimension(n_res_, n_res_, 0.0);
	for ( core::Size x = 1; x <= n_res_; ++x ) {
		VallResidueOP xr = current_chunk->at(x);

		for ( core::Size y = x - largest_fragment_ - 3; (y <= x + largest_fragment_ + 3) && (y <= n_res_); ++y ) {
			//for (core::Size y = 1; y <= n_res_; ++y) {
			if ( y > 0 ) {
				VallResidueOP yr = current_chunk->at(y);

				core::Real distance = sqrt( pow( xr->x() - yr->x(), 2) + pow( xr->y() - yr->y(), 2) + pow( xr->z() - yr->z(), 2) );
				chunk_ca_distances_(x,y) = distance;
			}
		}
	}
}

bool DisulfideDistance::score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {
	return cached_score( fragment, scores );
}

bool DisulfideDistance::cached_score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {

	std::string tmp = fragment->get_chunk()->chunk_key();
	if ( tmp.compare(cached_scores_id_) != 0 ) {
		do_caching(fragment->get_chunk());
		cached_scores_id_ = tmp;
	}

	core::Size offset_q = fragment->get_first_index_in_query() - 1;
	core::Size offset_v = fragment->get_first_index_in_vall() - 1;
	core::Real score = 0.0;

	//std::cout << "CACHING_FOR_DISULFIDE_DISTANCES" << std::endl;

	for ( core::Size i = 1; i < fragment->get_length(); ++i ) {

		if ( (i+offset_q) <= disulfide_data_.size() ) {
			if ( disulfide_data_[i+offset_q] != 0 ) {

				core::Size res1 = i+offset_q;//disulfide_[i+offset_q].first;
				core::Size res2 = disulfide_data_[i+offset_q];//.second;

				core::Size v1 = i+offset_v;

				core::Size v2(0);
				core::Size seq_sep(0);

				if ( res2 > res1 ) {
					seq_sep = res2 - res1;
					v2 = v1 + seq_sep;
				} else {
					seq_sep = res1 - res2;
					if ( v1 > (res1 - res2) ) {
						v2 = v1 - seq_sep;
					}
				}

				//std::cout << "CACHED " << res1 << " " << res2 << " " << n_res_ << std::endl;

				if ( (res2 != 0) && (res2 <= (3 + offset_q + fragment->get_length())) && (res2 >= offset_q - 3)
						&& (v2 > 0) && (v2 <= n_res_) ) {

					if ( seq_sep >= 4 ) {
						if ( (chunk_ca_distances_(v1,v2) < 3.6) || (chunk_ca_distances_(v1,v2) > 7.0) ) {
							score += 10.0;
						}
					} else {
						if ( seq_sep >= 2 ) {
							if ( (chunk_ca_distances_(v1,v2) < 4.6) || (chunk_ca_distances_(v1,v2) > 6.7) ) {
								score += 10.0;
							}
						} else {
							if ( (chunk_ca_distances_(v1,v2) > 4.0) ) {
								score += 10.0;
							}
						}
					}
				}

				//std::cout << "DISULFIDE_SCORING " << res1 << " " << res2 << " " << v1 << " " << v2 << " " << offset_q << " " << offset_v << " " << i << " " << score << " " << chunk_ca_distances_(v1,v2) << std::endl;
			}
		}
	}

	scores->set_score_component(score, id_);
	PROF_STOP( basic::FRAGMENTPICKING_PHIPSI_SCORE );
	if ( (score > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}

	return true;
}

void DisulfideDistance::clean_up() {
}

/// @brief Creates a DisulfideDistance scoring method
/// @param priority - priority of the scoring method. The higher value the earlier the score
///  will be evaluated
/// @param lowest_acceptable_value - if a calculated score is higher than this value,
///  fragment will be neglected
/// @param FragmentPickerOP object - not used
/// @param line - the relevant line extracted from the scoring configuration file that defines this scoring method
///   It could look like: "DisulfideDistance                140     -5.0     100.0 additional_string"
///  where 140, -5.0 && 100.0 are priority, weight && treshold, respectively.
///  The additional string may be:
///  - empty: then the maker tries to create a scoring object from a TALOS file
///   trying in::file::talos_phi_psi flag. If fails, will try to use a pose from in::file::s
///  - a pdb file, pdb extension is necessary. This will create a pose && steal Phi && Psi
///  - a TALOS file with Phi/Psi prediction (tab extension is necessary)
FragmentScoringMethodOP MakeDisulfideDistance::make(core::Size priority,
	core::Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker//picker
	, std::string )
{

	if ( option[in::fix_disulf].user() ) {

		core::io::raw_data::DisulfideFile ds_file( option[ in::fix_disulf ]() );

		utility::vector1< std::pair<core::Size,core::Size> > disulfides_in_file;

		ds_file.disulfides(disulfides_in_file);

		core::Size largest_number(0);

		for ( core::Size i = 1; i <= disulfides_in_file.size(); ++i ) {

			core::Size l = disulfides_in_file[i].first;
			core::Size u = disulfides_in_file[i].second;

			if ( u <= l ) {
				utility_exit_with_message("[ERROR] Disulfide File Format: res2 must be > res1");
			}

			if ( u > largest_number ) {
				largest_number = u;
			}
		}

		utility::vector1< core::Size > disulfide_data(largest_number,0);

		for ( core::Size i = 1; i <= disulfides_in_file.size(); ++i ) {

			auto l = static_cast< core::Size > ( disulfides_in_file[i].first );
			auto u = static_cast< core::Size > ( disulfides_in_file[i].second );

			disulfide_data[u] = l;
			disulfide_data[l] = u;
		}

		for ( core::Size i = 1; i <= largest_number; ++i ) {

			std::cout << "DISULFIDE_DATA " << i << " " << disulfide_data[i] << std::endl;

		}

		core::Size largest_fragment(0);
		for ( core::Size i = 1; i <= picker->frag_sizes_.size(); ++i ) {
			if ( picker->frag_sizes_[i] > largest_fragment ) largest_fragment = picker->frag_sizes_[i];
		}

		runtime_assert( largest_fragment > 0 );

		//disulfide_data_ = disulfide_data;
		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new DisulfideDistance(priority, lowest_acceptable_value, use_lowest,
			disulfide_data, largest_fragment) );
	}

	utility_exit_with_message(
		"Can't read disulfide connectivity file. Provide a connectivity file with -in::fix_disulf <file>\n");

	return nullptr;

}

} // scores
} // frag_picker
} // protocols

