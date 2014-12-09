// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/FragmentCrmsd.cc
/// @brief  Scores fragments by disulfide-linke Calpha distances
/// @author Robert Vernon

#include <protocols/frag_picker/scores/PCS_FragDistance.hh>

#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/FragmentPicker.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// AUTO-REMOVED #include <core/io/raw_data/DisulfideFile.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
//#include <core/util/Tracer.hh>

#include <utility/io/izstream.hh>

// utils
#include <ObjexxFCL/FArray1D.hh>
//#include <core/util/prof.hh>

// AUTO-REMOVED #include <boost/tuple/tuple.hpp>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

	PCS_FragDistance::PCS_FragDistance(Size priority, Real lowest_acceptable_value, bool use_lowest,
																		 ObjexxFCL::FArray2D_double target_ca_dev,
																		 ObjexxFCL::FArray2D_double target_ca_dist, Size largest_fragment,
																		 Size max_length) :
		CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "PCS_FragDistance"),
		target_ca_dev_(target_ca_dev),
		target_ca_dist_(target_ca_dist),
	  largest_fragment_(largest_fragment),
		max_length_(max_length)
{
}

void PCS_FragDistance::do_caching(VallChunkOP current_chunk) {

	n_res_ = current_chunk->size();

	chunk_ca_distances_.redimension(n_res_, n_res_, 0.0);
	for (Size x = 1; x <= n_res_; ++x) {
		VallResidueOP xr = current_chunk->at(x);

		for (Size y = x - largest_fragment_; (y <= x + largest_fragment_) && (y <= n_res_); ++y) {
			//for (Size y = 1; y <= n_res_; ++y) {
			if (y > 0) {
				VallResidueOP yr = current_chunk->at(y);

				Real distance = sqrt( pow( xr->x() - yr->x(), 2) + pow( xr->y() - yr->y(), 2) + pow( xr->z() - yr->z(), 2) );
				chunk_ca_distances_(x,y) = distance;
			}
		}
	}
}

bool PCS_FragDistance::score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {
	return cached_score( fragment, scores );
}

bool PCS_FragDistance::cached_score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {

	std::string tmp = fragment->get_chunk()->chunk_key();
	if (tmp.compare(cached_scores_id_) != 0) {
		do_caching(fragment->get_chunk());
		cached_scores_id_ = tmp;
	}

	//Size offset_q = fragment->get_first_index_in_query() - 1;
	//Size offset_v = fragment->get_first_index_in_vall() - 1;

	//Goes from first res - 2 to last res + 2
	Size offset_q = fragment->get_first_index_in_query() - 1;
	Size offset_v = fragment->get_first_index_in_vall() - 1;
	Real score = 0.0;

	//for (Size ix = 1; ix < fragment->get_length(); ++ix) {
	//	for (Size iy = ix+1; iy < fragment->get_length(); ++iy) {

	for (Size ix = 1; ix <= fragment->get_length()+3; ++ix) {
		for (Size iy = ix+1; iy <= fragment->get_length()+4; ++iy) {

			Size res1 = ix+offset_q;
			Size res2 = iy+offset_q;

			Size v1 = ix+offset_v;
			Size v2 = iy+offset_v;

			if ( ( res1 >= 3) && (res2 <= max_length_ + 2) &&
					 ( v1 >= 3) && (v2 <= n_res_ + 2) ) {

				res1 -= 2;
				res2 -= 2;
				v1 -= 2;
				v2 -= 2;

				//std::cout << "HEYO " << res1 << " " << res2 << " " << v1 << " " << v2 << " " << max_length_ << " " << n_res_ << std::endl;

				if ( !((target_ca_dist_(res1,res2) == 0.0) && (target_ca_dev_(res1,res2) == 0.0)) ) {

					Real dev(target_ca_dev_(res1,res2));
					Real diff( std::abs(target_ca_dist_(res1,res2) - chunk_ca_distances_(v1,v2)));

					Real sig_function = ( 1 / ( 1 + exp(-2*(diff/dev) + 5 )));// + (1 / ( 1 + exp(+2*(diff/dev) + 5 )));

					//std::cout << "HEYO " << res1 << " " << res2 << " " << target_ca_dev_(res1,res2) << " " << target_ca_dev_(res1,res2) << " " << chunk_ca_distances_(v1,v2) << " " << diff << " " << dev << " " << sig_function << std::endl;

					score += sig_function;
				}
			}
		}
	}

	score /= (Real) fragment->get_length();

	scores->set_score_component(score, id_);
	//PROF_STOP( core::util::FRAGMENTPICKING_PHIPSI_SCORE );
	if ((score > lowest_acceptable_value_) && (use_lowest_ == true))
		return false;

	return true;
}

void PCS_FragDistance::clean_up() {
}

/// @brief Creates a PCS_FragDistance scoring method
/// @param priority - priority of the scoring method. The higher value the earlier the score
///		will be evaluated
/// @param lowest_acceptable_value - if a calculated score is higher than this value,
///		fragment will be neglected
/// @param FragmentPickerOP object - not used
/// @param line - the relevant line extracted from the scoring configuration file that defines this scoring method
/// 		It could look like: "PCS_FragDistance                140     -5.0     100.0 additional_string"
///		where 140, -5.0 && 100.0 are priority, weight && treshold, respectively.
///		The additional string may be:
///		- empty: then the maker tries to create a scoring object from a TALOS file
///			trying in::file::talos_phi_psi flag. If fails, will try to use a pose from in::file::s
///		- a pdb file, pdb extension is necessary. This will create a pose && steal Phi && Psi
///		- a TALOS file with Phi/Psi prediction (tab extension is necessary)
FragmentScoringMethodOP MakePCS_FragDistance::make(Size priority,
		Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP picker//picker
		, std::string )
{

	std::string target( picker->get_query_seq_string());

	Size max_length = target.length();

	ObjexxFCL::FArray2D_double target_ca_dev;
	ObjexxFCL::FArray2D_double target_ca_dist;

	target_ca_dev.redimension(max_length, max_length, 0.0);
	target_ca_dist.redimension(max_length, max_length, 0.0);


	if (option[in::file::PCS_frag_cst].user()) {

		utility::io::izstream data( option[ in::file::PCS_frag_cst ]() );

		// mjo commenting out 'max_data_res' because it is not used and casuses a warning
		//Size max_data_res(0);
		std::string line;
		while (!data.eof()) {
			getline(data, line);
			std::istringstream line_stream(line);

			Size resX, resY;
			Real dist, dev;

			line_stream >> resX >> resY >> dist >> dev;

			//if ( dev < 0.3 ) {
			//	dev = 0.3;
				//This was the lowest dev I observed during training.
				//The score has not been tested for lower values and
				//therefore you should not expect it to work for lower values!
			//}

			target_ca_dist(resX,resY) = dist;
			target_ca_dev(resX,resY) = dev;
		}

		Size largest_fragment(0);
		for ( Size i = 1; i <= picker->frag_sizes_.size(); ++i ) {
			if ( picker->frag_sizes_[i] > largest_fragment ) largest_fragment = picker->frag_sizes_[i];
		}

		runtime_assert( largest_fragment > 0 );

		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new PCS_FragDistance(priority, lowest_acceptable_value, use_lowest,
																													target_ca_dev, target_ca_dist, largest_fragment, max_length) );
	}

		utility_exit_with_message(
			"Can't read PCS_frag_cst file. Provide a connectivity file with -in::PCS_frag_cst <file>\n");

	return NULL;

}

} // scores
} // frag_picker
} // protocols

