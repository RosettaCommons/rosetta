// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/FragmentCrmsd.cc
/// @brief  Object that scores a fragment by root mean square deviation of Phi && Psi dihedrals
/// @author Nikolas Sgourakis sgourn@u.w.edu


#include <protocols/frag_picker/scores/JCoupling.hh>

#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#ifdef WIN32
#include <protocols/frag_picker/FragmentPicker.hh>
#endif

#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <numeric/conversions.hh>

// option key includes
#include <basic/options/keys/OptionKeys.hh>

#include <basic/Tracer.hh>

// utils

#include <utility/vector1.hh>


// C++

// Boost

namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer trJCoupling(
	"protocols.frag_picker.scores.JCoupling");

JCoupling::JCoupling(core::Size priority, core::Real lowest_acceptable_value, bool use_lowest,
	JCouplingIO& reader) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "JCoupling"), data_(reader) {

	len_ = data_.get_length();
	A_ = data_.get_parameters()[1];
	B_ = data_.get_parameters()[2];
	C_ = data_.get_parameters()[3];
	THETA_ = data_.get_parameters()[4];
}

void JCoupling::do_caching(VallChunkOP current_chunk) {
	std::string ctmp = current_chunk->chunk_key();
	if ( ctmp.compare("change to 'cached_scores_id_' when ready") != 0 ) {
		return; // CACHING NOT BUILT IN YET
	}
}

bool JCoupling::cached_score(FragmentCandidateOP fragment,
	FragmentScoreMapOP scores) {


	return score( fragment, scores );

	//  std::string ctmp = fragment->get_chunk()->chunk_key();
	//  if (ctmp.compare(cached_scores_id_) != 0) {
	//   do_caching(fragment->get_chunk());
	//   cached_scores_id_ = ctmp;
	//  }

	//  PROF_START( basic::FRAGMENTPICKING_PHIPSI_SCORE );
	//  core::Size offset_q = fragment->get_first_index_in_query() - 1;
	//  core::Size offset_v = fragment->get_first_index_in_vall() - 1;
	//  core::Real score = 0.0;
	//  core::Real tmp = 0.0;
	//  for (core::Size i = 1; i < fragment->get_length(); ++i) {
	//   if (!existing_data_[i + offset_q])
	//    continue;
	//   core::Real d = 0.0;

	//   tmp = std::abs(chunk_phi_(i + offset_v) - query_phi_(i + offset_q));

	//   if ( tmp > 180.0 ) tmp = std::abs(tmp - 360.0);

	//   if ( tmp > query_d_phi_(i + offset_q) ) {
	//    tmp = tmp - query_d_phi_(i + offset_q);
	//   } else {
	//    tmp = 0.0;
	//   }

	//   d += tmp * tmp;

	//   tmp = std::abs(chunk_psi_(i + offset_v) - query_psi_(i + offset_q));

	//   if ( tmp > 180.0 ) tmp = std::abs(tmp - 360.0);

	//   if ( tmp > query_d_psi_(i + offset_q) ) {
	//    tmp = tmp - query_d_psi_(i + offset_q);
	//   } else {
	//    tmp = 0.0;
	//   }

	//   d += tmp * tmp;

	//   score += std::sqrt(d);
	//  }

	//  score = score / ((core::Real) fragment->get_length());
	//  PROF_STOP( basic::FRAGMENTPICKING_PHIPSI_SCORE );

	//  scores->set_score_component(score, id_);
	//  if ((score > lowest_acceptable_value_) && (use_lowest_ == true))
	//   return false;

	//  return true;
}


bool JCoupling::score(FragmentCandidateOP fragment,
	FragmentScoreMapOP scores) {

	core::Size offset_q = fragment->get_first_index_in_query();
	core::Size offset_v = fragment->get_first_index_in_vall();

	VallChunkOP chunk = fragment->get_chunk();
	core::Real score = 0.0;

	for ( core::Size i = 1; i <= fragment->get_length(); ++i ) {

		core::Size query_res = offset_q + i - 1;
		core::Size vall_res = offset_v + i - 1;

		VallResidueOP r = chunk->at( vall_res );
		core::Real cos_phi( cos( numeric::conversions::radians(r->phi()) + THETA_ ) );

		bool has_data(false);

		std::pair< core::Real, core::Real > datum( data_.get_data( query_res, has_data ) );

		if ( has_data ) {
			core::Real val = datum.first;
			core::Real dev = datum.second;

			core::Real tmp((val -  (A_ * cos_phi * cos_phi + B_ * cos_phi + C_)) / dev);

			//   std::cout << "COMPUTE " << query_res << " " << val << " " << dev << " " << A_ << " " << B_ << " " << C_ << " " << THETA_ << " " << r->phi() << " " << r->psi() << " " << tmp << " " << (numeric::conversions::radians(r->phi()) + THETA_) << " " << cos( numeric::conversions::radians(r->phi()) + THETA_ ) << std::endl;

			score += tmp*tmp;

		}// else {
		// std::cout << "NO DATA " << query_res << std::endl;
		//}

	}

	score = score / ((core::Real) fragment->get_length());

	scores->set_score_component(score, id_);
	if ( (score > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}

	return true;
}

void JCoupling::clean_up() {
}

/// @brief Creates a JCoupling scoring method
/// @param priority - priority of the scoring method. The higher value the earlier the score
///  will be evaluated
/// @param lowest_acceptable_value - if a calculated score is higher than this value,
///  fragment will be neglected
/// @param FragmentPickerOP object - not used
/// @param extras - additional parameters to create a new object. Allowed values are:
///  - empty: then the maker tries to create a scoring object from a TALOS file
///   trying in::file::talos_phi_psi flag. If fails, will try to use a pose from in::file::s
///  - a pdb file, pdb extension is necessary. This will create a pose && steal Phi && Psi
///  - a TALOS file with Phi/Psi prediction (tab extension is necessary)
FragmentScoringMethodOP MakeJCoupling::make(core::Size priority,
	core::Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP //picker
	, std::string input_file) {

	if ( input_file != "" ) {

		trJCoupling
			<< "Experimental Data for JCoupling scoring loaded from file: "
			<< input_file << std::endl;
		JCouplingIO reader(input_file);
		//in.write(trJCoupling.Debug);
		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new JCoupling(priority,
			lowest_acceptable_value, use_lowest, reader) );
	}

	utility_exit_with_message(
		"Can't read JCoupling file.");

	return NULL;
}

} // scores
} // frag_picker
} // protocols


