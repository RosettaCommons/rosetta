// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/TorsionBin.hh
/// @brief  scores a fragment based on scores for a probability over torsion bins.
/// @author James Thompson

#ifndef INCLUDED_protocols_frag_picker_scores_LAMBEGO_Similarity_hh
#define INCLUDED_protocols_frag_picker_scores_LAMBEGO_Similarity_hh

// type headers
#include <core/types.hh>

#include <protocols/frag_picker/VallChunk.fwd.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <utility/vector1_bool.hh>
#include <iostream>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

typedef utility::vector1< utility::vector1 < core::Real > > Matrix;

/// @brief scores a fragment by torsion bin similarity
class LAMBEGO_Similarity: public CachingScoringMethod {
public:

	LAMBEGO_Similarity(
		core::Size priority,
		core::Real lowest_acceptable_value,
		bool use_lowest,
		utility::vector1< utility::vector1< core::Real > > query_bin_probs,
		core::Size sequence_length,
		core::Size longest_vall_chunk
	) :
		CachingScoringMethod(
		priority, lowest_acceptable_value, use_lowest,
		"LAMBEGO"
		),
		query_len_( sequence_length ),
		query_bin_probs_( query_bin_probs )
	{
		for ( core::Size i = 1; i <= query_len_; ++i ) {
			utility::vector1< core::Real > row(longest_vall_chunk);
			scores_.push_back(row);
		}
	}

	void do_caching(VallChunkOP);

	void clean_up() {};

	/// @brief Computes the score
	virtual bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	char torsion2big_bin_(
		core::Real const phi,
		core::Real const psi,
		core::Real const omega,
		char const ss
	) const;

	core::Size bin_index_( char const bin_name ) const;

protected:
	Matrix scores_;

private:
	std::string name_;
	core::Size query_len_;
	utility::vector1< utility::vector1< core::Real > > query_bin_probs_;
	std::string cached_scores_id_;
}; // LAMBEGO_Similarity

/// @brief  Maker class that produces a new TorsionBin object
class MakeLAMBEGO_Similarity: public MakeFragmentScoringMethod {
public:

	MakeLAMBEGO_Similarity() :
		MakeFragmentScoringMethod("LAMBEGO_Similarity")
	{}

	FragmentScoringMethodOP make(
		core::Size priority,
		core::Real lowest_acceptable_value,
		bool use_lowest,
		FragmentPickerOP picker,
		std::string /* prediction_id */
	);
};

} // scores
} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_scores_TorsionBin_HH */
