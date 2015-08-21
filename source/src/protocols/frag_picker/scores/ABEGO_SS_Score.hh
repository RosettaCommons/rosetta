// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/ABEGO_SS_Score.hh
/// @brief  scores a fragment by secondary structure similarity
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_ABEGO_SS_Score_hh
#define INCLUDED_protocols_frag_picker_scores_ABEGO_SS_Score_hh

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>

#include <protocols/frag_picker/quota/ABEGO_SS_Config.hh>
#include <protocols/frag_picker/quota/ABEGO_SS_Map.hh>


#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

typedef utility::vector1<utility::vector1<Real> > Matrix;

class ABEGO_SS_Score;
typedef utility::pointer::shared_ptr< ABEGO_SS_Score > ABEGO_SS_ScoreOP;
typedef utility::pointer::shared_ptr< ABEGO_SS_Score const > ABEGO_SS_ScoreCOP;

class ABEGO_SS_Score: public CachingScoringMethod {
public:

	ABEGO_SS_Score(Size priority, Real lowest_acceptable_value, bool use_lowest,
		std::string prediction_file_name,Size longest_vall_chunk) :
		CachingScoringMethod(priority, lowest_acceptable_value, use_lowest,
		"ABEGO_SS_Score") {

		quota::ABEGO_SS_Config prediction_file(prediction_file_name);
		query_len_ = prediction_file.size();
		n_classes_ = prediction_file.n_columns();
		for ( Size i = 1; i <= query_len_; ++i ) {
			utility::vector1<Real> row(longest_vall_chunk);
			scores_.push_back(row);
		}
		for ( Size iseq=1; iseq<=query_len_; iseq++ ) {
			utility::vector1<Real> row;
			for ( Size ibin=1; ibin<=n_classes_; ibin++ ) {
				row.push_back(prediction_file.probability(iseq,ibin));
			}
			ratios_.push_back( row );
		}
		for ( Size ibin=1; ibin<=n_classes_; ibin++ ) {
			maps_.push_back( utility::pointer::shared_ptr<class protocols::frag_picker::quota::ABEGO_SS_Map>( new quota::ABEGO_SS_Map(prediction_file.get_pool_bins(ibin)) ) );
		}
	}

	~ABEGO_SS_Score() {}

	void do_caching(VallChunkOP);
	bool cached_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map);
	void clean_up() {}

	/// @brief Computes the score
	virtual bool score(FragmentCandidateOP, FragmentScoreMapOP);

private:
	Size query_len_;
	Size n_classes_;
	utility::vector1<quota::ABEGO_SS_MapOP> maps_;
	utility::vector1< utility::vector1<Real> > ratios_;
	std::string cached_scores_id_;
	Matrix scores_;
};

/// @brief  Maker class that produces a new ABEGO_SS_Score object
class MakeABEGO_SS_Score: public MakeFragmentScoringMethod {
public:

	MakeABEGO_SS_Score() :
		MakeFragmentScoringMethod("ABEGO_SS_Score") {
	}

	FragmentScoringMethodOP make(Size, Real, bool,
		FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_ABEGO_SS_Score_HH */
