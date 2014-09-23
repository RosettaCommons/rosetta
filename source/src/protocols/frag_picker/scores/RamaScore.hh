// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/RamaScore.hh
/// @brief  a base class for fragment scoring
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_scores_RamaScore_hh
#define INCLUDED_protocols_frag_picker_scores_RamaScore_hh

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/FragmentPicker.hh>

#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

#include <core/fragment/SecondaryStructure.fwd.hh>

namespace protocols {
namespace frag_picker {
namespace scores {


class RamaScore;
typedef utility::pointer::shared_ptr< RamaScore > RamaScoreOP;
typedef utility::pointer::shared_ptr< RamaScore const > RamaScoreCOP;

/// @brief  RamaScore score counts identical residues
/// @detailed Resulting score is the number of identical residues
/// on corresponding positions in a vall fragment and a query sequence
class RamaScore: public CachingScoringMethod {
public:

	RamaScore(Size, Real, bool, std::string&,std::string);

	RamaScore(Size, Real, bool, std::string&, core::fragment::SecondaryStructureOP,std::string);

	void SetupRamaTables();

	void do_caching(VallChunkOP);
	void clean_up();

	bool score(FragmentCandidateOP f, FragmentScoreMapOP empty_map);
	bool cached_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map);

	/// @brief prints a detailed explanation how a fragment score has been computed
	/// @detailed besides extensive output, the method should return the same result as score()
	bool describe_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map,
			std::ostream& out);

	std::string & get_prediction_name() { return prediction_name_; }

private:
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// Real minScoreAllowed_;
	std::string& query_;
	core::fragment::SecondaryStructureOP query_ss_;
	std::string prediction_name_;

	std::string cached_scores_id_;

	utility::vector1< utility::vector1< Real > > scores_;

	//Sequence Specific Ramachandran Surfaces (keep only one copy of this so we don't have to keep reading it from the database for each instance)
	static utility::vector1< utility::vector1< utility::vector1< Real > > > sequence_rama_tables_;

};

/// @brief  Maker class that produces a new RamaScore object
class MakeRamaScore: public MakeFragmentScoringMethod {
public:

	MakeRamaScore() :
		MakeFragmentScoringMethod("RamaScore") {
	}

	FragmentScoringMethodOP make(Size priority, Real lowest_acceptable_value, bool use_lowest,
			FragmentPickerOP picker, std::string prediction_id ) {

		if ( prediction_id == "" ) {
			core::fragment::SecondaryStructureOP default_ss( new core::fragment::SecondaryStructure );
			Size query_len = picker->get_query_seq_string().size();
			default_ss->extend( query_len );
			for ( Size i = 1; i <= query_len; ++i ) {
				default_ss->set_fractions(i, 0.333, 0.333, 0.334);
			}

			//			std::cout << "CO_PREDICTION_ID " << prediction_id << std::endl;
			std::string default_prediction_name("uniform_default");
			return (FragmentScoringMethodOP) FragmentScoringMethodOP( new RamaScore(priority,
				lowest_acceptable_value, use_lowest, picker->get_query_seq_string(), default_ss,default_prediction_name) );
		} else {

			//std::cout << "PREDICTION_ID " << prediction_id << std::endl;
			core::fragment::SecondaryStructureOP query_prediction( picker->get_query_ss(prediction_id) );
			if( ! query_prediction ) {
				utility_exit_with_message( "Unable to find secondary structure prediction for " + prediction_id );
			}
			return (FragmentScoringMethodOP) FragmentScoringMethodOP( new RamaScore(priority,
				lowest_acceptable_value, use_lowest, picker->get_query_seq_string(), query_prediction,prediction_id) );
		}
	}
};

} // scores
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_scores_RamaScore_HH */
