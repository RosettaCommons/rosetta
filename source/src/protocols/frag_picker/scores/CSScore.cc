// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/CSScore.cc
/// @brief  Object that scores a fragment by target-observed/vall-predicted chemical shift distances
/// @author Robert Vernon rvernon@u.washington.edu

#include <protocols/frag_picker/scores/CSScore.hh>

#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/CSTalosIO.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <utility/io/ozstream.hh>


// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>


// Boost
#include <boost/tuple/tuple.hpp>

// project headers
#include <basic/Tracer.hh>

#include <protocols/frag_picker/CS2ndShift.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer trCSScore(
	"protocols.frag_picker.scores.CSScore");


//CSScore Constructor
// The Talos file reader is passed in as an object, and secondary shifts are calculated during CSScore construction
// (Secondary shifts are shift deviations from random coil, where random coil values are defined according to
// the combination of atom type, residue type, previous residue type, and next residue type.
CSScore::CSScore(core::Size priority, core::Real lowest_acceptable_value, bool use_lowest,
	CSTalosIO& reader) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "CSScore")
{

	//outfile_ = utility::io::ozstream tmp("allcomparisons.out");
	//outfile_.open("allcomparisons.out");

	trCSScore.Debug << "READING SHIFTS!" << std::endl;
	CS2ndShift secondary_shift_calculator(reader, true);
	trCSScore.Debug << "SHOULD BE DONE WRITING 2nd SHIFTS" << std::endl;

	target_shifts_ = secondary_shift_calculator.shifts();
}

//Residue-Residue scores are cached for current vall chunk
//This is where the CSScore equation lives
void CSScore::do_caching(VallChunkOP current_chunk) {

	////clip_factor is used to define the maximum shift difference
	////larger differences are adjusted down to clip_factor*v_shift
	////(the new score is sigmoidal, so clip factors are not required)
	////ONLY USED IN THE OLD MFR VERSION OF THE SCORE
	//core::Real const clip_factor(3.0);

	//bool vall_data(false);
	trCSScore.Debug << "caching CS score for " << current_chunk->get_pdb_id()
		<< " of size " << current_chunk->size() << std::endl;

	//Check to see if the cache needs to be recalculated
	std::string & tmp = current_chunk->chunk_key();
	if ( tmp.compare(cached_scores_id_) == 0 ) {
		return;
	}
	cached_scores_id_ = tmp;

	//Initialize empty 2D table, vall-length x target-length
	core::Size query_sequence_length = target_shifts_.size();
	std::pair< core::Real, core::Real > empty(0,0);
	utility::vector1< utility::vector1< std::pair< core::Real, core::Real> > > temp( current_chunk->size(),
		utility::vector1<std::pair< core::Real, core::Real> > (query_sequence_length, empty ) );
	//runtime_assert( target_shifts_.size() > 0 );

	//SIGMOID CONSTANTS - Should be set in constructor, not command line flags
	core::Real a( option[frags::sigmoid_cs_A]() ); // default = 4
	core::Real b( option[frags::sigmoid_cs_B]() ); // default = 5


	//Loop logic is "For each target x vall residue comparison, sum up total of
	//all shift differences"
	for ( core::Size r = 1; r <= target_shifts_.size(); ++r ) {
		utility::vector1< std::pair< core::Size, core::Real > > query_residue_shifts(target_shifts_[r]);
		for ( core::Size i = 1; i <= current_chunk->size(); ++i ) {
			core::Real tmp = 0.0;
			core::Real count = 0.0;
			for ( core::Size d = 1; d <= query_residue_shifts.size(); ++d ) {

				//q_shift_type is target atom type, q_shift is that atom's secondary shift
				// 1 = N
				// 2 = HA (HA3 for Gly)
				// 3 = C
				// 4 = CA
				// 5 = CB (HA2 for Gly)
				// 6 = HN
				core::Size q_shift_type(query_residue_shifts[d].first);
				core::Real q_shift(query_residue_shifts[d].second);

				//v_shift is the vall atom's secondary shift, v_sigma is the average deviation
				//on v_shifts for that type of atom at the vall residue's specific phi/psi location
				// (Think of v_shift as a phi/psi dependent and atom type dependent weight constant)
				VallResidueOP res = current_chunk->at(i);

				if ( res->secondary_shifts().size() < q_shift_type*2 ) {
					trCSScore.Debug << "Chunk has not enough secondary shifts to perform this query at position "
						<< i << " " << std::endl
						<< "pdb_id: " << current_chunk->get_pdb_id() << std::endl
						<< "chain_id: " << current_chunk->get_chain_id() << std::endl
						<< "sequence: " << current_chunk->get_sequence() << std::endl;
					continue;
				}

				core::Real v_shift(res->secondary_shifts()[(q_shift_type*2)-1]);
				//q_shift_type*2-1 because the array of 12 numbers goes shift1, sigma1, shift2, sigma2...
				core::Real v_sigma(res->secondary_shifts()[ q_shift_type*2 ]);

				//v_sigma is only 0.0 for atoms that don't exist in the vall. CB on glycine, for example.
				if ( v_sigma > 0.0 ) {

					core::Real sig_diff(std::abs((q_shift - v_shift) / v_sigma ));
					core::Real sigmoid_diff( 1 / ( 1 + exp((-a*sig_diff)+b) ) );

					tmp += sigmoid_diff;
					count += 1;
					//vall_data = true;  set but never used ~Labonte


					//THIS IS WHAT THE ORIGINAL CSROSETTA CS SCORE FUNCTION LOOKED LIKE:
					//core::Real c1_weight(1.0); //Reweight hydrogen and nitrogen values by 0.9
					//if ((q_shift_type == 1) || (q_shift_type == 6)) {// or (q_shift_type == 3)) {
					// c1_weight = 0.9;
					//}
					//core::Real diff(q_shift - v_shift);
					//if ( std::abs(diff) > (clip_factor*v_sigma) ) {
					// diff = clip_factor*v_sigma;
					//}
					//tmp += c1_weight*(diff/v_sigma)*(diff/v_sigma);
				}
			}

			//Arbitrarily high score for vall residues that don't have any CS data
			//(ie: residues immediately adjacent to missing density or termini)
			if ( ( count == 0 ) && ( query_residue_shifts.size() != 0 ) ) {
				tmp = 9999.9;
			} else {
				//Sigma6: scores don't use /N_shifts to normalize
				// This reweights each residue based on its number of shifts instead
				// so that over gaps in the data the CS score decreases in power and other scores
				// can take over.
				if ( count != 0 ) {
					tmp = ( tmp / count ) * query_residue_shifts.size();
				}
			}

			temp[i][r].first = tmp;
			temp[i][r].second = count;
		}
	}

	//runtime_assert(vall_data == true);//Make sure the vall had some chemical shift data in it

	scores_ = temp;

	trCSScore.Debug << "caching CS score for " << current_chunk->get_pdb_id()
		<< " of size " << current_chunk->size()
		<< ". The matrix is: "<<scores_.size()<<" x "<<scores_[1].size()<<std::endl;
}

bool CSScore::score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {
	return cached_score( fragment, scores );
}

bool CSScore::cached_score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {

	std::string & tmp = fragment->get_chunk()->chunk_key();

	if ( tmp.compare(cached_scores_id_) != 0 ) {
		do_caching(fragment->get_chunk());
		cached_scores_id_ = tmp;
	}

	//core::Size offset_q = fragment->get_first_index_in_query() - 1;
	//core::Size offset_v = fragment->get_first_index_in_vall() - 1;

	core::Real totalScore = 0.0;
	core::Real totalCount = 0.0;

	for ( core::Size i = 1; i <= fragment->get_length(); i++ ) {
		runtime_assert(fragment->get_first_index_in_vall() + i - 1 <= scores_.size());
		runtime_assert(fragment->get_first_index_in_query() + i - 1 <= scores_[1].size());


		std::pair< core::Real, core::Real> tmp = scores_[fragment->get_first_index_in_vall() + i - 1]
			[fragment->get_first_index_in_query() + i - 1];

		//tmp.first is the score for that residue comparison
		//tmp.second is the number of chemical shifts

		totalScore += tmp.first;
		totalCount += tmp.second;
	}

	// runtime_assert( totalScore != NULL );

	totalScore /= (core::Real) fragment->get_length();

	scores->set_score_component(totalScore, id_);

	if ( (totalScore < lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

void CSScore::clean_up() {
}

FragmentScoringMethodOP MakeCSScore::make(core::Size priority,
	core::Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP //picker
	, std::string // line
) {

	if ( option[in::file::talos_cs].user() ) {
		CSTalosIO in(option[in::file::talos_cs]());
		//  in.write(std::cerr);
		return (FragmentScoringMethodOP) FragmentScoringMethodOP( new CSScore(priority,
			lowest_acceptable_value, use_lowest,in) );
	}

	utility_exit_with_message(
		"Can't read CS data. Provide a chemical shifts file in TALOS format.");

	return nullptr;
}

} // scores
} // frag_picker
} // protocols
