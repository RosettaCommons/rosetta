// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/AmbigCSScore.cc
/// @brief  Object that scores a fragment by target-observed/vall-predicted chemical shift distances
/// @author Robert Vernon rvernon@u.washington.edu

#include <protocols/frag_picker/scores/AmbigCSScore.hh>

#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/CSTalosIO.hh>
// AUTO-REMOVED #include <protocols/frag_picker/scores/FragmentCrmsd.hh>
// AUTO-REMOVED #include <protocols/frag_picker/scores/FragmentScoreManager.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
// AUTO-REMOVED #include <protocols/frag_picker/VallProvider.hh>
#include <utility/io/ozstream.hh>


// option key includes
// AUTO-REMOVED #include <core/init/init.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>

// AUTO-REMOVED #include <basic/prof.hh>

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

static thread_local basic::Tracer trAmbigCSScore(
		"protocols.frag_picker.scores.AmbigCSScore");


//AmbigCSScore Constructor
// The Talos file reader is passed in as an object, and secondary shifts are calculated during AmbigCSScore construction
// (Secondary shifts are shift deviations from random coil, where random coil values are defined according to
// the combination of atom type, residue type, previous residue type, and next residue type.
AmbigCSScore::AmbigCSScore(Size priority, Real lowest_acceptable_value, bool use_lowest,
								 CSTalosIO& readerA, CSTalosIO& readerB) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "AmbigCSScore")
{

	//outfile_ = utility::io::ozstream tmp("allcomparisons.out");
	//outfile_.open("allcomparisons.out");

	trAmbigCSScore << "READING SHIFTS!" << std::endl;
	CS2ndShift secondary_shift_calculatorA(readerA, false);
	CS2ndShift secondary_shift_calculatorB(readerB, false);
	trAmbigCSScore << "SHOULD BE DONE WRITING 2nd SHIFTS" << std::endl;

	target_Ashifts_ = secondary_shift_calculatorA.shifts();
	target_Bshifts_ = secondary_shift_calculatorB.shifts();
}

//Residue-Residue scores are cached for current vall chunk
//This is where the AmbigCSScore equation lives
void AmbigCSScore::do_caching(VallChunkOP current_chunk) {

	//ONLY USED IN THE OLD VERSION OF THE SCORE
	//clip_factor is used to define the maximum shift difference
	//larger differences are adjusted down to clip_factor*v_shift
	//Real const clip_factor(3.0);

	//bool vall_data(false);
	trAmbigCSScore << "caching CS score for " << current_chunk->get_pdb_id()
			<< " of size " << current_chunk->size() << std::endl;

	//Check to see if the cache needs to be recalculated
	std::string & tmp = current_chunk->chunk_key();
	if (tmp.compare(cached_scores_id_) == 0)
		return;
	cached_scores_id_ = tmp;

	//Initialize empty 2D table, vall-length x target-length
	Size query_sequence_length = target_Ashifts_.size();
	runtime_assert(query_sequence_length == target_Bshifts_.size());
	std::pair< Real, Real > empty(0,0);
	utility::vector1< utility::vector1< std::pair< Real, Real> > > temp( current_chunk->size(),
	utility::vector1<std::pair< Real, Real> > (query_sequence_length, empty ) );
	runtime_assert( target_Ashifts_.size() > 0 );

	//SIGMOID CONSTANTS - Should be set in constructor, not command line flags
	Real a( option[frags::sigmoid_cs_A]() ); // default = 4
	Real b( option[frags::sigmoid_cs_B]() ); // default = 5

	//Loop logic is "For each target x vall residue comparison, sum up total of
	//all shift differences"
	for (Size r = 1; r <= target_Ashifts_.size(); ++r) {
		utility::vector1< std::pair< Size, Real > > query_residue_shiftsA(target_Ashifts_[r]);
		utility::vector1< std::pair< Size, Real > > query_residue_shiftsB(target_Bshifts_[r]);

		if (query_residue_shiftsA.size() != query_residue_shiftsB.size()) {
			utility_exit_with_message("ERROR: -in::file::ambig_talos_cs_A file does not have the same number of shifts as -in::file::ambig_talos_cs_B file, check your formatting, aside from the shifts themselves the files must be identical");
		}

		for (Size i = 1; i <= current_chunk->size(); ++i) {
			Real tmp = 0.0;
			Real count = 0.0;
			for (Size d = 1; d <= query_residue_shiftsA.size(); ++d) {

				//q_shift_type is target atom type, q_shift is that atom's secondary shift
				Size q_shift_typeA(query_residue_shiftsA[d].first);
				Real q_shiftA(query_residue_shiftsA[d].second);

				Size q_shift_typeB(query_residue_shiftsB[d].first);
				Real q_shiftB(query_residue_shiftsB[d].second);

				if (q_shift_typeA != q_shift_typeB) {
					utility_exit_with_message("ERROR: -in::file::ambig_talos_cs_A file does not match -in::file::ambig_talos_cs_B file, check your formatting, aside from the shifts themselves the files must be identical, even the order matters");
				}

				//v_shift is the vall atom's secondary shift, v_sigma is the average deviation
				//on v_shifts for that type of atom at the vall residue's specific phi/psi location
				// (Think of v_shift as a phi/psi dependent and atom type dependent weight constant)
				VallResidueOP res = current_chunk->at(i);
				Real v_shift(res->secondary_shifts()[(q_shift_typeA*2)-1]);
				Real v_sigma(res->secondary_shifts()[ q_shift_typeA*2 ]);

				//v_sigma is only 0.0 for atoms that don't exist in the vall. CB on glycine, for example.
				if (v_sigma > 0.0) {

					Real sig_diffA(std::abs((q_shiftA - v_shift) / v_sigma ));
					Real sig_diffB(std::abs((q_shiftB - v_shift) / v_sigma ));

					//Always use the lowest diff of the two.
					if ( sig_diffB < sig_diffA ) {
						sig_diffA = sig_diffB;
					}

					Real sigmoid_diff( 1 / ( 1 + exp((-a*sig_diffA)+b) ) );


					tmp += sigmoid_diff;
					count += 1;
					//vall_data = true;  set but never used ~Labonte


					//THIS IS WHAT THE ORIGINAL CSROSETTA CS SCORE FUNCTION LOOKED LIKE:
					//Real c1_weight(1.0); //Reweight hydrogen and nitrogen values by 0.9
					//if ((q_shift_type == 1) || (q_shift_type == 6)) {// or (q_shift_type == 3)) {
					//	c1_weight = 0.9;
					//}
					//Real diff(q_shift - v_shift);
					//if ( std::abs(diff) > (clip_factor*v_sigma) ) {
						//	diff = clip_factor*v_sigma;
						//}
					//tmp += c1_weight*(diff/v_sigma)*(diff/v_sigma);
				}
			}

			//Arbitrarily high score for vall residues that don't have any CS data
			//(ie: residues immediately adjacent to missing density or termini)
			if ( ( count == 0 ) && ( query_residue_shiftsA.size() != 0 ) ) {
				tmp = 9999.9;
			} else {
				//Sigma6: scores don't use /N_shifts to normalize
				// This reweights each residue based on its number of shifts instead
				// so that over gaps in the data the CS score decreases in power and other scores
				// can take over.
				if ( count != 0 )
					tmp = ( tmp / count ) * query_residue_shiftsA.size();
			}

			temp[i][r].first = tmp;
			temp[i][r].second = count;
		}
	}

	//runtime_assert(vall_data == true);//Make sure the vall had some chemical shift data in it

	scores_ = temp;

	trAmbigCSScore << "caching CS score for " << current_chunk->get_pdb_id()
	                << " of size " << current_chunk->size()
	                << ". The matrix is: "<<scores_.size()<<" x "<<scores_[1].size()<<std::endl;
}

bool AmbigCSScore::score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {
	return cached_score( fragment, scores );
}

bool AmbigCSScore::cached_score(FragmentCandidateOP fragment, FragmentScoreMapOP scores) {

	std::string & tmp = fragment->get_chunk()->chunk_key();

	if (tmp.compare(cached_scores_id_) != 0) {
		do_caching(fragment->get_chunk());
		cached_scores_id_ = tmp;
	}

	//Size offset_q = fragment->get_first_index_in_query() - 1;
	//Size offset_v = fragment->get_first_index_in_vall() - 1;

	Real totalScore = 0.0;
	Real totalCount = 0.0;

	for (Size i = 1; i <= fragment->get_length(); i++) {
		runtime_assert(fragment->get_first_index_in_vall()	+ i - 1 <= scores_.size());
		runtime_assert(fragment->get_first_index_in_query() + i - 1 <= scores_[1].size());


		std::pair< Real, Real> tmp = scores_[fragment->get_first_index_in_vall() + i - 1]
			                [fragment->get_first_index_in_query()	+ i - 1];

		//tmp.first is the score for that residue comparison
		//tmp.second is the number of chemical shifts

		totalScore += tmp.first;
		totalCount += tmp.second;
	}

//	runtime_assert( totalScore != NULL );

	scores->set_score_component(totalScore, id_);

	if ((totalScore < lowest_acceptable_value_) && (use_lowest_ == true))
		return false;
	return true;
}

void AmbigCSScore::clean_up() {
}

FragmentScoringMethodOP MakeAmbigCSScore::make(Size priority,
		Real lowest_acceptable_value, bool use_lowest, FragmentPickerOP //picker
		, std::string // line
) {

	if (option[in::file::ambig_talos_cs_A].user() &&
			option[in::file::ambig_talos_cs_B].user()) {
	  CSTalosIO inA(option[in::file::ambig_talos_cs_A]());
		CSTalosIO inB(option[in::file::ambig_talos_cs_B]());
	  inA.write(std::cerr);
	  inB.write(std::cerr);
	  return (FragmentScoringMethodOP) FragmentScoringMethodOP( new AmbigCSScore(priority,
																								 lowest_acceptable_value,
																								 use_lowest,inA,inB) );
	}

	utility_exit_with_message(
			"Can't read ambiguous CS data. Provide two chemical shifts file in TALOS format using flags -in::file::ambig_talos_cs_A and in::file::ambig_talos_cs_B");

	return NULL;
}

} // scores
} // frag_picker
} // protocols
