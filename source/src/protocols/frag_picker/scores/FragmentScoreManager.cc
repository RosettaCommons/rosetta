// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/FragmentScoreManager.hh
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


// package headers
#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>

#include <protocols/frag_picker/scores/SequenceIdentity.hh>
#include <protocols/frag_picker/scores/ConstScore.hh>
#include <protocols/frag_picker/scores/BFactor.hh>
#include <protocols/frag_picker/scores/DisulfideIdentity.hh>
#include <protocols/frag_picker/scores/DisulfideDistance.hh>
#include <protocols/frag_picker/scores/PCS_FragDistance.hh>
#include <protocols/frag_picker/scores/ProlinePhiScore.hh>
#include <protocols/frag_picker/scores/RamaScore.hh>
#include <protocols/frag_picker/scores/RDCScore.hh>
#include <protocols/frag_picker/scores/SecondaryIdentity.hh>
#include <protocols/frag_picker/scores/SecondarySimilarity.hh>
#include <protocols/frag_picker/scores/TalosSSSimilarity.hh>
#include <protocols/frag_picker/scores/PartialSecondarySimilarity.hh>
#include <protocols/frag_picker/scores/TorsionBinSimilarity.hh>
#include <protocols/frag_picker/scores/ProfileScore.hh>
#include <protocols/frag_picker/scores/ProfileScoreL1.hh>
#include <protocols/frag_picker/scores/ProfileScoreBlosum62.hh>
#include <protocols/frag_picker/scores/ProfileScoreSubMatrix.hh>
#include <protocols/frag_picker/scores/ProfileScoreDistWeight.hh>
#include <protocols/frag_picker/scores/FragmentCrmsd.hh>
#include <protocols/frag_picker/scores/FragmentChunkCrms.hh>
#include <protocols/frag_picker/scores/FragmentAllAtomCrmsd.hh>
#include <protocols/frag_picker/scores/MidPsiOut.hh>
#include <protocols/frag_picker/scores/MidPhiOut.hh>
#include <protocols/frag_picker/scores/PhiPsiRmsd.hh>
#include <protocols/frag_picker/scores/JCoupling.hh>
#include <protocols/frag_picker/scores/PhiPsiSquareWell.hh>
#include <protocols/frag_picker/scores/AtomPairConstraintsScore.hh>
#include <protocols/frag_picker/scores/DihedralConstraintsScore.hh>
#include <protocols/frag_picker/scores/InterbondAngleScore.hh>
#include <protocols/frag_picker/scores/CSScore.hh>
#include <protocols/frag_picker/scores/AmbigCSScore.hh>
#include <protocols/frag_picker/scores/GunnCostScore.hh>
#include <protocols/frag_picker/scores/ABEGO_SS_Score.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/scores/FragmentDME.hh>
#include <protocols/frag_picker/scores/HydrophobicitySimilarity.hh>
#include <protocols/frag_picker/scores/HydrophobicityProfileSimilarity.hh>
#include <protocols/frag_picker/scores/SolventAccessibility.hh>
#include <protocols/frag_picker/scores/Phi.hh>
#include <protocols/frag_picker/scores/Psi.hh>
#include <protocols/frag_picker/scores/ProfileScoreStructL1.hh>
#include <protocols/frag_picker/scores/FragmentCrmsdResDepth.hh>

#include <basic/Tracer.hh>

#include <utility/io/izstream.hh>

// C++
#include <algorithm>
#include <map>
#include <utility>

#include <protocols/frag_picker/CommonFragmentComparators.hh>
#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer tr( "protocols.frag_picker.FragmentScoreManager" );

namespace protocols {
namespace frag_picker {
namespace scores {

/// @details Auto-generated virtual destructor
FragmentScoreManager::~FragmentScoreManager() {}

/// @brief creates an empty score map
/// @details this is the recommended way to create FragmentScoreMap objects since FragmentScoreManager knows exactly
/// what is the correct size of the map i.e. how many scoring terms have been registered.
FragmentScoreMapOP FragmentScoreManager::create_empty_map() {

	return FragmentScoreMapOP( new FragmentScoreMap(score_weights_.size()) );
}


bool sort_scores(FragmentScoringMethodOP elem1, FragmentScoringMethodOP elem2) {
	return elem1->get_priority() > elem2->get_priority();
}

void FragmentScoreManager::add_scoring_method(
	FragmentScoringMethodOP scoring_term,
	core::Real weight) {

	std::map<FragmentScoringMethodOP, core::Real> weights;
	scores_.push_back(scoring_term);
	score_weights_.push_back(weight);

	for ( core::Size i = 1; i <= scores_.size(); i++ ) {
		std::pair<FragmentScoringMethodOP, core::Real> p(scores_[i],
			score_weights_[i]);
		weights.insert(p);
	}

	sort(scores_.begin(), scores_.end(), sort_scores);

	for ( core::Size i = 1; i <= scores_.size(); i++ ) {
		scores_[i]->set_id(i);
		score_weights_[i] = weights.find(scores_[i])->second;
	}
	width_.insert(std::pair<FragmentScoringMethodOP, core::Size>(scoring_term,
		default_width_));
	precision_.insert(std::pair<FragmentScoringMethodOP, core::Size>(scoring_term,
		default_precision_));
}

core::Real FragmentScoreManager::total_score(FragmentScoreMapOP f) {

	if ( !f->was_modified() ) {
		return f->get_most_recent_total_score();
	}

	core::Real total = 0;
	utility::vector1<core::Real> &s = f->get_score_components();

	debug_assert (s.size() == score_weights_.size());
	for ( core::Size i = 1; i <= score_weights_.size(); i++ ) {
		total += score_weights_[i] * s[i];
	}
	f->recent_total_ = total;
	f->was_modified_ = false;

	return total;
}

void FragmentScoreManager::do_caching(VallChunkOP chunk) {

	for ( core::Size iScore = 1; iScore <= scores_.size(); ++iScore ) {
		if ( ( zeros_score_later_ ) && ( fabs(score_weights_[iScore]) < 0.000001 ) ) continue;
		CachingScoringMethodOP score =
			utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::CachingScoringMethod > ( scores_[iScore] );
		if ( score != 0 ) {
			score->do_caching(chunk);
		}
	}
}

void FragmentScoreManager::clean_up() {

	for ( core::Size iScore = 1; iScore <= scores_.size(); ++iScore ) {
		CachingScoringMethodOP score =
			utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::CachingScoringMethod > ( scores_[iScore] );
		if ( score != 0 ) {
			score->clean_up();
		}
	}
}

void FragmentScoreManager::show_scoring_methods(std::ostream& out) {

	out
		<< "Fragment scoring scheme used:\n-----------------------------------------\n";
	out << "id    fragment score method name   weight\n";
	for ( core::Size i = 1; i <= scores_.size(); i++ ) {
		out << std::setw(3) << scores_[i]->get_id() << " " << std::setw(30)
			<< scores_[i]->get_score_name() << std::setw(5)
			<< score_weights_[i] << "\n";
	}
	out << "-----------------------------------------\n" << std::endl;
}

void FragmentScoreManager::describe_fragments(utility::vector1<std::pair<
	FragmentCandidateOP, FragmentScoreMapOP> > const& pairs,
	std::ostream& out) {

	using namespace ObjexxFCL::format;
	if ( pairs.size() == 0 ) return;
	bool if_quota = false;
	if ( pairs[1].second->get_quota_score() < 999.98 ) if_quota = true;
	out << "#" << RJ(10, "query_pos ");
	out << RJ(10, "vall_pos ");
	out << RJ(6, "pdbid");
	out << " c ss ";
	utility::vector1<core::Size> w(scores_.size()+4);
	for ( core::Size i = 1; i <= scores_.size(); i++ ) {
		w[i] = width_.find(scores_[i])->second;
		out << " " << scores_[i]->get_score_name();
		if ( scores_[i]->get_score_name().length() > w[i] ) {
			w[i] = scores_[i]->get_score_name().length();
		}
	}
	if ( if_quota ) {
		w[scores_.size()+1] = 10; // 10 characters for the total quota score
		w[scores_.size()+2] = 9; // 9 characters for the total score
		w[scores_.size()+3] = 10; // 10 characters for the pool name
		out << " QUOTA_TOT TOTAL POOL_NAME FRAG_ID"<<std::endl;
	} else {
		w[scores_.size()+1] = 9; // 9 characters for the total score
		out << "  TOTAL  FRAG_ID"<<std::endl;
	}
	for ( core::Size iF = 1; iF <= pairs.size(); ++iF ) {
		if ( !pairs[iF].first || !pairs[iF].second ) {
			tr.Warning << "final_frag candidate " << iF << " is corrupted. skipping... " << std::endl;
			continue;
		}
		FragmentCandidateOP fr = pairs[iF].first;
		FragmentScoreMapOP sc = pairs[iF].second;
		VallResidueOP r = fr->get_residue(1);
		out << " " << I(10, fr->get_first_index_in_query());
		//  out << " " << I(10, fr->get_first_index_in_vall());
		out << " " << I(10, r->resi());
		out << " " << RJ(5, fr->get_pdb_id());
		out << " " << fr->get_chain_id();
		out << " " << fr->get_middle_ss();

		for ( core::Size i = 1; i <= scores_.size(); i++ ) {
			core::Size p = precision_.find(scores_[i])->second;
			out << " " << F(w[i], p, sc->get_score_components()[i]);
		}
		if ( if_quota ) {
			out << " " << F(w[scores_.size()+1],TOTAL_PRECISION,sc->get_quota_score());
			out << " " << F(w[scores_.size()+1],TOTAL_PRECISION,total_score(sc));
			out << " " << std::setw(w[scores_.size()+3])<<fr->get_pool_name();
		} else {
			out << F(w[scores_.size()+1],TOTAL_PRECISION,total_score(sc));
		}

		debug_assert ( fr->key() > 0 );
		debug_assert ( fr->key() < 4000000 ); // Put your Vall's size here

		out << I(10, fr->key() ) << std::endl;

	}
}

bool FragmentScoreManager::score_fragment(FragmentCandidateOP candidate,
	FragmentScoreMapOP empty_map) {

	for ( core::Size iScore = 1; iScore <= scores_.size(); iScore++ ) {
		if ( ( zeros_score_later_ ) && ( fabs(score_weights_[iScore]) < 0.000001 ) ) {
			continue;
		}
		if ( !scores_[iScore]->score(candidate, empty_map) ) {
			return false;
		}
	}
	return true;
}

bool FragmentScoreManager::score_zero_scores(
	FragmentCandidateOP candidate,
	FragmentScoreMapOP empty_map)
{
	for ( core::Size iScore = 1; iScore <= scores_.size(); iScore++ ) {
		if ( fabs(score_weights_[iScore]) < 0.000001 ) {
			scores_[iScore]->score(candidate, empty_map);
		}
	}
	return true;
}


bool FragmentScoreManager::score_fragment_from_cache(
	FragmentCandidateOP candidate,
	FragmentScoreMapOP empty_map)
{
	for ( core::Size iScore = 1; iScore <= scores_.size(); iScore++ ) {
		if ( ( zeros_score_later_ ) && ( fabs(score_weights_[iScore]) < 0.000001 ) ) {
			continue;
		}
		CachingScoringMethodOP s =
			utility::pointer::dynamic_pointer_cast< protocols::frag_picker::scores::CachingScoringMethod > ( scores_[iScore] );
		if ( s != 0 ) {
			if ( !s->cached_score(candidate, empty_map) ) {
				return false;
			}
		} else {
			if ( !scores_[iScore]->score(candidate, empty_map) ) {
				return false;
			}
		}
	}
	return true;
}

void FragmentScoreManager::register_score_maker(
	MakeFragmentScoringMethodOP scoring_term_maker) {
	registered_makers_[scoring_term_maker->get_score_name()]
		= scoring_term_maker;
}

FragmentScoreManager::FragmentScoreManager() {

	default_precision_ = 2;
	default_width_ = 6;
	zeros_score_later_ = true;

	register_score_maker(MakeFragmentScoringMethodOP( new MakeSecondaryIdentity() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeSequenceIdentity() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeBFactor() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeDisulfideIdentity() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeDisulfideDistance() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakePCS_FragDistance() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeProlinePhiScore() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeRamaScore() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeRDCScore() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeProfileScore() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeProfileScoreL1() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeProfileScoreBlosum62() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeProfileScoreSubMatrix() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeProfileScoreDistWeight() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeSecondarySimilarity() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeTalosSSSimilarity() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakePartialSecondarySimilarity() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeTorsionBinSimilarity() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeFragmentCrmsd() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeFragmentChunkCrms() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeMidPhiOut() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeMidPsiOut() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakePhiPsiRmsd() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakePhiPsiSquareWell() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeJCoupling() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeAtomPairConstraintsScore() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeDihedralConstraintsScore() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeCSScore() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeAmbigCSScore() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeTorsionBinSimilarity() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeInterbondAngleScore() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeGunnCostScore() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeABEGO_SS_Score() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeFragmentAllAtomCrmsd() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeConstScore() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeFragmentDME() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeHydrophobicitySimilarity() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeHydrophobicityProfileSimilarity() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeSolventAccessibility() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakePhi() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakePsi() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeProfileScoreStructL1() ));
	register_score_maker(MakeFragmentScoringMethodOP( new MakeFragmentCrmsdResDepth() ));

}

void FragmentScoreManager::create_scoring_method(
	std::string const & score_name, core::Size priority, core::Real weight,
	core::Real lowest, bool use_lowest, FragmentPickerOP picker, std::string config_line) {

	std::map<std::string, MakeFragmentScoringMethodOP>::iterator m =
		registered_makers_.find(score_name);
	if ( m == registered_makers_.end() ) {
		utility_exit_with_message("[ERROR]: unknown score type " + score_name
			+ "!");
	}

	add_scoring_method(m->second->make(priority, lowest, use_lowest, picker, config_line),
		weight);
}

void FragmentScoreManager::create_scores(std::string const & file_name,
	FragmentPickerOP picker) {

	utility::io::izstream input(file_name.c_str());
	if ( !input ) {
		utility_exit_with_message("[ERROR]: can't open file " + file_name + "!");
	}

	std::string line;
	std::string score_name;
	core::Size priority;
	core::Real weight;
	std::string low_string;
	core::Real lowest;
	bool use_lowest;
	std::string remark("");

	while ( getline(input, line) ) {
		if ( line[0] == '#' ) {
			continue;
		}

		std::istringstream line_stream(line);
		line_stream >> score_name >> priority >> weight >> low_string;

		if ( low_string == "-" ) {
			lowest = 0.0;
			use_lowest = false;
		} else {
			std::istringstream low_stream(low_string);
			//lowest = static_cast<core::Real>(low_string);
			low_stream >> lowest;
			use_lowest = true;
		}

		remark = "";
		if ( line_stream.good() ) {
			line_stream >> remark;
		}

		//  std::cout << line << std::endl;
		//  std::cout << score_name << " " << priority << " " << weight << " " << low_string << " " << use_lowest << " " << remark << std::endl;

		create_scoring_method(score_name, priority, weight, lowest, use_lowest, picker, remark);
	}
}

} // scores
} // frag_picker
} // protocols


