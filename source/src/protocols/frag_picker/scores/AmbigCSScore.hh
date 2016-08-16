// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/AmbigCSScore.hh
/// @brief  Object that scores a fragment by its crmsd to the native
/// @author Robert Vernon

#ifndef INCLUDED_protocols_frag_picker_scores_AmbigCSScore_hh
#define INCLUDED_protocols_frag_picker_scores_AmbigCSScore_hh

// package headers
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/CSTalosIO.hh>

#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <utility/io/ozstream.hh>

// mini
#include <core/types.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

typedef utility::vector1<utility::vector1<Real> > Matrix;

//typedef utility::vector1<utility::vector1<std::pair<Real,Real> > > Score_Matrix;

/// @brief  scores a fragment by the root mean square deviation of Phi and Psi angles.
class AmbigCSScore: public CachingScoringMethod {
public:

	/// @brief  creates a Phi-Psi-based scoring function.
	/// @details Phi-Psi angles from a fragment will be compared to relevant angles in a given pose, which should have the same number of residues a the query sequence

	AmbigCSScore(Size, Real, bool, CSTalosIO&, CSTalosIO&);
	void do_caching(VallChunkOP);
	void clean_up();
	bool score(FragmentCandidateOP, FragmentScoreMapOP);
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);

	/// @brief  Print debugging informations about a score for a given fragment
	//bool describe_score(FragmentCandidateOP, FragmentScoreMapOP, std::ostream&);

private:

	utility::io::ozstream outfile_;

	//Shift Data. Should be one vector per residue, and residues without shift data should
	//have empty vectors.
	utility::vector1< utility::vector1< std::pair< Size, Real > > > target_Ashifts_;
	utility::vector1< utility::vector1< std::pair< Size, Real > > > target_Bshifts_;

	utility::vector1< utility::vector1< std::pair< Real, Real> > > scores_;

	std::string cached_scores_id_;

	utility::vector1<bool> existing_data_;
	void read_talos_phi_psi(std::string const&);
}; // AmbigCSScore

/// @brief  Maker class that produces a new AmbigCSScore object
class MakeAmbigCSScore: public MakeFragmentScoringMethod {
public:

	MakeAmbigCSScore() :
		MakeFragmentScoringMethod("AmbigCSScore") {
	}

	FragmentScoringMethodOP make(Size, Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols

#endif // INCLUDED_protocols_frag_picker_scores_AmbigCSScore_HH
