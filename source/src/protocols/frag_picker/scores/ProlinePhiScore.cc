// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/ProlinePhi.cc
/// @brief  a base class for fragment scoring
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/ProlinePhiScore.hh>

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

using namespace ObjexxFCL::format;

bool ProlinePhiScore::score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map) {

	// describe_score(f, empty_map, std::cerr);

	core::Real totalScore = 0;
	for ( core::Size i = 1; i <= f->get_length(); i++ ) {
		VallResidueOP r = f->get_residue(i);

		//std::cout << "DISULF " << i << " " << f->get_first_index_in_query() << " " << query_[i + f->get_first_index_in_query() - 2] << std::endl;

		core::Size query_res = i + f->get_first_index_in_query() - 2;

		core::Real phi = r->phi();
		core::Real psi = r->psi();

		if ( query_[query_res] == 'P' ) {
			if ( (phi < -103) || (phi > -33) ) {
				return false; // Never possible
			}

			if ( phi < -59.0 ) {
				totalScore += 0.00197753 * (phi + 55.7777)*(phi + 55.7777);
			} else {
				totalScore += 6.08188 * (-0.282087*phi)
					* (-0.0133411*phi)*(-0.0133411*phi)
					* (-0.000115456*phi)*(-0.000115456*phi)*(-0.000115456*phi);
			}
		}

		if ( query_[query_res] != 'G' ) {

			if ( phi > 75 ) {
				return false;
			}

			if ( (psi < -75) && (psi > -170) ) {
				return false;
			}
		}

		if ( (query_[query_res] == 'I')
				|| (query_[query_res] == 'V')
				|| (query_[query_res] == 'T') ) {
			if ( phi > -40 ) {
				return false;
			}
		}

		if ( (query_res + 1) < query_.size() ) {
			if ( (query_[query_res+1] == 'P') && (query_[query_res] != 'G') ) {
				if ( (psi < 40) && ((psi > -25) || (phi < -90)) ) {
					return false;
				}
			}
		}

	}
	totalScore /= (core::Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

bool ProlinePhiScore::describe_score(FragmentCandidateOP f,
	FragmentScoreMapOP empty_map, std::ostream& out) {

	core::Real totalScore = 0;

	out << f->get_chunk()->get_pdb_id() << "  " << I(5,
		f->get_first_index_in_vall()) << " ";
	for ( core::Size i = 1; i <= f->get_length(); i++ ) {
		out << f->get_residue(i)->aa();
	}
	out << std::endl << "            ";
	for ( core::Size i = 1; i <= f->get_length(); i++ ) {
		VallResidueOP r = f->get_residue(i);

		if ( query_[i + f->get_first_index_in_query() - 2] == 'P' ) {
			core::Real phi = r->phi();

			if ( (phi < -120) || (phi > -10) ) {
				return false; // Never possible
			}

			if ( phi < -59.0 ) {
				totalScore += 0.00197753 * (phi + 55.7777)*(phi + 55.7777);
			} else {
				totalScore += 6.08188 * (-0.282087*phi)
					* (-0.0133411*phi)*(-0.0133411*phi)
					* (-0.000115456*phi)*(-0.000115456*phi)*(-0.000115456*phi);
			}

		}
	}
	out << "\nquery " << I(5, f->get_first_index_in_query()) << " ";
	for ( core::Size i = 1; i <= f->get_length(); i++ ) {
		out << query_[i + f->get_first_index_in_query() - 2];
	}
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		out << "\nTotal score " << F(5, 3, totalScore) << " ACCEPTED"
			<< std::endl;
	} else {
		out << "\nTotal score " << F(5, 3, totalScore) << " REJECTED"
			<< std::endl;
	}

	totalScore /= (core::Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ( (totalScore > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

} // scores
} // frag_picker
} // protocols


