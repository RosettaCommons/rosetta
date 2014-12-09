// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/DisulfideIdentity.cc
/// @brief  a base class for fragment scoring
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/DisulfideIdentity.hh>

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

bool DisulfideIdentity::score(FragmentCandidateOP f,
		FragmentScoreMapOP empty_map) {

	// describe_score(f, empty_map, std::cerr);

	Real totalScore = 0;
	for (Size i = 1; i <= f->get_length(); i++) {


		//std::cout << "DISULF " << i << " " << f->get_first_index_in_query() << " " << query_[i + f->get_first_index_in_query() - 2] << std::endl;
		if (query_[i + f->get_first_index_in_query() - 2] == 'c') {

			VallResidueOP r = f->get_residue(i);
			utility::vector1<Real> tmplt_prof_row = r->profile();

			for (Size k = 1; k <= 20; k++){

				if ( k == 1 ) {
					totalScore += 1.0 - tmplt_prof_row[k];
				} else {
					totalScore += tmplt_prof_row[k];
				}
			}
		}



		//if (query_[i + f->get_first_index_in_query() - 2] == 'c') {
		//	if (r->aa() != 'c')
		//		totalScore += 1;
		//} else {
		//	if (r->aa() == 'c')
		//		totalScore += 1;
		//}
	}
	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
		return false;
	return true;
}

bool DisulfideIdentity::describe_score(FragmentCandidateOP f,
		FragmentScoreMapOP empty_map, std::ostream& out) {

	Real totalScore = 0;

	out << f->get_chunk()->get_pdb_id() << "  " << I(5,
			f->get_first_index_in_vall()) << " ";
	for (Size i = 1; i <= f->get_length(); i++)
		out << f->get_residue(i)->aa();
	out << std::endl << "            ";
	for (Size i = 1; i <= f->get_length(); i++) {
		VallResidueOP r = f->get_residue(i);
		if (query_[i + f->get_first_index_in_query() - 2] == 'c') {

			if (r->aa() != query_[i + f->get_first_index_in_query() - 2]) {
				totalScore += 1;
				out << "-";
			} else {
				if (r->aa() == 'c')
					totalScore += 1;
				out << "-";
			}

		}
	}
	out << "\nquery " << I(5, f->get_first_index_in_query()) << " ";
	for (Size i = 1; i <= f->get_length(); i++)
		out << query_[i + f->get_first_index_in_query() - 2];
	if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
		out << "\nTotal score " << F(5, 3, totalScore) << " ACCEPTED"
				<< std::endl;
	else
		out << "\nTotal score " << F(5, 3, totalScore) << " REJECTED"
				<< std::endl;

	totalScore /= (Real) f->get_length();
	empty_map->set_score_component(totalScore, id_);
	if ((totalScore > lowest_acceptable_value_) && (use_lowest_ == true))
		return false;
	return true;
}

} // scores
} // frag_picker
} // protocols


