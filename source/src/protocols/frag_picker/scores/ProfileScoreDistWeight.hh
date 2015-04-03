// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/ProfileScore.hh
/// @brief  scores a fragment by weighting L1 profile distances by residue type
/// @author Robert Vernon

#ifndef INCLUDED_protocols_frag_picker_scores_ProfileScoreDistWeight_hh
#define INCLUDED_protocols_frag_picker_scores_ProfileScoreDistWeight_hh

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/scores/CachingScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <core/fragment/SecondaryStructure.hh>
#include <core/sequence/SequenceProfile.hh>

// utility headers
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <string>
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

typedef utility::vector1<utility::vector1<Real> > Matrix;

/// @brief  a fragment candidate
class ProfileScoreDistWeight: public CachingScoringMethod {
public:

	ProfileScoreDistWeight(
		Size priority,
		Real lowest_acceptable_value,
		bool use_lowest,
		sequence::SequenceProfileOP query_profile,
		core::fragment::SecondaryStructureOP query_ss_prediction,
		std::string query_sequence, Size longest_vall_chunk
	) :
		CachingScoringMethod(
			priority, lowest_acceptable_value, use_lowest, "ProfileScoreDistWeight"
	) {

		query_profile_ = query_profile;
		query_ss_ = query_ss_prediction;
		query_sequence_ = query_sequence;

		for (Size i = 1; i <= query_profile->length(); ++i) {
			utility::vector1<Real> row(longest_vall_chunk);
			scores_.push_back(row);
		}


		utility::vector1< utility::vector1< utility::vector1 <Real> > > temp(
			3, utility::vector1< utility::vector1 <Real> > (
				20, utility::vector1<Real> (
					20, 0.0
				)
			)
		);
		distance_weights_ = temp;

		utility::io::izstream data("distances.txt");

		if (!data)
			utility_exit_with_message("[ERROR] Unable to open distances.txt");

		std::string line;
		char res1;
		char res2;
		char ss1;
		char ss2;
		Real dist;

		std::map<char,Size> ss_type_temp;
		ss_type_temp.insert(std::make_pair('H',1));
		ss_type_temp.insert(std::make_pair('E',2));
		ss_type_temp.insert(std::make_pair('L',3));

		ss_type_map_ = ss_type_temp;

		std::map<char,Size> aa_order_tmp;
		aa_order_tmp.insert(std::make_pair('A',1));
		aa_order_tmp.insert(std::make_pair('C',2));
		aa_order_tmp.insert(std::make_pair('D',3));
		aa_order_tmp.insert(std::make_pair('E',4));
		aa_order_tmp.insert(std::make_pair('F',5));
		aa_order_tmp.insert(std::make_pair('G',6));
		aa_order_tmp.insert(std::make_pair('H',7));
		aa_order_tmp.insert(std::make_pair('I',8));
		aa_order_tmp.insert(std::make_pair('K',9));
		aa_order_tmp.insert(std::make_pair('L',10));
		aa_order_tmp.insert(std::make_pair('M',11));
		aa_order_tmp.insert(std::make_pair('N',12));
		aa_order_tmp.insert(std::make_pair('P',13));
		aa_order_tmp.insert(std::make_pair('Q',14));
		aa_order_tmp.insert(std::make_pair('R',15));
		aa_order_tmp.insert(std::make_pair('S',16));
		aa_order_tmp.insert(std::make_pair('T',17));
		aa_order_tmp.insert(std::make_pair('V',18));
		aa_order_tmp.insert(std::make_pair('W',19));
		aa_order_tmp.insert(std::make_pair('Y',20));

		aa_order_map_ = aa_order_tmp;

		while (getline(data, line)) {
			std::istringstream line_stream(line);
			line_stream >> res1 >> ss1 >> res2 >> ss2 >> dist;


			Size res_type( 0 );
			if (ss1 == ss2) {
				res_type = ss_type_map_.find(ss1)->second;
			}

			//std::cout << "AAAA: " << res1 << " " << res2 << " " << ss1 << " " << ss2 << " " << dist << " " << res_type << std::endl;

			if ( res_type != 0 ) {


				Size i_res1, i_res2;

				i_res1 = aa_order_map_.find(res1)->second;
				i_res2 = aa_order_map_.find(res2)->second;

				distance_weights_[res_type][i_res1][i_res2] = dist;
				//std::cout << "BBBB: " << res1 << " " << res2 << " " << ss1 << " " << ss2 << " " << res_type << " " << i_res1 << " " << i_res2 << " " << dist << std::endl;
			}
		}
	}

	void do_caching(VallChunkOP);
	void clean_up() {
	}
	bool cached_score(FragmentCandidateOP, FragmentScoreMapOP);
	//bool describe_score(FragmentCandidateOP f, FragmentScoreMapOP empty_map, std::ostream& out);

protected:
	Matrix scores_;

private:
	utility::vector1< utility::vector1< utility::vector1 <Real> > > distance_weights_;

	std::string query_sequence_;

	std::map<char,Size> ss_type_map_;
	std::map<char,Size> aa_order_map_;

	sequence::SequenceProfileOP query_profile_;
	core::fragment::SecondaryStructureOP query_ss_;
	std::string cached_scores_id_;
	void clear();
};

class MakeProfileScoreDistWeight: public MakeFragmentScoringMethod {
public:

	MakeProfileScoreDistWeight() :
		MakeFragmentScoringMethod("ProfileScoreDistWeight") {
	}

	FragmentScoringMethodOP make(Size, Real, bool, FragmentPickerOP, std::string);
};

} // scores
} // frag_picker
} // protocols

#endif /* INCLUDED_protocols_frag_picker_scores_ProfileScoreDistWeight_HH */
