// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/scores/Rama.cc
/// @brief  Ramachandran surface scores for fragment picking
/// @author Robert Vernon (rvernon@u.washington.edu)

#include <protocols/frag_picker/scores/RamaScore.hh>

// type headers
#include <core/types.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>

// utility headers
#include <basic/database/open.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <ObjexxFCL/string.functions.hh>

namespace protocols {
namespace frag_picker {
namespace scores {

using namespace basic::options;
using namespace basic::options::OptionKeys;

RamaScore::RamaScore(Size priority, Real lowest_acceptable_value, bool use_lowest, std::string & fastaQuery, std::string prediction_name ) :
	CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "RamaScore"), query_(fastaQuery),prediction_name_(prediction_name) {

	core::fragment::SecondaryStructureOP default_ss;
	std::cout << "QUERY" << fastaQuery << " " << query_.size() << std::endl;
	default_ss->extend(query_.size());

	for ( Size i = 1; i <= query_.size(); ++i ) {
		default_ss->set_fractions(i, 1.0, 1.0, 1.0 );
	}

	query_ss_ = default_ss;

	SetupRamaTables();
}


RamaScore::RamaScore(Size priority, Real lowest_acceptable_value, bool use_lowest, std::string & fastaQuery, core::fragment::SecondaryStructureOP query_prediction, std::string prediction_name ) :
		CachingScoringMethod(priority, lowest_acceptable_value, use_lowest, "RamaScore"), query_(fastaQuery), query_ss_(query_prediction), prediction_name_(prediction_name)
{
	SetupRamaTables();
}

void
RamaScore::SetupRamaTables()
{

	// skip if static data is already set up
	if (sequence_rama_tables_.size() > 1) return;

	utility::vector1< utility::vector1< utility::vector1< float > > > temp( query_.size(), utility::vector1< utility::vector1< float > > ( 37, utility::vector1< float > ( 37, 0 ) ) );

	for( Size i = 2; i <= query_.size() -1; ++i ) {
		std::string curr_aa, next_aa, aa_type;//, ss_type;
		utility::vector1< core::Real > ss_weight( 3, 0 );
		utility::vector1< std::string > ss_types( 3, "" );

		ss_weight[1] = query_ss_->helix_fraction(i);
		ss_weight[2] = query_ss_->strand_fraction(i);
		ss_weight[3] = query_ss_->loop_fraction(i);

		if ( ( ss_weight[1] + ss_weight[2] + ss_weight[3]) == 0) {
			ss_weight[1] = 1.0;
			ss_weight[2] = 1.0;
			ss_weight[3] = 1.0;
		}

		ss_types[1] = "H";
		ss_types[2] = "E";
		ss_types[3] = "L";

		curr_aa = query_[i-1];
		next_aa = query_[i];

		//(C) is a non-disulfide cysteine, and (c) is a disulfide cysteine
		// but then the database needs to refer to them as "C" and "dc" because macs
		// ignore capitalization.
		if ( curr_aa == "c" ) {
			curr_aa = "dc";
		}

		//There are 42 amino acid types.
		// 21 for the normal set w. cysteins split into non-disulfide (C) and disulfide (c)
		// 21 for the above set but where the next residue is proline
		if ( next_aa == "P" ) {
			aa_type = curr_aa+next_aa;
		} else {
			aa_type = curr_aa;
		}

		//The ramachandran fragment score is a sigmoid function based on sequence and secondary structure specific phi/psi counts
		//However secondary structure is at this point just a weighted probability, so first we have to combine the
		//sequence specific H, E & L counts to create a secondary structure weighted ramachandran table.
		//(note: the counts here were extracted from the vall and have been gaussian smoothed to blur out the noise)
		for( Size s = 1; s <= 3; ++s ) {
			if ( ss_weight[s] > 0.0 ) {
				//std::string db_location("/work/rvernon/fragpicking_tests/vall/final/"+ss_types[s]+"_"+aa_type+".counts");
				std::string db_location("sampling/fragpicker_rama_tables/"+ss_types[s]+"_"+aa_type+".counts");
				utility::io::izstream table_file;
				basic::database::open(table_file, db_location);
				//table_file.open(db_location);

				std::string line;

				while ( getline(table_file, line) ) {
					if ( line.length() != 0 ) {
						std::istringstream line_stream(line);
						Size x, y;
						float count;

						line_stream >> x >> y >> count;
						temp[i][x+1][y+1] += ss_weight[s] * count;
						//std::cout << "HEYO " << i << " " << x << " " << y << " " << ss_weight[s] << " " << count << std::endl;
					}
				}
			}
		}

		utility::io::ozstream outtable;
		if (option[frags::write_rama_tables].user()) {
			std::string res = ObjexxFCL::string_of( i );
			outtable.open("res"+res+"_"+aa_type+".rama_table");
		}

		//Post-Normalization Multiplier (puts things back roughly into the magnitude of the raw counts
		//If zero then don't bother normalizing, just use raw counts
		float const A( option[frags::rama_norm] );
		if (A > 0.0) {
			Real total(0.0);
			for( Size x = 1; x <= 37; ++x ) {
				for( Size y = 1; y <= 37; ++y ) {
					total += temp[i][x][y];
				}
			}
			runtime_assert( total != 0.0 );
			for( Size x = 1; x <= 37; ++x ) {
				for( Size y = 1; y <= 37; ++y ) {
					temp[i][x][y] = (temp[i][x][y] / total) * A;
				}
			}
		}

		float const C( option[frags::rama_C] ); //default 0.0 <- Sigmoid inflection point adjuster
		float const B( option[frags::rama_B] ); //default 1.0 <- Sigmoid slope adjuster
		//Now we convert the count tables into sigmoid function score tables. The score goes from 1 (no counts) to
		//0 (many counts). Because this is a sigmoid there is a sharp transition between 1 and 0, this transition
		//takes place at an arbitrary point defined by me. It can be changed by adding in a constant to the final exp.
		for( Size x = 1; x <= 37; ++x ) {
			for( Size y = 1; y <= 37; ++y ) {
				temp[i][x][y] += 0.000000000000000000000000001;
				temp[i][x][y] = std::log(temp[i][x][y]);
				temp[i][x][y] = 1.0 / ( 1 + std::exp( C + B*temp[i][x][y] ) );

				if (option[frags::write_rama_tables].user()) {
					float xf( static_cast< float >( x ));
					float yf( static_cast< float >( y ));
					outtable << ((xf-1)*10)-175 << " " << ((yf-1)*10)-175 << " " << temp[i][x][y] << std::endl;
				}
			}
			//This blank line is so I can plot the tables in gnuplot. Don't judge me! -rv
			if (option[frags::write_rama_tables].user()) outtable << std::endl;
		}
	}

	sequence_rama_tables_ = temp;

}

void RamaScore::do_caching(VallChunkOP current_chunk) {

	std::string & tmp = current_chunk->chunk_key();
	if (tmp.compare(cached_scores_id_) == 0)
		return;
	cached_scores_id_ = tmp;

	Size query_sequence_length = query_.size();

	utility::vector1< utility::vector1< Real > > temp( current_chunk->size(),
	utility::vector1< Real > (query_sequence_length, 0 ) );

  runtime_assert( query_sequence_length > 0 );

	for (Size r = 2; r <= query_sequence_length - 1; ++r) {

		for (Size i = 1; i <= current_chunk->size(); ++i) {
			VallResidueOP res = current_chunk->at(i);

			Real phi = res->phi();
			Real psi = res->psi();

			//Frigging vall...
			if ( phi > 180 ) phi = -180 + (phi - 180);
			if ( psi > 180 ) psi = -180 + (psi - 180);
			if ( phi < -180 )	phi = 180 + (phi + 180);
			if ( psi < -180 )	psi = 180 + (psi + 180);

			Size i_phi = static_cast< Size > (((phi + 180)/10)+1);
			Size i_psi = static_cast< Size > (((psi + 180)/10)+1);

			runtime_assert( (i_phi >= 1) && (i_phi <= 37) );
			runtime_assert( (i_psi >= 1) && (i_psi <= 37) );

			temp[i][r] = sequence_rama_tables_[r][i_phi][i_psi];
		}
	}

	scores_ = temp;
}

void RamaScore::clean_up() {
}

bool RamaScore::score(FragmentCandidateOP fragment,FragmentScoreMapOP scores) {
	return cached_score( fragment, scores );
}

bool RamaScore::cached_score(FragmentCandidateOP fragment,
		FragmentScoreMapOP scores) {


	std::string & tmp = fragment->get_chunk()->chunk_key();
	if (tmp.compare(cached_scores_id_) != 0)
		do_caching(fragment->get_chunk());


	Real totalScore = 0.0;

	for (Size i = 1; i <= fragment->get_length(); i++) {
//		runtime_assert(fragment->get_first_index_in_vall()	+ i - 1 <= scores_.size());
//		runtime_assert(fragment->get_first_index_in_query() + i - 1 <= scores_[1].size());

		Real tmp = scores_[fragment->get_first_index_in_vall() + i - 1]
			                [fragment->get_first_index_in_query()	+ i - 1];
		totalScore += tmp;
	}

	//std::cout << "TOTALSCORE " << totalCount << " " << totalScore << " " << totalScore / (Real) fragment->get_length();
	totalScore /= (Real) fragment->get_length();

	scores->set_score_component(totalScore, id_);
	if ((totalScore < lowest_acceptable_value_) && (use_lowest_ == true))
		return false;
	return true;
}

bool RamaScore::describe_score(FragmentCandidateOP,
															 FragmentScoreMapOP, std::ostream&)
{
	return true;
}
															 //
															 //    return true;
	//}


	//bool RamaScore::describe_score(FragmentCandidateOP f,
															 //		FragmentScoreMapOP empty_map, std::ostream& out) {
															 //
															 //    return true;
	//}


utility::vector1< utility::vector1< utility::vector1< Real > > > RamaScore::sequence_rama_tables_;

} // scores
} // frag_picker
} // protocols


