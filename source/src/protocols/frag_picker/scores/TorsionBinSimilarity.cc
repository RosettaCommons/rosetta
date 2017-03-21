// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/TorsionBinSimilarity.cc
/// @brief
/// @author James Thompson

#include <protocols/frag_picker/scores/TorsionBinSimilarity.hh>

#include <core/types.hh>
#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/TorsionBinIO.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

namespace protocols {
namespace frag_picker {
namespace scores {

static THREAD_LOCAL basic::Tracer tr(
	"protocols.frag_picker.scores.TorsionBinSimilarity"
);

void TorsionBinSimilarity::do_caching( VallChunkOP chunk ) {

	std::string tmp = chunk->chunk_key();
	if ( tmp.compare(cached_scores_id_) == 0 ) {
		return;
	}
	cached_scores_id_ = tmp;

	tr.Debug << "caching score for " << chunk->get_pdb_id()
		<< " of size " << chunk->size() << std::endl;

	for ( core::Size i = 1; i <= query_len_; ++i ) {
		for ( core::Size j = 1; j <= chunk->size(); ++j ) {
			core::Real const phi  ( chunk->at(j)->phi()   );
			core::Real const psi  ( chunk->at(j)->psi()   );
			core::Real const omega( chunk->at(j)->omega() );
			char const torsion_bin( torsion2big_bin_( phi, psi, omega ) );
			//scores_[i][j] = std::log( query_bin_probs_[i][bin_index_(torsion_bin)] );
			scores_[i][j] = 1 - query_bin_probs_[i][bin_index_(torsion_bin)];
		}
	}

	tr.Debug << "precomputed matrix of scores " << scores_.size()
		<< "x" << chunk->size() << std::endl;
}

bool TorsionBinSimilarity::cached_score(
	FragmentCandidateOP f,
	FragmentScoreMapOP empty_map
) {


	std::string tmp = f->get_chunk()->chunk_key();
	if ( tmp.compare(cached_scores_id_) != 0 ) {
		do_caching(f->get_chunk());
	}

	core::Real total_score = 0;
	for ( core::Size i = 1; i <= f->get_length(); ++i ) {
		debug_assert( f->get_first_index_in_query() + i - 1 <= scores_.size()    );
		debug_assert( f->get_first_index_in_vall() + i - 1 <= scores_[1].size() );
		total_score += scores_[f->get_first_index_in_query() + i - 1][f->get_first_index_in_vall() + i - 1];
	}
	//std::cout << "total_score = " << total_score << std::endl;

	total_score /= (core::Real) f->get_length();

	empty_map->set_score_component(total_score, id_);
	if ( (total_score > lowest_acceptable_value_) && (use_lowest_ == true) ) {
		return false;
	}
	return true;
}

char
TorsionBinSimilarity::torsion2big_bin_(
	core::Real const phi,
	core::Real const psi,
	core::Real const omega
) const {
	if ( std::abs( omega ) < 90 ) {
		return 'O'; // cis-omega
	} else if ( phi >= 0.0 ) {
		if ( -100 < psi && psi <= 100 ) {
			return 'G'; // alpha-L
		} else {
			return 'E'; // E
		}
	} else {
		if ( -125 < psi && psi <= 50 ) {
			return 'A'; // helical
		} else {
			return 'B'; // extended
		}
	}
	return 'X';
}

core::Size
TorsionBinSimilarity::bin_index_( char const bin_name ) const
{
	switch( bin_name ) {
	case 'A' : return 1;
	case 'B' : return 2;
	case 'E' : return 3;
	case 'G' : return 4;
	case 'O' : return 5;
	default :
		std::string const msg( "Error: don't recognize bin" + std::string(1, bin_name) );
		utility_exit_with_message(msg);
		break;
	}

	return 'X'; // satisfy compiler
}

FragmentScoringMethodOP MakeTorsionBinSimilarity::make(
	core::Size priority,
	core::Real lowest_acceptable_value,
	bool use_lowest,
	FragmentPickerOP picker,
	std::string /* prediction_id */
) {
	using core::Size;
	using core::Real;
	using utility::vector1;
	using utility::io::izstream;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::frag_picker::TorsionBinIO reader;
	if ( !option[ in::file::torsion_bin_probs ].user() ) {
		utility_exit_with_message("Error: no file specified");
	}
	reader.read( izstream( option[ in::file::torsion_bin_probs ]() ) );
	if ( reader.nrows() == 0 ) {
		utility_exit_with_message( "Error: didn't read any torsions!" );
	}
	core::Size const sequence_length( reader.nrows() );
	core::Size const vall_max_len( picker->get_vall()->get_largest_chunk_size() );

	utility::vector1< utility::vector1< core::Real > > probs = reader.matrix();

	return (
		FragmentScoringMethodOP( new TorsionBinSimilarity(
		priority,
		lowest_acceptable_value,
		use_lowest,
		probs,
		sequence_length,
		vall_max_len
		) )
	);
}

} // scores
} // frag_picker
} // protocols


