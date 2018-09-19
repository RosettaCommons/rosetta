// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/sequence/Blosum62Map.hh
/// @brief Creates a std::unordered_map that maps a char pair of name1s to their BLOSUM62 score
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_core_sequence_Blosum62Map_hh
#define INCLUDED_core_sequence_Blosum62Map_hh

#include <unordered_map>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <utility>

#ifndef NDEBUG
#include <basic/Tracer.hh>
static basic::Tracer blosum_62_map_tr( "core.sequence.Blosum62Map" );
#endif

namespace core {
namespace sequence {

struct CharPairHash {
public:
	std::size_t operator()( std::pair< char, char > const & x ) const {
		//simple bit-pack
		std::size_t first = x.first;
		first <<= 8;
		std::size_t second = x.second;
		return first | second;
	}
};

std::unordered_map< std::pair< char, char >, int, CharPairHash > create_map_for_Blosum62Map();


struct Blosum62Map {
public:
	Blosum62Map() :
		map_( create_map_for_Blosum62Map() )
	{}

	~Blosum62Map() = default;

	int score_for_aa_pair( std::pair< char, char > name1s ) const {
		debug_assert( map_.find( name1s ) != map_.end() );
		return map_.at( name1s );
	}

	int score_for_aa_pair( char aa1, char aa2 ) const {
		return score_for_aa_pair( std::make_pair( aa1, aa2 ) );
	}

private:
	std::unordered_map< std::pair< char, char >, int, CharPairHash > map_;
};

inline
std::unordered_map< std::pair< char, char >, int, CharPairHash >
create_map_for_Blosum62Map(){

	std::unordered_map< std::pair< char, char >, int, CharPairHash > score_for_aa_pair;

	utility::vector1< char > aas = { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X' };

	//Feel free to call me a lazy programmer, but this is just so much easier than reading from the database.
	utility::vector1< utility::vector1< int > > values = {
		// A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X
		{  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0 },//A
		{ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1 },//R
		{ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1 },//N
		{ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1 },//D
		{  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2 },//C
		{ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1 },//Q
		{ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1 },//E
		{  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1 },//G
		{ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1 },//H
		{ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1 },//I
		{ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1 },//L
		{ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1 },//K
		{ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1 },//M
		{ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1 },//F
		{ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2 },//P
		{  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0 },//S
		{  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0 },//T
		{ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2 },//W
		{ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1 },//Y
		{  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1 },//V
		{ -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1 },//B
		{ -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1 },//Z
		{  0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1 },//X
		};

	for ( core::Size i = 1; i <= aas.size(); ++i ) {
		for ( core::Size j = 1; j <= aas.size(); ++j ) {
			std::pair< char, char > char_pair = { aas[ i ], aas[ j ] };
			int value = values[ i ][ j ];
			std::pair< std::pair< char, char >, int > map_element = std::make_pair( char_pair, value );
			score_for_aa_pair.insert( map_element );
		}
	}

#ifndef NDEBUG
	//Count number of collisions in map
	core::Size max_bin_size = 0;
	core::Size num_occupied_bins = 0;
	for ( core::Size bin = 0; bin < score_for_aa_pair.bucket_count(); ++bin ) {
		auto const size = score_for_aa_pair.bucket_size( bin );
		if ( size > 0 ) {
			++num_occupied_bins;
			if ( size > max_bin_size ) max_bin_size = size;
		}
	}
	blosum_62_map_tr << score_for_aa_pair.size() << " elements stored in " << num_occupied_bins << " bins with a maximum bin size of " << max_bin_size << std::endl;
#endif

	return score_for_aa_pair;
}

} // sequence
} // core

#endif
