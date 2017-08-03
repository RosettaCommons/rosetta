// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/pose/rna/RNA_SecStruct.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/rna/RNA_SecStruct.hh>
#include <core/sequence/util.hh>
#include <utility/io/izstream.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/vector1.hh>
#include <utility/stream_util.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>
#include <set>

static basic::Tracer TR( "core.pose.rna.RNA_SecStruct" );

using namespace core;
using utility::tools::make_vector1;
using utility::vector1;
using ObjexxFCL::format::I;
using ObjexxFCL::format::A;

//////////////////////////////////////////////////////////////////////////////////////////
// @details
//
//  RNA_SecStruct
//
//  - originally written to hold a text representation (based on dot-bracket notation)
//     of RNA secondary structure. E.g.,
//
//         (((....)))  is a hairpin
//         (((.[[.)))]] is a pseudoknot
//
// TODO:
//  - need to upgrade to use vector of base pairs as primary representation, not a string
//  - need to upgrade to handle brackets beyond (,[,{. E.g., identify aaa...aaa as stems
//  - cache as RNA_SecStructInfo in pose
//  - actually make use of spacer_positions_ instead of input cutpoint_open in get_all_stems()
//  - set up cool use case where we use sec_struct information during FARNA/stepwise runs.
//  - 'promote' out of protocols:farna into core:pose:rna
//
//    Originally written by Rhiju Das, June 2016.
//
// AMW: did all of the above except the cool use case + caching (since nothing actually
// uses the currently-Pose-cached RNA_SecStructLegacyInfo string).
//
//////////////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace pose {
namespace rna {


//Constructor
RNA_SecStruct::RNA_SecStruct( std::string const & secstruct,
	std::string const & secstruct_file /* = "" */,
	std::string const & sequence /* = "" */ )
{
	// AMW todo: expand greatly, map on string or add NC map on string.
	rna_complement_.clear();
	rna_complement_[ 'a' ] = make_vector1( 'u', 't' );
	rna_complement_[ 'u' ] = make_vector1( 'a','g' );
	rna_complement_[ 't' ] = make_vector1( 'a','g' );
	rna_complement_[ 'c' ] = make_vector1( 'g' );
	rna_complement_[ 'g' ] = make_vector1( 'c','u' );

	set_secstruct( secstruct );
	if ( secstruct_file.size() > 0 ) {
		runtime_assert( secstruct_.size() == 0 ); // can only supply secstruct_file or secstruct, not both.
		read_secstruct_from_file( secstruct_file );
	}

	if ( secstruct_.size() == 0 ) blank_secstruct( sequence );
	// We don't want to depend on the sequence here: could be a
	// 'general' secstruct.
	set_basepairs_from_secstruct();

	check_balanced_secstruct();

	figure_out_stems();
}

//Destructor
RNA_SecStruct::~RNA_SecStruct()
{}

void
RNA_SecStruct::set_basepairs_from_secstruct() {

	// We don't care which BPs are in stems with each other yet
	//vector1< vector1< std::pair< Size, Size > > > & stems,
	//std::string const & line /*secstruct*/,
	//std::string const & sequence_for_fasta,
	//vector1< Size > const & chainbreak_pos,
	//char const & left_bracket_char, char const & right_bracket_char ) const

	// Go through pairs of left and right bracket chars.
	// After those base pairs are assigned, find matching e.g. aaaaa....aaaaa
	// Note that some issues may come up if we do not carefully force a 'break'
	// of ... or a chainbreak character.

	std::string const lefts  = "([{<";
	std::string const rights = ")]}>";

	for ( Size ii = 0; ii < lefts.size(); ++ii ) {
		char left_bracket_char = lefts[ii], right_bracket_char = rights[ii];

		Size count( 0 );
		vector1< Size > left_brackets;

		for ( Size i = 0; i < secstruct_.size(); i++ ) {
			if ( core::sequence::spacers.has_value( secstruct_[i] ) ) continue;
			count++;
			if ( secstruct_[i] == left_bracket_char ) left_brackets.push_back( count );
			if ( secstruct_[i] == right_bracket_char ) {
				if ( left_brackets.size() == 0 ) {
					utility_exit_with_message( "Number of right brackets does not match left brackets" );
				}
				Size const res1 = left_brackets[ left_brackets.size() ];
				Size const res2 = count;
				left_brackets.remove_back();
				base_pairs_.emplace_back( res1, res2 );
			}
		}
		if ( left_brackets.size() > 0 ) {
			utility_exit_with_message( "Number of right brackets does not match left brackets" );
		}
	}

	// UNICHAR
	// We need to be able to assign new base pairs, without
	// creating unrealistic stems. That is,
	// ...aa.a....aaa
	// needs to act like
	// ...((.(....)))
	// not
	// ...((.)....)((
	// which is a danger if we get greedy and try to close pairs ASAP once
	// we hit a gap of N residues or more.
	// So, we have to have a rule. These non-paired characters get to be used once
	// per SS string -- we have enough of them! -- so we can say that the first
	// and final appearance of each one matches, and so on inward.

	// AMW: at some point, generate this string sanely
	std::string const single_chars = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

	for ( char const single_char : single_chars ) {

		// Obtain set of all occurrences of each character.
		utility::vector1< Size > occurrences;
		for ( Size i = 0; i < secstruct_.size(); i++ ) {
			if ( secstruct_[i] == single_char ) {
				occurrences.push_back( i + 1 );
			}
		}
		if ( occurrences.empty() ) continue;

		// Iterate inward, pairing first with last, second with penultimate, etc.
		//auto bp1 = occurrences.begin(), bp2 = occurrences.end();
		//--bp2;
		for ( Size ii = 1; ii < occurrences.size() - ii; ++ii ) {
			//while ( std::distance( bp1, bp2 ) > 0 ) {
			base_pairs_.emplace_back( occurrences[ ii ], occurrences[ occurrences.size() - ii + 1 ] );
			//++bp1; --bp2;
		}
	}
}

//////////////////////////////////////////////////////////
void
RNA_SecStruct::set_secstruct( std::string const & secstruct )
{
	// If the final character is a space, it DOES NOT count. Strip it off, I say.
	// I'm not sure how this state is obtained now but couldn't have been earlier,
	// but there you go.
	secstruct_ = secstruct;
	if ( secstruct_[ secstruct.size() - 1 ] == ' ' ) {
		secstruct_ = secstruct_.substr( 0, secstruct.size() - 1 );
	}
	spacer_positions_ = core::sequence::strip_spacers( secstruct_, false /*annotations_in_brackets*/ );
}

//////////////////////////////////////////////////////////////////////////////////////
// @details Pretty naive at the moment -- just reads in first line from text file.
void
RNA_SecStruct::read_secstruct_from_file( std::string const & filename )
{
	TR << "Reading RNA secstruct file: " << filename << std::endl;

	utility::io::izstream data_stream( filename );
	if ( !data_stream ) {
		data_stream.close();
		utility_exit_with_message(  "Data file? " + filename );
	}

	std::string line, secstruct;
	getline( data_stream, line );

	std::istringstream line_stream( line );
	while ( !line_stream.fail() ) {
		std::string secstruct_string;
		line_stream >> secstruct_string;
		if ( secstruct.size() > 0 ) secstruct += " " ;
		secstruct += secstruct_string;
	}
	data_stream.close();

	set_secstruct( secstruct ); // handles spaces, etc.
}

//////////////////////////////////////////////////////////
void
RNA_SecStruct::blank_secstruct( std::string const & sequence_in )
{
	std::string sequence = sequence_in;
	spacer_positions_ = core::sequence::strip_spacers( sequence );
	secstruct_ = "";
	for ( Size n = 1; n <= sequence.size(); n++ ) secstruct_ += '.';
}

//////////////////////////////////////////////////////////
bool
RNA_SecStruct::check_balanced_secstruct() const
{
	// AMW: not really necessary now that the base pairs are created/checked at
	// construction
	runtime_assert( std::count( secstruct_.begin(), secstruct_.end(),'(') ==
		std::count( secstruct_.begin(), secstruct_.end(),')') );
	runtime_assert( std::count( secstruct_.begin(), secstruct_.end(),'[') ==
		std::count( secstruct_.begin(), secstruct_.end(),']') );
	runtime_assert( std::count( secstruct_.begin(), secstruct_.end(),'{') ==
		std::count( secstruct_.begin(), secstruct_.end(),'}') );
	runtime_assert( std::count( secstruct_.begin(), secstruct_.end(),'<') ==
		std::count( secstruct_.begin(), secstruct_.end(),'>') );
	return true;
}

//////////////////////////////////////////////////////////
void
RNA_SecStruct::check_compatible_with_sequence( std::string const & sequence_in,
	bool const check_complementarity /* = true */ ) const
{
	std::string sequence = sequence_in;
	// This is a no-op for any sequences read in from fasta or that have gone
	// through denovo setup; the spacers have already been stripped.
	// Therefore, don't throw a fit later when this is an empty vector.
	vector1< Size > spacer_positions_sequence = core::sequence::strip_spacers( sequence );

	if ( sequence.size() != secstruct_.size() ) {
		TR << "'" << sequence  << "'" << std::endl;;
		TR << "'" << secstruct_ << "'" << std::endl;
		utility_exit_with_message( "Length of sequence & secstruct do not match: " +  ObjexxFCL::format::I(5,sequence.size()) + " " + ObjexxFCL::format::I(5,secstruct_.size()) );
	}

	//if ( spacer_positions_sequence != spacer_positions_ ) {
	// TR << "Spacer positions from sequence:  " << spacer_positions_sequence << std::endl;
	// TR << "Spacer positions from secstruct: " << spacer_positions_ << std::endl;
	//}

	if ( check_complementarity ) {
		for ( auto const & pair : base_pairs_ ) {
			Size const pos1 = pair.first;
			Size const pos2 = pair.second;

			auto it =  rna_complement_.find( sequence[pos2-1] );
			if ( it == rna_complement_.end() ||
					( !( it->second ).has_value( sequence[pos1-1] ) ) ) {
				utility_exit_with_message( "Not complementary at positions " + A(1,sequence[pos1-1]) + I(5,pos1) + " and " + A(1,sequence[pos2-1]) + I(5,pos2) );
			}

		}
	}
}

///////////////////////////////////////////////////////////////////
// @details figures out stems based on base pairs.
void
RNA_SecStruct::figure_out_stems()
{
	vector1< vector1< std::pair< Size, Size > > > stems;

	Size numres = 0;
	std::map< Size, Size > pair_map;
	for ( auto const & elem : base_pairs_ ) {
		if ( elem.first > numres  ) numres = elem.first;
		if ( elem.second > numres ) numres = elem.second;
		pair_map[ elem.first ] = elem.second;
		pair_map[ elem.second ] = elem.first;
	}

	// Parse out stems
	vector1< bool > already_in_stem( numres, false );

	Size stem_count( 0 );
	for ( Size i = 1; i <= numres; i++ ) {
		if ( pair_map.find( i ) != pair_map.end() && !already_in_stem[ i ] ) {
			// In a base pair
			Size k = i;
			stem_count += 1;
			vector1< std::pair< Size, Size > > stem_res;

			stem_res.push_back( std::make_pair( k, pair_map[k] ) );
			already_in_stem[ k ] = true;
			already_in_stem[ pair_map[k] ] = true;

			// Can we extend in one direction?
			while ( pair_map.find( k + 1 ) != pair_map.end() &&
					pair_map[ k + 1 ] == pair_map[ k ] - 1  &&
					!already_in_stem[ k + 1 ] &&
					!spacer_positions_.has_value( k ) &&
					!spacer_positions_.has_value( pair_map[ k + 1 ] ) ) {
				k++;
				stem_res.push_back( std::make_pair( k, pair_map[k] ) );
				already_in_stem[ k ] = true;
				already_in_stem[ pair_map[k] ] = true;
			}
			stems.push_back( stem_res );
		}
	}

	stems_ = stems;
}

////////////////////////////////////////////////////////////////
void
RNA_SecStruct::remove_pair( std::pair< Size, Size > pair ) {
	Size pos1 = pair.first;
	Size pos2 = pair.second;
	secstruct_[ pos1 - 1 ] = '.';
	secstruct_[ pos2 - 1 ] = '.';

	// Also remove base pairs from base_pairs_
	// AMW: conceivably one could trigger the other, or we could just stop
	// maintaining secstruct_ after construction.
	// But, it's nicer to use the same string representation you were built with
	// if ever you are displayed.

	// Standard erase-remove idiom
	base_pairs_.erase( std::remove( base_pairs_.begin(), base_pairs_.end(), pair ), base_pairs_.end());
}

////////////////////////////////////////////////////////////////
bool
RNA_SecStruct::blank() const {
	for ( Size n = 1; n <= secstruct_.size(); n++ ) {
		if ( secstruct_[ n-1 ] != '.' ) return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////
bool
RNA_SecStruct::in_helix( core::Size const & i ) const
{
	for ( auto const & pair : base_pairs_ ) {
		if ( i == pair.first ) return true;
		if ( i == pair.second ) return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////
bool
RNA_SecStruct::in_same_helix( core::Size const & i, core::Size const & j ) const
{
	for ( auto const & stem : stems_ ) {
		bool found_i( false ), found_j( false );
		for ( auto const & pair: stem ) {
			if ( i == pair.first || i == pair.second ) found_i = true;
			if ( j == pair.first || j == pair.second ) found_j = true;
		}
		if ( found_i && found_j ) return true;
	}
	return false;
}

} //secstruct
} //farna
} //protocols
