// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farna/secstruct/RNA_SecStruct.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/secstruct/RNA_SecStruct.hh>
#include <core/sequence/util.hh>
#include <utility/io/izstream.hh>
#include <utility/tools/make_vector1.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.farna.secstruct.RNA_SecStruct" );

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
//////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace farna {
namespace secstruct {


//Constructor
RNA_SecStruct::RNA_SecStruct( std::string const & secstruct,
	std::string const & secstruct_file /* = "" */,
	std::string const & sequence /* = "" */ )
{
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
	check_balanced_secstruct();
}

//Destructor
RNA_SecStruct::~RNA_SecStruct()
{}

//////////////////////////////////////////////////////////
void
RNA_SecStruct::set_secstruct( std::string const & secstruct )
{
	secstruct_ = secstruct;
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
	runtime_assert( std::count( secstruct_.begin(), secstruct_.end(),'(') ==
		std::count( secstruct_.begin(), secstruct_.end(),')') );
	runtime_assert( std::count( secstruct_.begin(), secstruct_.end(),'[') ==
		std::count( secstruct_.begin(), secstruct_.end(),']') );
	runtime_assert( std::count( secstruct_.begin(), secstruct_.end(),'{') ==
		std::count( secstruct_.begin(), secstruct_.end(),'}') );
	return true;
}

//////////////////////////////////////////////////////////
void
RNA_SecStruct::check_compatible_with_sequence( std::string const & sequence_in,
	bool const check_complementarity /* = true */ ) const
{
	std::string sequence = sequence_in;
	vector1< Size > spacer_positions_sequence = core::sequence::strip_spacers( sequence );
	if ( sequence.size() != secstruct_.size() ) {
		TR << "'" << sequence  << "'" << std::endl;;
		TR << "'" << secstruct_ << "'" << std::endl;
		utility_exit_with_message( "Length of sequence & secstruct do not match: " +  ObjexxFCL::format::I(5,sequence.size()) + " " + ObjexxFCL::format::I(5,secstruct_.size()) );
	}

	/////////////////////////////////////////////////
	/////////////////////////////////////////////////
	// check that spacer positions don't disagree.
	/////////////////////////////////////////////////
	/////////////////////////////////////////////////

	if ( check_complementarity ) {
		/////////////////////////////////////////////////
		/////////////////////////////////////////////////
		// check for complementarity at matched bases:
		/////////////////////////////////////////////////
		/////////////////////////////////////////////////
	}
}

//////////////////////////////////////////////////////////
void
RNA_SecStruct::get_stems(
	vector1< vector1< std::pair< Size, Size > > > & stems,
	std::string const & line /*secstruct*/,
	std::string const & sequence_for_fasta,
	vector1< Size > const & chainbreak_pos,
	char const & left_bracket_char, char const & right_bracket_char ) const
{
	Size count( 0 );
	vector1< Size > left_brackets;
	std::map< Size, Size > pair_map;
	vector1< std::pair< Size, Size > > all_pairs;

	for ( Size i = 0; i < line.size(); i++ ) {
		if ( core::sequence::spacers.has_value( line[i] ) ) continue;
		count++;
		if ( line[i] == left_bracket_char ) left_brackets.push_back( count );
		if ( line[i] == right_bracket_char ) {
			if ( left_brackets.size() == 0 ) {
				utility_exit_with_message( "Number of right brackets does not match left brackets" );
			}
			Size const res1 = left_brackets[ left_brackets.size() ];
			Size const res2 = count;
			left_brackets.remove_back();
			pair_map[ res1 ] = res2;
			pair_map[ res2 ] = res1;
			all_pairs.push_back( std::make_pair( res1, res2 ) );
			if ( sequence_for_fasta.size() > 0 ) {
				std::map< char, utility::vector1< char > >::const_iterator it =  rna_complement_.find( sequence_for_fasta[res2-1] );
				if ( it == rna_complement_.end() ||
						( !( it->second ).has_value( sequence_for_fasta[res1-1] ) ) ) {
					utility_exit_with_message( "Not complementary at positions " + A(1,sequence_for_fasta[res1-1]) + I(5,res1) + " and " + A(1,sequence_for_fasta[res2-1]) + I(5,res2) );
				}
			}
		}
	}
	if ( left_brackets.size() > 0 ) {
		utility_exit_with_message( "Number of right brackets does not match left brackets" );
	}
	Size numres( count );

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
					!chainbreak_pos.has_value( k ) &&
					!chainbreak_pos.has_value( pair_map[ k + 1 ] ) ) {
				k++;
				stem_res.push_back( std::make_pair( k, pair_map[k] ) );
				already_in_stem[ k ] = true;
				already_in_stem[ pair_map[k] ] = true;
			}
			stems.push_back( stem_res );
		}
	}
}

///////////////////////////////////////////////////////////////////
// @details figures out stems based on (),[],{} bracket markers.
//   accepts cutpoint_open (chainbreaks) to help decide on stem boundaries
//   in principle could/should use spacer_positions_ to do this.
vector1< vector1< std::pair< Size, Size > > >
RNA_SecStruct::get_all_stems( std::string const & sequence /* = "" */,
	vector1< Size > const & cutpoint_open /* = vector1< Size >() */ ) const
{
	vector1< vector1< std::pair< Size, Size > > > stems;
	get_stems( stems, secstruct_, sequence, cutpoint_open, '(', ')' );
	get_stems( stems, secstruct_, sequence, cutpoint_open, '[', ']' );
	get_stems( stems, secstruct_, sequence, cutpoint_open, '{', '}' );
	return stems;
}

////////////////////////////////////////////////////////////////
void
RNA_SecStruct::remove_pair( std::pair< Size, Size > pair ) {
	Size pos1 = pair.first;
	Size pos2 = pair.second;
	secstruct_[ pos1 - 1 ] = '.';
	secstruct_[ pos2 - 1 ] = '.';
}

////////////////////////////////////////////////////////////////
bool
RNA_SecStruct::blank() const {
	for ( Size n = 1; n <= secstruct_.size(); n++ ) {
		if ( secstruct_[ n-1 ] != '.' ) return false;
	}
	return true;
}

} //secstruct
} //farna
} //protocols
