// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/options/ResidueChainVectorOption.hh
/// @brief  Program integer vector option class
/// @author Modified by Rhiju Das (rhiju@stanford.edu)


// Unit headers
#include <utility/options/ResidueChainVectorOption.hh>
#include <utility/string_util.hh>

// ObjexxFCL headers
#include <ObjexxFCL/char.functions.hh>

namespace utility {
namespace options {

void
ResidueChainVectorOption::copy_from( Option const & other )
{
	debug_assert( dynamic_cast< ResidueChainVectorOption const * > ( & other ));

	ResidueChainVectorOption const & rcvect_opt_other =
		dynamic_cast< ResidueChainVectorOption const & > ( other );

	Super::operator = ( rcvect_opt_other ); // rely on parent's assignment operator
	value_strings_ = rcvect_opt_other.value_strings_;
}

// @brief Specialized function that converts tags like A:1-4 B:1-3 into a pair of  [1 2 3 4 1 2 3], [A A A A B B B]

VectorOption_T_< ResidueChainVectorOptionKey, int > &
ResidueChainVectorOption::cl_value( std::string const & value_str ){
	value_strings_.push_back( value_str ); // need this to reconstitute resnum_and_chain.
	return VectorOption_T_< ResidueChainVectorOptionKey, int >::cl_value( value_str );
}

// @brief Specialized function that converts tags like A:1-4 B:1-3 into a pair of  [1 2 3 4 1 2 3], [A A A A B B B]
std::pair< utility::vector1<int>, utility::vector1<char> >
ResidueChainVectorOption::resnum_and_chain() const {
	bool string_is_ok;
	utility::vector1<int> resnum;
	utility::vector1<char> chains;
	for ( Size n = 1; n <= value_strings_.size(); n++ ) {
		std::pair< std::vector< int >, std::vector< char > > const resnum_and_chain = get_resnum_and_chain( value_strings_[n], string_is_ok );
		runtime_assert( string_is_ok );
		for ( Size n = 0; n < resnum_and_chain.first.size(); n++ ) {
			resnum.push_back( resnum_and_chain.first[n] );
			chains.push_back( resnum_and_chain.second[n] );
		}
	}
	return std::make_pair( resnum, chains );
}

/// @brief Value of a string
int
ResidueChainVectorOption::value_of( std::string const & value_str ) const
{
	std::vector< int >  resnum;
	std::vector< char > chains;
	bool string_is_ok = get_resnum_and_chain_from_one_tag( value_str, resnum, chains );
	if ( ( value_str.empty() ) || ( ! string_is_ok ) || resnum.size() == 0 ) {
		std::cerr << "ERROR: Illegal value for resnum/chain option -" << id()
			<< " specified: " << value_str << std::endl;
		std::exit( EXIT_FAILURE );
	}

	return resnum[ 0 ];
}

/// @brief Value of a string
utility::vector1< int >
ResidueChainVectorOption::values_of( std::string const & value_str ) const
{
	bool string_is_ok( false );
	std::pair< std::vector< int >, std::vector< char > > const resnum_and_chain_info = get_resnum_and_chain( value_str, string_is_ok );
	if ( !string_is_ok ) {
		std::cerr << "ERROR: Illegal value for resnum/chain option -" << id()
			<< " specified: " << value_strings_ << std::endl;
		std::exit( EXIT_FAILURE );
	}

	//convert to utility vector1. This is a workaround to prevent ObjexxFCL from knowing about vector1.
	Values vector1_ints;
	std::vector< int > const & resnum = resnum_and_chain_info.first;
	for ( int n : resnum ) {
		vector1_ints.push_back( n );
	}
	return vector1_ints;
}

} // namespace options
} // namespace utility


