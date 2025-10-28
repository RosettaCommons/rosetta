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

#include <iostream>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace utility {
namespace options {

void
ResidueChainVectorOption::copy_from( Option const & other )
{
	debug_assert( dynamic_cast< ResidueChainVectorOption const * > ( & other ));

	auto const & rcvect_opt_other =
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
std::tuple< utility::vector1<int>, utility::vector1<std::string>, utility::vector1< std::string > >
ResidueChainVectorOption::resnum_and_chain() const {
	bool string_is_ok;
	utility::vector1<int> resnum;
	utility::vector1<std::string> chains;
	utility::vector1<std::string> segids;
	for ( Size n = 1; n <= value_strings_.size(); n++ ) {
		std::tuple< utility::vector1< int >, utility::vector1< std::string >, utility::vector1< std::string > > const resnum_and_chain = get_resnum_and_chain_and_segid( value_strings_[n], string_is_ok );
		runtime_assert( string_is_ok );
		for ( Size n2 = 1; n2 <= std::get< 0 >( resnum_and_chain ).size(); n2++ ) {
			resnum.push_back( std::get<0>(resnum_and_chain)[n2] );
			chains.push_back( std::get<1>(resnum_and_chain)[n2] );
			segids.push_back( std::get<2>(resnum_and_chain)[n2] );
		}
	}
	return std::make_tuple( resnum, chains, segids );
}

/// @brief Value of a string
int
ResidueChainVectorOption::value_of( std::string const & value_str ) const
{
	utility::vector1< int >  resnum;
	utility::vector1< std::string > chains;
	utility::vector1< std::string > segids;
	bool string_is_ok = get_resnum_and_chain_from_one_tag( value_str, resnum, chains, segids );
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
	std::tuple< utility::vector1< int >, utility::vector1< std::string >, utility::vector1< std::string > > const resnum_and_chain_info = get_resnum_and_chain_and_segid( value_str, string_is_ok );
	if ( !string_is_ok ) {
		std::cerr << "ERROR: Illegal value for resnum/chain option -" << id()
			<< " specified: {";
		for ( Size i=1, imax=value_strings_.size(); i<=imax; ++i ) {
			std::cerr << value_strings_[i];
			if ( i < imax ) std::cerr << ", ";
		}
		std::cerr << std::endl;
		std::exit( EXIT_FAILURE );
	}

	return std::get<0>(resnum_and_chain_info);
}

} // namespace options
} // namespace utility



#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
utility::options::ResidueChainVectorOption::save( Archive & arc ) const {
	arc( cereal::base_class< utility::options::VectorOption_T_<class utility::options::ResidueChainVectorOptionKey, int> >( this ) );
	arc( CEREAL_NVP( value_strings_ ) ); // utility::vector1<std::string>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
utility::options::ResidueChainVectorOption::load( Archive & arc ) {
	arc( cereal::base_class< utility::options::VectorOption_T_<class utility::options::ResidueChainVectorOptionKey, int> >( this ) );
	arc( value_strings_ ); // utility::vector1<std::string>
}

SAVE_AND_LOAD_SERIALIZABLE( utility::options::ResidueChainVectorOption );
CEREAL_REGISTER_TYPE( utility::options::ResidueChainVectorOption )

CEREAL_REGISTER_DYNAMIC_INIT( utility_options_ResidueChainVectorOption )
#endif // SERIALIZATION
