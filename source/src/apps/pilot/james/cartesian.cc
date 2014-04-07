// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author James Thompson

#include <iostream>
#include <set>
#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>

template< typename T >
utility::vector1< utility::vector1< T > >
cartesian_product(
	utility::vector1< std::set< T > > sets
) {
	using std::set;
	using utility::vector1;

	utility::vector1< utility::vector1< T > > product;

	typedef typename set< T >::const_iterator set_iter;
	typedef typename vector1< set< T > >::const_iterator iter;
	typedef typename vector1< vector1< T > >::const_iterator vec_iter;
	for ( iter it = sets.begin(), end = sets.end(); it != end; ++it ) {
		if ( it == sets.begin() ) {
			for ( set_iter s_it = it->begin(), s_end = it->end();
					s_it != s_end; ++s_it
			) {
				vector1< T > new_set( 1, *s_it );
				product.push_back( new_set );
			}
		} else {
			vector1< vector1< T > > new_product;
			for ( set_iter s_it = it->begin(), s_end = it->end();
					s_it != s_end; ++s_it
			) {
				for ( vec_iter v_it = product.begin(), v_end = product.end();
					v_it != v_end; ++v_it
				) {
					vector1< T > new_vector( *v_it );
					new_vector.push_back( *s_it );
					new_product.push_back( new_vector );
				}
			}
			product = new_product;
		}
	} // for sets

	return product;
} // cartesian_product

int
main( int /*argc*/, char* /*argv*/ [] ) {
	try {

	using std::set;
	using utility::vector1;
	set< char > set1;
	set< char > set2;

	set1.insert( 'A' );
	set1.insert( 'B' );
	set1.insert( 'C' );

	set2.insert( '1' );
	set2.insert( '2' );
	set2.insert( '3' );

	vector1< set< char > > sets;
	sets.push_back( set1 );
	sets.push_back( set2 );

	vector1< vector1< char > > product
		= cartesian_product< char > ( sets );

	typedef vector1< vector1< char > >::const_iterator iter;
	for ( iter it = product.begin(), end = product.end(); it != end; ++it ) {
		for ( vector1< char >::const_iterator it2 = it->begin(), end2 = it->end();
				it2 != end2; ++it2
		) {
			std::cout << *it2;
		}
		std::cout << std::endl;
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main( int argc, char * argv [] )
