// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/InvKinLigLoopDesign/std_extra.hh
///
/// @brief
/// @author

#ifndef DEVEL_INVKINLIGLOOPDESIGN_STD_EXTRA_HH
#define DEVEL_INVKINLIGLOOPDESIGN_STD_EXTRA_HH

#include <iostream>
#include <vector>
#include <map>
#include <utility/assert.hh>

namespace devel {

namespace inv_kin_lig_loop_design {

template< class K, class V >
V const& find_or_throw( std::map<K,V> const& m, K const& k ) {
	typename std::map<K,V>::const_iterator i = m.find(k);
	if ( i == m.end() ) {
		std::cout << "find_or_throw - couldn't find " << k << std::endl;
		assert( false );
	}
	return i->second;
}

template< class K, class V >
V const& find_or_default( std::map<K,V> const& m, K const& k, V const& v_default ) {
	typename std::map<K,V>::const_iterator i = m.find(k);
	if ( i != m.end() ) {
		return i->second;
	} else {
		return v_default;
	}
}

template<class T>
std::vector<T> operator+=(std::vector<T>& v, const std::vector<T>& other) {
	v.insert(v.end(),other.begin(),other.end());
	return v;
} // operator+=

} //

}

#endif // STD_EXTRA_HH
