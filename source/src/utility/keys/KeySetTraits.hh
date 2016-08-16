// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/KeySetTraits.hh
/// @brief  KeySet traits class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_keys_KeySetTraits_hh
#define INCLUDED_utility_keys_KeySetTraits_hh


// C++ headers
#include <set>


namespace utility {
namespace keys {


/// @brief KeySet traits class
template< typename K >
class KeySetTraits
{


private: // Types


	typedef  std::set< K >  Set;


public: // Types


	// STL/boost style
	typedef  K  key_type;
	typedef  typename Set::value_type  value_type;
	typedef  typename Set::reference  reference;
	typedef  typename Set::const_reference  const_reference;
	typedef  typename Set::pointer  pointer;
	typedef  typename Set::const_pointer  const_pointer;
	typedef  typename Set::iterator  iterator;
	typedef  typename Set::const_iterator  const_iterator;
	typedef  typename Set::reverse_iterator  reverse_iterator;
	typedef  typename Set::const_reverse_iterator  const_reverse_iterator;
	typedef  typename Set::size_type  size_type;
	typedef  typename Set::difference_type  difference_type;
	typedef  typename Set::allocator_type  allocator_type;

	// Project style
	typedef  K  Key;
	typedef  typename Set::value_type  Value;
	typedef  typename Set::reference  Reference;
	typedef  typename Set::const_reference  ConstReference;
	typedef  typename Set::pointer  Pointer;
	typedef  typename Set::const_pointer  ConstPointer;
	typedef  typename Set::iterator  Iterator;
	typedef  typename Set::const_iterator  ConstIterator;
	typedef  typename Set::reverse_iterator  ReverseIterator;
	typedef  typename Set::const_reverse_iterator  ConstReverseIterator;
	typedef  typename Set::size_type  Size;
	typedef  typename Set::difference_type  Difference;
	typedef  typename Set::allocator_type  Allocator;


}; // KeySetTraits


} // namespace keys
} // namespace utility


#endif // INCLUDED_utility_keys_KeySetTraits_HH
