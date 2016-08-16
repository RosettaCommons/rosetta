// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/keys/KeyVectorTraits.hh
/// @brief  KeyVector traits class
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)


#ifndef INCLUDED_utility_keys_KeyVectorTraits_hh
#define INCLUDED_utility_keys_KeyVectorTraits_hh


// Project headers
#include <utility/vector1.hh>


namespace utility {
namespace keys {


/// @brief KeyVector traits class
template< typename K >
class KeyVectorTraits
{


private: // Types


	typedef  vector1< K >  Vector;


public: // Types


	// STL/boost style
	typedef  K  key_type;
	typedef  typename Vector::value_type  value_type;
	typedef  typename Vector::reference  reference;
	typedef  typename Vector::const_reference  const_reference;
	typedef  typename Vector::pointer  pointer;
	typedef  typename Vector::const_pointer  const_pointer;
	typedef  typename Vector::iterator  iterator;
	typedef  typename Vector::const_iterator  const_iterator;
	typedef  typename Vector::reverse_iterator  reverse_iterator;
	typedef  typename Vector::const_reverse_iterator  const_reverse_iterator;
	typedef  typename Vector::size_type  size_type;
	typedef  typename Vector::index_type  index_type;
	typedef  typename Vector::difference_type  difference_type;
	typedef  typename Vector::allocator_type  allocator_type;

	// Project style
	typedef  K  Key;
	typedef  typename Vector::Value  Value;
	typedef  typename Vector::Reference  Reference;
	typedef  typename Vector::ConstReference  ConstReference;
	typedef  typename Vector::Pointer  Pointer;
	typedef  typename Vector::ConstPointer  ConstPointer;
	typedef  typename Vector::Iterator  Iterator;
	typedef  typename Vector::ConstIterator  ConstIterator;
	typedef  typename Vector::ReverseIterator  ReverseIterator;
	typedef  typename Vector::ConstReverseIterator  ConstReverseIterator;
	typedef  typename Vector::Size  Size;
	typedef  typename Vector::Index  Index;
	typedef  typename Vector::Difference  Difference;
	typedef  typename Vector::Allocator  Allocator;


}; // KeyVectorTraits


} // namespace keys
} // namespace utility


#endif // INCLUDED_utility_keys_KeyVectorTraits_HH
