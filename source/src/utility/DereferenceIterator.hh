// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/DereferenceIterator.hh
/// @brief  class for Iterating across a container of OPs if the function expects iterators of the dereferenced objects.
/// example: if you have a loop that iterates on decoys via SilentFileData::begin()... SilentFileData::end() it uses
/// the dereferenced objects: i.e., it->fill_pose()
/// SilentFileData sfd;
/// my_func( sfd.begin(), sfd.end() );

/// if you have a container of OPs: i.e., typedef std::list< SilentStructOP > SilentStructs;
/// you can now call this function with the SilentStructs container instead of the SilentFileData objects by
///    my_func(    DereferenceIterator< SilentStructs >( decoys().begin() ),
///                   DereferenceIterator< SilentStructs >( decoys().end() ) );


/// @author  Oliver Lange


#ifndef INCLUDED_utility_DereferenceIterator_hh
#define INCLUDED_utility_DereferenceIterator_hh

// Unit headers
//#include <utility/DereferenceIterator.fwd.hh>
#include <utility/vector1.hh>

namespace utility {
/// @brief const_iterator class for SilentFileData container.
template< typename Container >
class DereferenceIterator {
	typedef typename Container::value_type value_type;
	typedef typename Container::const_iterator const_iterator;
	typedef value_type* pointer;
	typedef value_type& reference;
	typedef std::ptrdiff_t       difference_type;
	typedef std::bidirectional_iterator_tag iterator_category;

public:
	/// @brief empty constructor
	DereferenceIterator() {}

	/// @brief Constructor, given an iterator into the Structure_Map.
	DereferenceIterator( const_iterator s_iter ) {
		it_ = s_iter;
	}

	~DereferenceIterator() {}

	DereferenceIterator& operator=( const DereferenceIterator& src ) {
		it_ = src.it_;
		return (*this);
	}

	bool operator==( const DereferenceIterator& other ) const {
		return ( it_ == other.it_ );
	}

	bool operator!=( const DereferenceIterator& other ) const {
		return ( it_ != other.it_ );
	}

	DereferenceIterator& operator++() {
		it_++;
		return (*this);
	}

	DereferenceIterator& operator--() {
		it_--;
		return (*this);
	}

	value_type operator->() const {
		return *it_;
	}

	value_type operator*() const {
		return *it_;
	}

private:
	const_iterator it_; // keeps track of my place in a Structure_Map
}; // class iterator

}

#endif
