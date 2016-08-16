// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file   core/fragment/picking_old/concepts/Book.hh
/// @brief  wrapper class demonstrating the Book concept
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_concepts_Book_hh
#define INCLUDED_core_fragment_picking_old_concepts_Book_hh

// unit headers
#include <core/fragment/picking_old/concepts/Book.fwd.hh>

// type headers
#include <core/types.hh>

#include <core/fragment/picking_old/vall/VallResidue.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace concepts {


/// @brief class demonstrating the Book concept
/// @remarks Class is usable as a concrete implementation of a Book.
template< typename Pages >
class Book {


public: // typedefs


	typedef core::Size Size;
	typedef typename Pages::Page Page;
	typedef typename Pages::const_iterator PageConstIterator;
	typedef typename Pages::iterator PageIterator;


public: // construct/destruct


	/// @brief default constructor
	Book() {}


	/// @brief Pages constructor
	Book( Pages const & pages ) : pages_( pages ) {}


	/// @brief default destructor
	~Book() {}


public: // iterators


	/// @brief return an iterator that points to the first Page in this book
	inline
	PageConstIterator begin() const {
		return pages_.begin();
	}


	/// @brief return an iterator that points to the first Page in this book
	inline
	PageIterator begin() {
		return pages_.begin();
	}


	/// @brief return an iterator that points just beyond the last Page of this book
	inline
	PageConstIterator end() const {
		return pages_.end();
	}


	/// @brief return an iterator that points just beyond the last Page of this book
	inline
	PageIterator end() {
		return pages_.end();
	}


public: // page management


	/// @brief number of pages in the book
	inline
	Size size() const {
		return pages_.size();
	}


	/// @brief clear the book
	inline
	void clear() {
		pages_.clear();
	}


protected: // page management


	/// @brief return the pages in this book
	inline
	Pages const & pages() const {
		return pages_;
	}


	/// @brief return the pages in this book
	inline
	Pages & pages() {
		return pages_;
	}


private: // data


	/// @brief continuous section of lines from a fragment library
	Pages pages_;

};


// Concrete version for PyRosetta
class Book_VallResidues : public Book<core::fragment::picking_old::vall::VallResidues> {};

} // concepts
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_concepts_Book_HH */
