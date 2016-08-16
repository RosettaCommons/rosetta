// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/concepts/Library.hh
/// @brief  class demonstrating the Library concept
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_concepts_Library_hh
#define INCLUDED_core_fragment_picking_old_concepts_Library_hh


// unit headers
#include <core/fragment/picking_old/concepts/Library.fwd.hh>

// project headers
#include <core/types.hh>

#include <core/fragment/picking_old/vall/VallSection.hh>

namespace core {
namespace fragment {
namespace picking_old {
namespace concepts {


/// @brief class demonstrating the Library concept
/// @remarks Class is useable as a concrete implementation of a Library.
template< typename Books >
class Library {


public: // typedefs


	typedef core::Size Size;
	typedef typename Books::Book Book;
	typedef typename Books::const_iterator BookConstIterator;
	typedef typename Books::iterator BookIterator;


public: // constructor


	/// @brief default constructor
	inline
	Library() {}


	/// @brief Books constructor
	inline
	Library( Books const & books ) : books_( books ) {}


	/// @brief default destructor
	inline
	~Library() {}


public: // iterators


	/// @brief return an iterator that points to the first Book in this library
	inline
	BookConstIterator begin() const {
		return books_.begin();
	}


	/// @brief return an iterator that points to the first Book in this library
	inline
	BookIterator begin() {
		return books_.begin();
	}


	/// @brief return an iterator that points beyond the last Book in this library
	inline
	BookConstIterator end() const {
		return books_.end();
	}


	/// @brief return an iterator that points beyond the last Book in this library
	inline
	BookIterator end()  {
		return books_.end();
	}


public: // book management


	/// @brief number of books in the library
	inline
	Size size() const {
		return books_.size();
	}


	/// @brief clear the library
	inline
	void clear() {
		books_.clear();
	}


protected: // book management


	/// @brief return the books in this library
	inline
	Books const & books() const {
		return books_;
	}


	/// @brief return the books in this library
	inline
	Books & books() {
		return books_;
	}


private: // data


	/// @brief all sections of lines from a fragment library
	Books books_;


};

// Concrete version for PyRosetta
class Library_VallSections: public Library<core::fragment::picking_old::vall::VallSections> {};


} // concepts
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_concepts_Library_HH */
