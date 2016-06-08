// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/vall/VallLibrary.hh
/// @brief  Vall fragment library
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_VallLibrary_hh
#define INCLUDED_core_fragment_picking_old_vall_VallLibrary_hh

// unit headers
#include <core/fragment/picking_old/vall/VallLibrary.fwd.hh>

// package headers
#include <core/fragment/picking_old/concepts/Library.hh>
#include <core/fragment/picking_old/vall/VallSection.hh>

#include <utility/vector1.hh>


// utility headers


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {


/// @brief Vall fragment library
class VallLibrary : public core::fragment::picking_old::concepts::Library< VallSections > {


private: // typedefs


	typedef core::fragment::picking_old::concepts::Library< VallSections > Super;


public: // typedefs


	typedef Super::Size Size;


public: // concept typedefs


	/// @brief typedef for Library concept
	typedef Super::BookConstIterator BookConstIterator;


	/// @brief typedef for Library concept
	typedef Super::BookIterator BookIterator;


public: // concept translation typedefs


	typedef BookConstIterator VallSectionConstIterator;
	typedef BookIterator VallSectionIterator;


public: // operators


	/// @brief access a specific VallSection within this section, 1-based indexing
	inline
	VallSection const & operator[]( Size idx ) const {
		return sections()[ idx ];
	}


public: // section management


	/// brief add a section to the end of this library
	inline
	void add_section( VallSection const & section ) {
		sections().push_back( section );
		sections().back().index( sections().size() );
	}


	/// @brief access the last section in thie library
	inline
	VallSection const & last_section() const {
		return sections().back();
	}


	/// @brief access the last section in thie library
	inline
	VallSection & last_section() {
		return sections().back();
	}


	/// @brief total number of residues in the library
	inline
	Size n_residues() const {
		Size n = 0;
		for ( VallSectionConstIterator i = begin(), ie = end(); i != ie; ++i ) {
			n += i->size();
		}
		return n;
	}


	/// @brief re-index VallSections within this library so they number 1 .. n
	inline
	void reindex_sections() {
		Size count = 0;
		for ( VallSectionIterator i = begin(), ie = end(); i != ie; ++i ) {
			i->index( ++count );
		}
	}


public: // memory


	/// @brief ensure storage currently allocated for VallSections is at least the
	///  given size
	inline
	void reserve( Size n ) {
		sections().reserve( n );
	}


	/// @brief tighten memory usage
	/// @param sections_only if true, only tightens memory per-section and skips
	///  the operation for the whole library
	/// @details Shrinks the capacity of the library by calling
	///  <tt>VallSection::tighten_memory()</tt> for each section.  If 'sections_only' is
	///  false, will additionally perform a swap of the library, effectively
	///  requiring a similar amount of memory as the current library object.
	/// @warning The Vall libraries can be very large, so set 'sections_only'
	///  false only if you won't be running on a constrained memory platform.
	inline
	void tighten_memory( bool const sections_only = true ) {
		for ( VallSectionIterator i = begin(), ie = end(); i != ie; ++i ) {
			i->tighten_memory();
		}

		// the following could require a lot of memory
		if ( !sections_only ) {
			VallSections( sections() ).swap( sections() );
		}
	}


private: // concept translation


	/// @brief return sections container
	inline
	VallSections const & sections() const {
		return books();
	}


	/// @brief return sections container
	inline
	VallSections & sections() {
		return books();
	}


};


} // vall
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_vall_VallLibrary_HH */
