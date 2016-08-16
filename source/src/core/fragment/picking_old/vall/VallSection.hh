// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/VallSection.hh
/// @brief  class implementing the Book concept for a continuous section of lines in the Vall library
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_VallSection_hh
#define INCLUDED_core_fragment_picking_old_vall_VallSection_hh

// unit headers
#include <core/fragment/picking_old/vall/VallSection.fwd.hh>

// package headers
#include <core/fragment/picking_old/concepts/Book.hh>
#include <core/fragment/picking_old/vall/VallResidue.hh>

#include <utility/vector1.hh>


// utility headers


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {


/// @brief  class implementing the Book concept for a continuous section of lines in the Vall library
class VallSection : public core::fragment::picking_old::concepts::Book< VallResidues >{


private: // typedefs


	typedef core::fragment::picking_old::concepts::Book< VallResidues > Super;


public: // concept typedefs


	/// @brief typedef for Book concept
	typedef Super::PageConstIterator PageConstIterator;


	/// @brief typedef for Book concept
	typedef Super::PageIterator PageIterator;


public: // concept translation typedefs


	typedef PageConstIterator VallResidueConstIterator;
	typedef PageIterator VallResidueIterator;


public: // construct/destruct


public: // operators


	/// @brief access a specific VallResidue within this section, 1-based indexing
	inline
	VallResidue const & operator[]( Size const idx ) const {
		return residues()[ idx ];
	}


public: // page management


	/// @brief append a residue to the end of this section
	inline
	void append_residue( VallResidue const & r ) {
		residues().push_back( r );
		residues().back().position_index( residues().size() );
		residues().back().section_index( index_ );
	}


	/// @brief re-index VallResidues within this library so they number 1 .. n
	inline
	void reindex_residues() {
		Size count = 0;
		for ( VallResidueIterator i = begin(), ie = end(); i != ie; ++i ) {
			i->position_index( ++count );
		}
	}


public: // accessors


	/// @brief stores the 1-based indexing for accessing this section
	///  via VallLibrary::operator []
	inline
	Size index() const {
		return index_;
	}


public: // mutators


	/// @brief stores the 1-based indexing for accessing this section
	///  via VallLibrary::operator []
	inline
	void index( Size const idx ) {
		index_ = idx;
		for ( VallResidueIterator r = begin(), re = end(); r != re; ++r ) {
			r->section_index( idx );
		}
	}

public: // memory


	/// @brief ensure storage currently allocated for residues is at least the
	///  given size
	inline
	void reserve( Size const n ) {
		residues().reserve( n );
	}


	/// @brief tighten memory usage
	/// @details Shrinks the capacity of the VallResidues by performing
	///  a swap, so this operation will effectively require a
	///  similar amount of memory as the current residue object.
	inline
	void tighten_memory() {
		if ( residues().size() < residues().capacity() ) {
			VallResidues( residues() ).swap( residues() );
		}
	}


private: // concept translation


	/// @brief return residues container
	inline
	VallResidues const & residues() const {
		return pages();
	}


	/// @brief return residues container
	inline
	VallResidues & residues() {
		return pages();
	}


private: // data


	/// @brief stores the 1-based indexing for accessing this section
	///  via VallLibrary::operator []
	Size index_;


};


/// @brief wrapper for a collection of VallSection
class VallSections : public utility::vector1< VallSection > {


public: // concept typedefs


	/// @brief Book typedef required to satisfy concept
	typedef VallSection Book;


};


} // vall
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_vall_VallSection_HH */
