// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/gen/LengthGen.hh
/// @brief  default constant length fragment VallExtentGenerator
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_gen_LengthGen_hh
#define INCLUDED_core_fragment_picking_old_vall_gen_LengthGen_hh

// unit headers
#include <core/fragment/picking_old/vall/gen/LengthGen.fwd.hh>
#include <core/fragment/picking_old/vall/gen/VallFragmentGen.hh>

#include <utility/vector1.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace gen {


/// @brief the default constant length fragment Vall ExtentGenerator
/// @remarks assumes that Pages in the Book are stored in a container
///  capable of returning a RandomAccessIterator, such as std::vector
class LengthGen : public VallFragmentGen {


private: // typedefs


	typedef VallFragmentGen Super;


public: // typedefs


	typedef Super::Size Size;


public: // concept typedefs


	/// @brief typedef for ExtentGenerator concept
	typedef Super::Extent Extent;


	/// @brief typedef for ExtentGenerator concept
	typedef Super::PageIterator PageIterator;


public: // concept translation typedefs


	typedef PageIterator VallResidueIterator;


public: // construct/destruct


	/// @brief default constructor
	LengthGen();


	/// @brief fragment size constructor
	/// @param[in] frag_size the desired length of the fragment
	LengthGen( Size const frag_size );


	/// @brief copy constructor
	LengthGen( LengthGen const & rval );


	/// @brief default destructor
	virtual
	~LengthGen();


public: // copy assignment


	/// @brief copy assignment
	LengthGen & operator =( LengthGen const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	VallFragmentGenOP clone() const;


public: // extent generation


	/// @brief return the desired fragment extent w/requested fragment size
	/// @return valid (true) extent if the end of the extent does not go past the
	///  section_end, invalid (false) extent otherwise
	/// @remarks we assume VallResidueIterator is a type of RandomAccessIterator, such as
	///  those used in std::vector
	virtual
	Extent operator ()( VallResidueIterator extent_begin, VallResidueIterator section_end ) const;


private: // data


	/// @brief the size of desired fragment
	Size frag_size_;

};


} // namespace gen
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core

#endif /* INCLUDED_core_fragment_picking_old_vall_gen_LengthGen_HH */
