// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/gen/VallFragmentGen.hh
/// @brief  base class Vall ExtentGenerator
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_gen_VallFragmentGen_hh
#define INCLUDED_core_fragment_picking_old_vall_gen_VallFragmentGen_hh

// unit headers
#include <core/fragment/picking_old/vall/gen/VallFragmentGen.fwd.hh>

// type headers
#include <core/types.hh>

// project headers
#include <core/fragment/picking_old/concepts/Extent.hh>
#include <core/fragment/picking_old/vall/VallSection.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace gen {


/// @brief  base class Vall ExtentGenerator
class VallFragmentGen : public utility::pointer::ReferenceCount {


private: // typedefs


	typedef utility::pointer::ReferenceCount Super;


public: // typedefs


	typedef core::Size Size;


public: // concept typedefs


	/// @brief typedef for ExtentGenerator concept
	typedef core::fragment::picking_old::concepts::Extent< VallSection::PageConstIterator > Extent;


	/// @brief typedef for ExtentGenerator concept
	typedef Extent::PageIterator PageIterator;


public: // concept translation typedefs


	typedef PageIterator VallResidueIterator;


public: // construct/destruct


	/// @brief default constructor
	VallFragmentGen();


	/// @brief copy constructor
	VallFragmentGen( VallFragmentGen const & rval );


	/// @brief default destructor
	virtual
	~VallFragmentGen();


public: // copy assignment


	/// @brief copy assignment
	VallFragmentGen & operator =( VallFragmentGen const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	VallFragmentGenOP clone() const = 0;


public: // extent generation


	/// @brief return the desired fragment extent
	/// @return valid (true) extent if extent will be evaluated, invalid
	///  (false) extent otherwise
	virtual
	Extent operator ()( VallResidueIterator extent_begin, VallResidueIterator section_end ) const = 0;


};


} // namespace gen
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core

#endif /* INCLUDED_core_fragment_picking_old_vall_gen_VallFragmentGen_HH */
