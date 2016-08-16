// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/gen/SecStructGen.hh
/// @brief  Generator that requires fragments to have a specific secondary
///         structure string.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_gen_SecStructGen_hh
#define INCLUDED_core_fragment_picking_old_vall_gen_SecStructGen_hh

// unit headers
#include <core/fragment/picking_old/vall/gen/SecStructGen.fwd.hh>
#include <core/fragment/picking_old/vall/gen/VallFragmentGen.hh>

#include <utility/vector1.hh>


// C++ headers


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {
namespace gen {


/// @brief Generator that requires fragments to have a specific secondary
///  structure string.
/// @remarks assumes that Pages in the Book are stored in a container
///  capable of returning a RandomAccessIterator, such as std::vector
class SecStructGen : public VallFragmentGen {


private: // typedefs


	typedef VallFragmentGen Super;


public: // typedefs


	typedef Super::Size Size;
	typedef std::string String;


public: // concept typedefs


	/// @brief typedef for ExtentGenerator concept
	typedef Super::Extent Extent;


	/// @brief typedef for ExtentGenerator concept
	typedef Super::PageIterator PageIterator;


public: // concept translation typedefs


	typedef PageIterator VallResidueIterator;


public: // construct/destruct


	/// @brief default constructor
	SecStructGen();


	/// @brief secondary structure string constructor
	/// @param[in] ss the required secondary structure string of the fragment
	SecStructGen( String const & ss );


	/// @brief copy constructor
	SecStructGen( SecStructGen const & rval );


	/// @brief default destructor
	virtual
	~SecStructGen();


public: // copy assignment


	/// @brief copy assignment
	SecStructGen & operator =( SecStructGen const & rval );


public: // virtual constructors


	/// @brief clone this object
	virtual
	VallFragmentGenOP clone() const;


public: // extent generation


	/// @brief return the desired fragment extent w/ length equal to the
	///  secondary structure string
	/// @return Valid (true) extent if the extent has exactly the required
	///  secondary structure string and the end of the extent does not go past
	///  section_end.  Invalid (false) extent otherwise.
	/// @remarks we assume VallResidueIterator is a type of RandomAccessIterator, such as
	///  those used in std::vector
	virtual
	Extent operator ()( VallResidueIterator extent_begin, VallResidueIterator section_end ) const;


private: // data


	/// @brief the required secondary structure of desired fragment
	String ss_;

};


} // namespace gen
} // namespace vall
} // namespace picking_old
} // namespace fragment
} // namespace core

#endif /* INCLUDED_core_fragment_picking_old_vall_gen_SecStructGen_HH */
