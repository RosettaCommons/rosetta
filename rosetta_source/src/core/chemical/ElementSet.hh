// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/ElementSet.hh
/// @author P. Douglas Renfrew (renfrew@unc.edu)


#ifndef INCLUDED_core_chemical_ElementSet_hh
#define INCLUDED_core_chemical_ElementSet_hh


// Unit headers
#include <core/chemical/ElementSet.fwd.hh>
// AUTO-REMOVED #include <core/chemical/Element.hh>

// Project headers

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/exit.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
// Commented by inclean daemon #include <string>
#include <map>

#include <utility/vector1_bool.hh>

#include <core/types.hh>
#include <core/chemical/Element.fwd.hh>




namespace core {
namespace chemical {


/// @brief A set of Elements
///
/// @details This class contains a vector of pointers each of which points to an
/// Element and the vector index is looked up by an element_name string
/// in a map.
///
class ElementSet : public utility::pointer::ReferenceCount {

public:
	ElementSet();
	~ElementSet();

	/// @brief Number of elements in the set
	Size
	n_elements() const
	{
		return elements_.size();
	}


	/// @brief Check if there is an element_type associated with an element_symbol string
	bool
	contains_element_type( std::string const & element_symbol ) const
	{
		std::map< std::string, int >::const_iterator
			iter( element_index_.find( element_symbol ) );
		return iter != element_index_.end();
	}


	/// @brief Lookup the element index by the element_symbol string
	int
	element_index( std::string const & element_symbol ) const
	{
		std::map< std::string, int >::const_iterator
			iter( element_index_.find( element_symbol ) );
		if ( iter == element_index_.end() ) {
			utility_exit_with_message("unrecognized element_symbol "+element_symbol);
		}
		return iter->second;
	}


	/// @brief Lookup an Element by 1-based indexing
	Element const &
	operator[] ( Size const index ) const
	{
		return *( elements_[ index ] );
	}


	/// @brief Load the ElementSet from a file
	void
	read_file( std::string const & filename );

	/// @brief Print all of the symbols of all of the Elements in the set. Usefull for debuging.
	void
	print_all_types();


	// data
private:

	/// @brief element_index_ lookup map
	///
	/// @details element_index_ allows lookup of the element by its symbol
	std::map< std::string, int > element_index_;

	/// @brief a collection of Elements,
	///
	/// @details Element has data of atom properties, and it can be
	/// looked up by element_index.
	utility::vector1< Element* > elements_;

};

} // chemical
} // core

#endif // INCLUDED_core_chemical_ElementSet_HH
