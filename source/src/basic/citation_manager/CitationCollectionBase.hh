// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/CitationCollectionBase.hh
/// @brief Base structure for storing citations.
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_basic_citation_manager_CitationCollectionBase_hh
#define INCLUDED_basic_citation_manager_CitationCollectionBase_hh

#include <basic/citation_manager/CitationCollectionBase.fwd.hh>

// Utility headers
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

namespace basic {
namespace citation_manager {

///@brief Base structure for storing a citation.
/// This is a base class for vitual inheritance. There's not much it actually does itself.
/// Its main purpose is to allow us to put both CitationCollections and UnpublishedModuleInfos into the same list.
class CitationCollectionBase : public utility::VirtualBase {

public:

	virtual
	bool
	operator==( CitationCollectionBase const & other ) const = 0;
};

///@brief A collection of CitationCollectionBases
class CitationCollectionList : public utility::VirtualBase {

public:

	/// @brief The list of citation collections as a vector
	utility::vector1< CitationCollectionBaseCOP > const &
	citations() const;

	/// @brief Are we empty of citations?
	bool
	empty() const;

	// The ability to clear a CitationCollectionList is deliberately omitted.

	/// @brief Add the citation to the list, being aware of duplicates
	void
	add( CitationCollectionBaseCOP const & citr );

	/// @brief Add all the citations in the other list
	void
	add( CitationCollectionList const & other );

	/// @brief Add all the citations in the other list
	template< class T >
	void
	add( utility::vector1< T > const & other ) {
		for ( T const & oentry: other ) {
			add( oentry );
		}
	}

	/// @brief Convenience function for safely getting citation info from a member OP
	template< class T >
	// Only enable template if we have the get_citation_info member (SFINAE)
	typename std::enable_if< std::is_member_function_pointer< decltype(&T::provide_citation_info) >::value, void>::type
	add( utility::pointer::shared_ptr< T > const & ptr ) {
		if ( ptr != nullptr ) {
			ptr->provide_citation_info( *this );
		}
	}

	/// @brief Convenience function for getting citation info from another object.
	template< class T >
	// Only enable template if we have the get_citation_info member (SFINAE)
	typename std::enable_if< std::is_member_function_pointer< decltype(&T::provide_citation_info) >::value, void>::type
	add( T const & obj ) {
		obj.provide_citation_info( *this );
	}

private:

	utility::vector1< CitationCollectionBaseCOP > entries_;
};

} //basic
} //citation_manager

#endif //INCLUDED_basic_citation_manager_CitationCollectionBase_hh

