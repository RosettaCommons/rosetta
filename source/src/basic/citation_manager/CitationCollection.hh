// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/CitationCollection.hh
/// @brief A class for keeping track of a collection of citations for a particular Rosetta module.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_basic_citation_manager_CitationCollection_hh
#define INCLUDED_basic_citation_manager_CitationCollection_hh

#include <basic/citation_manager/CitationCollection.fwd.hh>
#include <basic/citation_manager/Citation.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

namespace basic {
namespace citation_manager {

/// @brief A class for keeping track of a collection of citations for a particular Rosetta module.
class CitationCollection : public utility::VirtualBase {

public:

	/// @brief Default constructor, deleted.
	CitationCollection() = delete;

	/// @brief Options constructor.  The module type must not be CustomType.
	CitationCollection( std::string const & module_name, CitedModuleType const module_type );

	/// @brief Options constructor for custom module types.  A string must be provided for the
	/// name of the module type.
	CitationCollection( std::string const & module_name, std::string const & module_type_name );

	/// @brief Destructor.
	~CitationCollection() override = default;

	/// @brief Copy constructor.
	CitationCollection(CitationCollection const &) = default;

	/// @brief Assignment operator.
	CitationCollection & operator=( CitationCollection const & ) = default;

	/// @brief Clone operator: copy this object and return an owning pointer to the copy.
	CitationCollectionOP clone() const;

	/// @brief Comparison operator.
	/// @details Compares only module name and module type, not reference list.
	bool operator==( CitationCollection const & other ) const;

public:

	/// @brief Add a citation without cloning it.
	/// @details Does nothing if citation_in is nullptr.
	void add_citation( CitationCOP citation_in );

	/// @brief Get this module's name.
	inline std::string const & module_name() const { return module_name_; }

	/// @brief Get this module's type.
	std::string const & module_type() const;

	/// @brief Get the module type name.
	/// @details The module_type cannot be CitedModuleType::CustomType (or this will throw).
	static std::string const & get_enumerated_module_type_name( CitedModuleType const module_type );

	/// @brief Get the citations, written out prettily to a string stream.
	/// @details Multiple citations will be on multiple lines, but there's no
	/// newline at the start or the end of the output.
	void
	get_citations_formatted(
		std::ostream & outstream,
		CitationFormat const format = CitationFormat::DefaultStyle
	) const;

	/// @brief Does this collection have more than one citation?
	bool has_multiple_citations() const;

private:

	/// @brief Citations in this citation collection.
	utility::vector1< CitationCOP > citations_;

	/// @brief The name of the module for which citations are being provided.
	std::string module_name_;

	/// @brief The type of module for which citations are being provided.
	CitedModuleType module_type_ = CitedModuleType::Mover;

	/// @brief For custom module types, provide a string.
	std::string module_type_name_;

};

/// @brief Utility function to determine whether a particular CitationCollection is already in a vector.
bool
has_citation_collection(
	CitationCollection const & collection,
	utility::vector1< CitationCollectionCOP > const & vec
);

/// @brief Utility function that takes all of the citation collections from vector A that are NOT in vector
/// B and appends them to vector B.
void
merge_into_citation_collection_vector(
	utility::vector1< CitationCollectionCOP > const & source_vector,
	utility::vector1< CitationCollectionCOP > & destination_vector
);

} //basic
} //citation_manager



#endif //INCLUDED_basic_citation_manager_CitationCollection_hh





