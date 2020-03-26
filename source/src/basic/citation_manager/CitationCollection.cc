// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/CitationCollection.cc
/// @brief A class for keeping track of a collection of citations for a particular Rosetta module.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/Citation.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <map>

static basic::Tracer TR( "basic.citation_manager.CitationCollection" );


namespace basic {
namespace citation_manager {

/// @brief Options constructor.  The module_type must not be CustomType.
CitationCollection::CitationCollection(
	std::string const & module_name,
	CitedModuleType const module_type
) :
	utility::VirtualBase(),
	module_name_(module_name),
	module_type_(module_type)
{
	runtime_assert_string_msg( !module_name.empty(), "Error in constructor for CitationCollection: The name of the module for which citations are provided cannot be empty.");
	runtime_assert( module_type != CitedModuleType::CustomType );
}

/// @brief Options constructor for custom module types.  A string must be provided for the
/// name of the module type.
CitationCollection::CitationCollection(
	std::string const & module_name,
	std::string const & module_type_name
) :
	utility::VirtualBase(),
	module_name_( module_name ),
	module_type_( CitedModuleType::CustomType ),
	module_type_name_( module_type_name )
{
	runtime_assert_string_msg( !module_name.empty(), "Error in constructor for CitationCollection: The name of the module for which citations are provided cannot be empty.");
	runtime_assert_string_msg( !module_type_name.empty(), "Error in constructor for CitationCollection: The name of a custom module type cannot be empty.");
}

/// @brief Clone operator: copy this object and return an owning pointer to the copy.
CitationCollectionOP
CitationCollection::clone() const {
	return utility::pointer::make_shared< CitationCollection >( *this );
}

/// @brief Comparison operator.
/// @details Compares only module name and module type, not reference list.
bool
CitationCollection::operator==(
	CitationCollection const & other
) const {
	return ( (other.module_type_ == module_type_) && ( other.module_name_ == module_name_ ) && (other.module_type_name_ == module_type_name_ ) );
}

/// @brief Add a citation without cloning it.
/// @details Does nothing if citation_in is nullptr.
void
CitationCollection::add_citation(
	CitationCOP citation_in
) {
	if ( citation_in == nullptr ) return; //Do nothing if nullptr is passed in.

	for ( auto const & citation : citations_ ) {
		if ( citation.get() == citation_in.get() ) return; //Do nothing if the citation is already in the collection.
	}

	citations_.push_back( citation_in );
}

/// @brief Get this module's type.
std::string const &
CitationCollection::module_type() const {
	if ( module_type_ == CitedModuleType::CustomType ) return module_type_name_;
	return get_enumerated_module_type_name( module_type_ );
}

/// @brief Get the module type name.
/// @details The module_type cannot be CitedModuleType::CustomType (or this will throw).
/*static*/ std::string const &
CitationCollection::get_enumerated_module_type_name(
	CitedModuleType const module_type
) {
	runtime_assert( module_type != CitedModuleType::CustomType );
	static const std::map< CitedModuleType, std::string > module_type_strings{
		{ CitedModuleType::Mover, "Mover" },
		{ CitedModuleType::Filter, "Filter" },
		{ CitedModuleType::ScoreTerm, "ScoreTerm" },
		{ CitedModuleType::ResidueSelector, "ResidueSelector" },
		{ CitedModuleType::TaskOperation, "TaskOperation" },
		{ CitedModuleType::PackerPalette, "PackerPalette" },
		{ CitedModuleType::ScoreFunction, "ScoreFunction" },
		{ CitedModuleType::EnergyMethod, "EnergyMethod" },
		{ CitedModuleType::SimpleMetric, "SimpleMetric" },
		{ CitedModuleType::Singleton, "Singleton" },
		{ CitedModuleType::Application, "Application" }
		};
	debug_assert( module_type_strings.count( module_type ) != 0 );
	return module_type_strings.at( module_type );
}

/// @brief Get the citations, written out prettily to a string stream.
/// @details Multiple citations will be on multiple lines, but there's no
/// newline at the start or the end of the output.
void
CitationCollection::get_citations_formatted(
	std::ostream & outstream,
	CitationFormat const format /*= CitationFormat::DefaultStyle*/
) const {
	if ( citations_.empty() ) return; //Do nothing if we have no citations.
	platform::Size const total_citations( citations_.size() );
	platform::Size count(0);
	for ( auto const & citation : citations_ ) {
		++count;
		citation->get_formatted_citation( outstream, format );
		outstream << "\n";
		if ( count < total_citations ) outstream << "\n";
	}
}

/// @brief Does this collection have more than one citation?
bool
CitationCollection::has_multiple_citations() const {
	return (citations_.size() > 1);
}


//////////////////////////////// UTILITY FUNCTIONS ////////////////////////////////

/// @brief Utility function to determine whether a particular CitationCollection is already in a vector.
bool
has_citation_collection(
	CitationCollection const & collection,
	utility::vector1< CitationCollectionCOP > const & vec
) {
	for ( CitationCollectionCOP const & entry : vec ) {
		if ( (*entry) == collection ) return true;
	}
	return false;
}

/// @brief Utility function that takes all of the citation collections from vector A that are NOT in vector
/// B and appends them to vector B.
void
merge_into_citation_collection_vector(
	utility::vector1< CitationCollectionCOP > const & source_vector,
	utility::vector1< CitationCollectionCOP > & destination_vector
) {
	if ( source_vector.empty() ) return;
	for ( CitationCollectionCOP const & collection : source_vector ) { // const & so as not to increment reference count.
		if ( !has_citation_collection( *collection, destination_vector ) ) {
			destination_vector.push_back( collection ); // COP is copied here, so this is only place where reference count is incremented.
		}
	}
}


} //basic
} //citation_manager






