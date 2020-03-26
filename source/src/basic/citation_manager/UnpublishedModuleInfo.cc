// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/UnpublishedModuleInfo.cc
/// @brief Authorship information for an unpublished Rosetta module.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/UnpublishedModuleAuthor.hh>
#include <basic/citation_manager/CitationCollection.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "basic.citation_manager.UnpublishedModuleInfo" );


namespace basic {
namespace citation_manager {

/// @brief Initialization constructor.
/// @details Cannot be CustomType.
UnpublishedModuleInfo::UnpublishedModuleInfo(
	std::string const & module_name,
	CitedModuleType const & module_type
) :
	utility::VirtualBase(),
	module_name_(module_name),
	module_type_(module_type)
{
	runtime_assert( !module_name_.empty() );
	runtime_assert( module_type_ != CitedModuleType::CustomType );
}

/// @brief Initialization constructor for CustomTypes.
UnpublishedModuleInfo::UnpublishedModuleInfo(
	std::string const & module_name,
	std::string const & module_type_name
) :
	utility::VirtualBase(),
	module_name_(module_name),
	module_type_(CitedModuleType::CustomType),
	module_type_name_(module_type_name)
{
	runtime_assert( !module_name_.empty() );
	runtime_assert( !module_type_name_.empty() );
}

/// @brief Initialization constructor with single author.
/// @details Cannot be CustomType.
UnpublishedModuleInfo::UnpublishedModuleInfo(
	std::string const & module_name,
	CitedModuleType const & module_type,
	std::string const & author_name,
	std::string const & affiliation,
	std::string const & email
) :
	utility::VirtualBase(),
	authors_({ utility::pointer::make_shared< UnpublishedModuleAuthor >( author_name, affiliation, email ) }),
	module_name_(module_name),
	module_type_(module_type)
{
	runtime_assert( !module_name_.empty() );
	runtime_assert( module_type_ != CitedModuleType::CustomType );
}

/// @brief Initialization constructor with single author.
/// @details Cannot be CustomType.  This version has a field for notes.
UnpublishedModuleInfo::UnpublishedModuleInfo(
	std::string const & module_name,
	CitedModuleType const & module_type,
	std::string const & author_name,
	std::string const & affiliation,
	std::string const & email,
	std::string const & notes
) :
	utility::VirtualBase(),
	authors_({ utility::pointer::make_shared< UnpublishedModuleAuthor >( author_name, affiliation, email, notes ) }),
	module_name_(module_name),
	module_type_(module_type)
{
	runtime_assert( !module_name_.empty() );
	runtime_assert( module_type_ != CitedModuleType::CustomType );
}

/// @brief Comparison operator.
bool
UnpublishedModuleInfo::operator==(
	UnpublishedModuleInfo const & other
) const {
	return ( (other.module_name_ == this->module_name_ ) && (other.module_type_ == this->module_type_ ) && (other.module_type_name_ == this->module_type_name_) );
}

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
UnpublishedModuleInfoOP
UnpublishedModuleInfo::clone() const {
	return utility::pointer::make_shared< UnpublishedModuleInfo >( *this );
}

/// @brief Get the unpublished module's type.
std::string const &
UnpublishedModuleInfo::module_type() const {
	if ( module_type_ == CitedModuleType::CustomType ) return module_type_name_;
	return CitationCollection::get_enumerated_module_type_name( module_type_ );
}

/// @brief Add an author for this unpublished module.
/// @details The author name must be nonempty; the other fields may be empty.
void
UnpublishedModuleInfo::add_author(
	std::string const & author_name,
	std::string const & affiliation,
	std::string const & email
) {
	authors_.push_back( utility::pointer::make_shared< UnpublishedModuleAuthor >( author_name, affiliation, email ) );
}

/// @brief Add an author for this unpublished module.
/// @details The author name must be nonempty; the other fields may be empty.  This version has a field for notes.
void
UnpublishedModuleInfo::add_author(
	std::string const & author_name,
	std::string const & affiliation,
	std::string const & email,
	std::string const & notes
) {
	authors_.push_back( utility::pointer::make_shared< UnpublishedModuleAuthor >( author_name, affiliation, email, notes ) );
}

/// @brief List authors in the format:
/// Name, Institution <email>
/// Name, Institution <email>
/// Name, Institution <email>
/// ...
void
UnpublishedModuleInfo::get_author_list(
	std::ostream & outstream
) const {
	for ( auto const & author : authors_ ) {
		outstream << author->author_name();
		if ( !(author->affiliation().empty()) ) outstream << ", " << author->affiliation();
		if ( !(author->email_address().empty()) ) outstream << " <" << author->email_address() << ">";
		if ( !(author->notes().empty()) ) outstream << "  (" << author->notes() << ")";
		outstream << "\n";
	}
}


////////////////////////////////// HELPER FUNCTIONS //////////////////////////////////
///////////////////////////////// (Not part of class) ////////////////////////////////

/// @brief Helper function that compares two lists of UnpublishedModuleInfoCOPs and appends
/// all the ones from list A that aren't already in list B onto the end of list B.
void
merge_into_unpublished_collection_vector(
	utility::vector1< UnpublishedModuleInfoCOP > const & source,
	utility::vector1< UnpublishedModuleInfoCOP > & destination
) {
	if ( source.empty() ) return;
	for ( auto const & moduleinfo : source ) {
		bool already_present( false );
		for ( auto const & destinfo : destination ) {
			if ( (*moduleinfo) == (*destinfo) ) {
				already_present = true;
				break;
			}
		}
		if ( already_present ) continue;
		destination.push_back( moduleinfo );
	}
}

} //citation_manager
} //basic
