// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/UnpublishedModuleAuthor.cc
/// @brief Authorship information for a single author in an unpublished Rosetta module.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <basic/citation_manager/UnpublishedModuleAuthor.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "basic.citation_manager.UnpublishedModuleAuthor" );


namespace basic {
namespace citation_manager {

/// @brief Initialization constructor.
UnpublishedModuleAuthor::UnpublishedModuleAuthor(
	std::string const & author_name,
	std::string const & affiliation,
	std::string const & email_address
) :
	utility::VirtualBase(),
	author_name_( author_name ),
	affiliation_( affiliation ),
	email_address_( email_address )
{
	runtime_assert_string_msg( !author_name.empty(), "Error in UnpublishedModuleAuthor constructor: An author must have a name!" );
}

/// @brief Initialization constructor with notes.
UnpublishedModuleAuthor::UnpublishedModuleAuthor(
	std::string const & author_name,
	std::string const & affiliation,
	std::string const & email_address,
	std::string const & notes
) :
	utility::VirtualBase(),
	author_name_( author_name ),
	affiliation_( affiliation ),
	email_address_( email_address ),
	notes_(notes)
{
	runtime_assert_string_msg( !author_name.empty(), "Error in UnpublishedModuleAuthor constructor: An author must have a name!" );
}

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
UnpublishedModuleAuthorOP
UnpublishedModuleAuthor::clone() const {
	return utility::pointer::make_shared< UnpublishedModuleAuthor >( *this );
}

} //citation_manager
} //basic
