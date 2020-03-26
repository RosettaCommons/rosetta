// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/UnpublishedModuleAuthor.hh
/// @brief Authorship information for a single author in an unpublished Rosetta module.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_basic_citation_manager_UnpublishedModuleAuthor_hh
#define INCLUDED_basic_citation_manager_UnpublishedModuleAuthor_hh

#include <basic/citation_manager/UnpublishedModuleAuthor.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

// C++ headers
#include <string>

namespace basic {
namespace citation_manager {

/// @brief Authorship information for a single author in an unpublished Rosetta module.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class UnpublishedModuleAuthor : public utility::VirtualBase {

public:

	/// @brief Default constructor, deleted.
	UnpublishedModuleAuthor() = delete;

	/// @brief Initialization constructor.
	UnpublishedModuleAuthor( std::string const & author_name, std::string const & affiliation, std::string const & email_address );

	/// @brief Initialization constructor with notes.
	UnpublishedModuleAuthor( std::string const & author_name, std::string const & affiliation, std::string const & email_address, std::string const & notes );

	/// @brief Copy constructor.
	UnpublishedModuleAuthor(UnpublishedModuleAuthor const &) = default;

	/// @brief Destructor.
	~UnpublishedModuleAuthor() override = default;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	UnpublishedModuleAuthorOP clone() const;

public: //Accessors

	/// @brief Get the author's name.
	inline std::string const & author_name() const { return author_name_; }

	/// @brief Get his or her institutional affiliation.
	inline std::string const & affiliation() const { return affiliation_; }

	/// @brief Get his or her e-mail address.
	inline std::string const & email_address() const { return email_address_; }

	/// @brief Get the notes, if any.
	inline std::string const & notes() const { return notes_; }

private:

	/// @brief The author's name.
	std::string author_name_;

	/// @brief His or her institutional affiliation.
	std::string affiliation_;

	/// @brief His or her e-mail address.
	std::string email_address_;

	/// @brief An optional notes field (e.g. "Expanded functionality for noncanonicals.").  Notes will
	/// be printed in parentheses (which should NOT be included).
	std::string notes_;

};

} //citation_manager
} //basic

#endif //INCLUDED_basic_citation_manager_UnpublishedModuleAuthor_hh
