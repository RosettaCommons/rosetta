// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/UnpublishedModuleInfo.hh
/// @brief Authorship information for an unpublished Rosetta module.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_basic_citation_manager_UnpublishedModuleInfo_hh
#define INCLUDED_basic_citation_manager_UnpublishedModuleInfo_hh

#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>
#include <basic/citation_manager/UnpublishedModuleAuthor.fwd.hh>
#include <basic/citation_manager/CitationCollection.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

namespace basic {
namespace citation_manager {

/// @brief Authorship information for an unpublished Rosetta module.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class UnpublishedModuleInfo : public utility::VirtualBase {

public:

	/// @brief Default constructor, explicitly deleted.
	UnpublishedModuleInfo() = delete;

	/// @brief Initialization constructor.
	/// @details Cannot be CustomType.
	UnpublishedModuleInfo( std::string const & module_name, CitedModuleType const & module_type );

	/// @brief Initialization constructor with single author.
	/// @details Cannot be CustomType.
	UnpublishedModuleInfo(
		std::string const & module_name,
		CitedModuleType const & module_type,
		std::string const & author_name,
		std::string const & affiliation,
		std::string const & email
	);

	/// @brief Initialization constructor with single author.
	/// @details Cannot be CustomType.  This version has a field for notes.
	UnpublishedModuleInfo(
		std::string const & module_name,
		CitedModuleType const & module_type,
		std::string const & author_name,
		std::string const & affiliation,
		std::string const & email,
		std::string const & notes
	);

	/// @brief Initialization constructor for CustomTypes.
	UnpublishedModuleInfo( std::string const & module_name, std::string const & module_type_name );

	/// @brief Copy constructor.
	UnpublishedModuleInfo(UnpublishedModuleInfo const &) = default;

	/// @brief Destructor.
	~UnpublishedModuleInfo() override = default;

	/// @brief Comparison operator.
	bool operator==( UnpublishedModuleInfo const & other ) const;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	UnpublishedModuleInfoOP clone() const;

public: //Accessors

	/// @brief Get the unpublished module's name.
	inline std::string const & module_name() const { return module_name_; }

	/// @brief Get the unpublished module's type.
	std::string const & module_type() const;

	/// @brief Add an author for this unpublished module.
	/// @details The author name must be nonempty; the other fields may be empty.
	void add_author( std::string const & author_name, std::string const & affiliation, std::string const & email );

	/// @brief Add an author for this unpublished module.
	/// @details The author name must be nonempty; the other fields may be empty.  This version has a field for notes.
	void add_author( std::string const & author_name, std::string const & affiliation, std::string const & email, std::string const & notes );

	/// @brief List authors in the format:
	/// Name, Institution <email>
	/// Name, Institution <email>
	/// Name, Institution <email>
	/// ...
	void get_author_list( std::ostream & outstream ) const;

private:

	/// @brief Information about all of the authors.
	utility::vector1< UnpublishedModuleAuthorCOP > authors_;

	/// @brief The name of the unpublished module.
	std::string module_name_;

	/// @brief The type of the unpublished module.
	CitedModuleType module_type_;

	/// @brief The name of the unpublished module type, if it's a CustomType.
	std::string module_type_name_;

};

/// @brief Helper function that compares two lists of UnpublishedModuleInfoCOPs and appends
/// all the ones from list A that aren't already in list B onto the end of list B.
void
merge_into_unpublished_collection_vector(
	utility::vector1< UnpublishedModuleInfoCOP > const & source,
	utility::vector1< UnpublishedModuleInfoCOP > & destination
);

} //citation_manager
} //basic

#endif //INCLUDED_basic_citation_manager_UnpublishedModuleInfo_hh
