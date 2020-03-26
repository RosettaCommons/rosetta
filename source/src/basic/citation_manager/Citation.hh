// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/Citation.hh
/// @brief Data structure for storing a citation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_basic_citation_manager_Citation_hh
#define INCLUDED_basic_citation_manager_Citation_hh

#include <basic/citation_manager/Citation.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

// STL headers
#include <string>

namespace basic {
namespace citation_manager {

/// @brief Data structure for author names.
class AuthorNames : public utility::VirtualBase {

public:

	/// @brief Basic constructor -- deleted.
	AuthorNames() = delete;

	/// @brief Initialization constructor.
	AuthorNames( std::string const &given_names, std::string const &surname, std::string const &initials );

	AuthorNames( AuthorNames const & ) = default;
	~AuthorNames() override = default;
	AuthorNamesOP clone() const;

	/// @brief Assignment operator.
	AuthorNames & operator=( AuthorNames const & ) = default;

public: //Accessors

	/// @brief Get the given name(s).
	inline std::string const & given_names() const { return given_names_; }

	/// @brief Get the surname(s).
	inline std::string const & surname() const { return surname_; }

	/// @brief Get the initial(s).
	inline std::string const & initials() const { return initials_; }

private:

	std::string given_names_;

	std::string surname_;

	std::string initials_;

};

///@brief Data structure for storing a citation.
class Citation : public utility::VirtualBase {

public:

	/// @brief Constructor.
	Citation() = default;

	/// @brief Initialization constructor.
	/// @details Initializes from input string.  Below is an example of the format:
	/// [BEGIN_CITATION]
	///     [BEGIN_PRIMARY_AUTHORS]
	///         "Bobo" "Dang" "B"
	///         "Haifan" "Wu" "H"
	///         "Vikram Khipple" "Mulligan" "VK"
	///     [END_PRIMARY_AUTHORS]
	///     [BEGIN_COAUTHORS]
	///         "Marco" "Mravic" "M"
	///         "Yibing" "Wu" "Y"
	///         "Thomas" "Lemmin" "T"
	///         "Alexander" "Ford" "A"
	///         "Daniel-Adriano" "Silva" "D-A"
	///         "David" "Baker" "D"
	///     [END_COAUTHORS]
	///     [BEGIN_SENIOR_AUTHORS]
	///         "William" "DeGrado" "WF"
	///     [END_SENIOR_AUTHORS]
	///     [BEGIN_YEAR]
	///         2017
	///     [END_YEAR]
	///     [BEGIN_TITLE]
	///         De novo design of covalently constrained mesosize protein scaffolds with unique tertiary structures
	///     [END_TITLE]
	///     [BEGIN_JOURNAL]
	///         Proc Natl Acad Sci USA
	///     [END_JOURNAL]
	///     [BEGIN_VOLUME_ISSUE_PAGES]
	///         114(41):10852–10857
	///     [END_VOLUME_ISSUE_PAGES]
	///     [BEGIN_DOI]
	///         10.1073/pnas.1710695114
	///     [END_DOI]
	/// [END_CITATION]
	Citation( std::string const & citation_as_string );

	/// @brief Copy constructor.
	Citation(Citation const &) = default;

	~Citation() override = default;

	CitationOP clone() const;

public: //Setters

	/// @brief Completely clear the data in this citation.
	void clear();

	/// @brief Initializes from a string formatted as in the Rosetta database.  See initialization
	/// constructor's description for an example.
	/// @details Calls clear() before initializing.
	void initialize_from_string( std::string const & citation_as_string );

	/// @brief Add a primary author.
	/// @details Must be added in order of appearance.
	/// @note If I were adding "James Tiberius Kirk", I would set given_names="James Tiberius", surname="Kirk", initials="JT".
	void add_primary_author( std::string const &given_names, std::string const &surname, std::string const &initials );

	/// @brief Add a co-author.
	/// @details Must be added in order of appearance.
	/// @note If I were adding "John Fitzgerald Kennedy Jr.", I would set given_names="John Fitzgerald", surname="Kennedy", initials="JF Jr.".
	void add_co_author( std::string const &given_names, std::string const &surname, std::string const &initials );

	/// @brief Add a primary author.
	/// @details Must be added in order of appearance.
	/// @note If I were adding "Mary-Kate Olsen", I would set given_names="Mary-Kate", surname="Olsen", initials="M-K".
	void add_senior_author( std::string const &given_names, std::string const &surname, std::string const &initials );

	/// @brief Set the year, article title, journal title, and the string for the volume/issue/pages.
	void set_year_and_article( int const year_in, std::string const &article_title, std::string const &journal_title, std::string const &volume_issue_pages, std::string const &doi);

public: //Accessors

	/// @brief Get this citation as a string, in a particular format.
	void get_formatted_citation( std::ostream & outstream, CitationFormat const format = CitationFormat::DefaultStyle ) const;

	/// @brief Access the list of primary authors.
	inline utility::vector1< AuthorNames > const & primary_authors() const { return primary_authors_; }

	/// @brief Access the list of coauthors.
	inline utility::vector1< AuthorNames > const & co_authors() const { return co_authors_; };

	/// @brief Access the list of senior authors.
	inline utility::vector1< AuthorNames > const & senior_authors() const { return senior_authors_; };

	/// @brief Get the year.
	inline int year() const { return year_; }

	/// @brief Get the article title.
	inline std::string const & article_title() const { return article_title_; }

	/// @brief Get the journal title.
	inline std::string const & journal_title() const { return journal_title_; }

	/// @brief Get the volume, issue, and pages.
	inline std::string const & volume_issue_pages() const { return volume_issue_pages_; }

	/// @brief Get the DOI.
	inline std::string const & doi() const { return doi_; }

private: //Functions

	/// @brief Get the formatted citation, in default format.
	void get_default_formatted_citation( std::ostream & outstream ) const;

	/// @brief Get the formatted citation, in Nature format.
	void get_nature_formatted_citation( std::ostream & outstream ) const;

	/// @brief Given a list of lines, parse out a list of author names from a block of the format:
	/// [BEGIN_SECTION_NAME]
	///      "Given names" "Surname" "Initials"
	///      "Given names" "Surname" "Initials"
	///      "Given names" "Surname" "Initials"
	/// [END_SECTION_NAME]
	/// @details Stores the result in destination.
	void parse_authors( std::string const & section_name, utility::vector1< AuthorNames > & destination, utility::vector1< std::string > const & lines );

	/// @brief Given a list of lines, parse out an integer from a block of the format:
	/// [BEGIN_YEAR]
	///      1982
	/// [END_YEAR]
	/// @details Stores the result in year_.
	void parse_year( utility::vector1< std::string > const & lines );

	/// @brief Given a list of lines, parse out a string from a block of the format:
	/// [BEGIN_SECTION_NAME]
	///      The string to parse out
	/// [END_SECTION_NAME]
	/// @details Stores the result in destination.
	void parse_string_entry( std::string const & section_name, std::string & destination, utility::vector1< std::string > const & lines );

private: //Data

	utility::vector1< AuthorNames > primary_authors_;
	utility::vector1< AuthorNames > co_authors_;
	utility::vector1< AuthorNames > senior_authors_;
	int year_ = 0;
	std::string article_title_ = "";
	std::string journal_title_ = "";
	std::string volume_issue_pages_ = "";
	std::string doi_ = "";

};


} //basic
} //citation_manager

#endif //INCLUDED_basic_citation_manager_Citation_hh

