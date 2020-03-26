// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/Citation.cc
/// @brief Data structure for storing a citation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#include <basic/citation_manager/Citation.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/memory.hh>
#include <utility/string_util.hh>

static basic::Tracer TR( "basic.citation_manager.Citation" );

namespace basic {
namespace citation_manager {

///////////////////////////////////// class AuthorNames /////////////////////////////////////////////////

/// @brief Initialization constructor.
AuthorNames::AuthorNames(
	std::string const &given_names,
	std::string const &surname,
	std::string const &initials
) :
	utility::VirtualBase(),
	given_names_(given_names),
	surname_(surname),
	initials_(initials)
{
	runtime_assert_string_msg( !surname_.empty(), "Error in AuthorNames class constructor: An author must have a surname.");
}

AuthorNamesOP
AuthorNames::clone() const {
	return utility::pointer::make_shared< AuthorNames >(*this);
}

///////////////////////////////////// class Citation /////////////////////////////////////////////////

/// @brief Initialization constructor.
Citation::Citation(
	std::string const & citation_as_string
) :
	utility::VirtualBase()
{
	initialize_from_string( citation_as_string );
}

CitationOP
Citation::clone() const {
	return utility::pointer::make_shared< Citation >( *this );
}

/// @brief Completely clear the data in this citation.
void
Citation::clear() {
	primary_authors_.clear();
	co_authors_.clear();
	senior_authors_.clear();
	year_ = 0;
	article_title_ = "";
	journal_title_ = "";
	volume_issue_pages_ = "";
	doi_ = "";
}

/// @brief Initializes from a string formatted as in the Rosetta database.  See initialization
/// constructor's description for an example.
/// @details Calls clear() before initializing.
void
Citation::initialize_from_string(
	std::string const & citation_as_string
) {
	std::string const errmsg( "Error in Citation::initialize_from_string():  " );
	clear();
	utility::vector1< std::string > const lines( utility::split_by_newlines( citation_as_string ) );

	runtime_assert_string_msg( lines.size() > 2 && utility::strip(lines[1]) == "[BEGIN_CITATION]" && utility::strip(lines[lines.size()]) == "[END_CITATION]",
		errmsg + "A block of text beginning with \"[BEGIN_CITATION]\" and ending with \"[END_CITATION]\" must be passed to this function.  Received:\n" + citation_as_string );

	parse_authors( "PRIMARY_AUTHORS", primary_authors_, lines );
	parse_authors( "COAUTHORS", co_authors_, lines );
	parse_authors( "SENIOR_AUTHORS", senior_authors_, lines );
	parse_year( lines );
	parse_string_entry( "TITLE", article_title_, lines );
	parse_string_entry( "JOURNAL", journal_title_, lines );
	parse_string_entry( "VOLUME_ISSUE_PAGES", volume_issue_pages_, lines );
	parse_string_entry( "DOI", doi_, lines );

	runtime_assert_string_msg( !doi_.empty(), errmsg + "No DOI was found!  Original string:\n\n" + citation_as_string + "\n\nInput must include a DOI!" );
	runtime_assert_string_msg( !primary_authors_.empty(), errmsg + "No primary authors were found!  Original string:\n\n" + citation_as_string + "\n\nInput must include at least one primary author!" );
}

/// @brief Add a primary author.
/// @details Must be added in order of appearance.
/// @note If I were adding "James Tiberius Kirk", I would set given_names="James Tiberius", surname="Kirk", initials="JT".
void
Citation::add_primary_author(
	std::string const &given_names,
	std::string const &surname,
	std::string const &initials
) {
	primary_authors_.emplace_back( given_names, surname, initials );
}

/// @brief Add a co-author.
/// @details Must be added in order of appearance.
/// @note If I were adding "John Fitzgerald Kennedy Jr.", I would set given_names="John Fitzgerald", surname="Kennedy", initials="JF Jr.".
void
Citation::add_co_author(
	std::string const &given_names,
	std::string const &surname,
	std::string const &initials
) {
	co_authors_.emplace_back( given_names, surname, initials );
}

/// @brief Add a primary author.
/// @details Must be added in order of appearance.
/// @note If I were adding "Mary-Kate Olsen", I would set given_names="Mary-Kate", surname="Olsen", initials="M-K".
void
Citation::add_senior_author(
	std::string const &given_names,
	std::string const &surname,
	std::string const &initials
) {
	senior_authors_.emplace_back( given_names, surname, initials );
}

/// @brief Set the year, article title, journal title, and the string for the volume/issue/pages.
void
Citation::set_year_and_article(
	int const year_in,
	std::string const &article_title,
	std::string const &journal_title,
	std::string const &volume_issue_pages,
	std::string const &doi
) {
	runtime_assert_string_msg( year_in != 0, "Error in Citation::set_year(): There was no year 0!  Please use positive numbers for the common era/anno domini, and negative numbers for before the common era/before Christ.");
	year_ = year_in;
	article_title_ = article_title;
	journal_title_ = journal_title;
	volume_issue_pages_ = volume_issue_pages;
	doi_ = doi;
}


/// @brief Get this citation as a string, in a particular format.
void
Citation::get_formatted_citation(
	std::ostream & outstream,
	CitationFormat const format /* = CitationFormat::DefaultStyle*/
) const {
	switch( format ) {
	case CitationFormat::DefaultStyle :
		get_default_formatted_citation( outstream );
		break;
	case CitationFormat::NatureStyle :
		get_nature_formatted_citation( outstream );
		break;
	default :
		utility_exit_with_message( "Error in Citation::get_formatted_citation(): Unknown format specified!");
	}
}

/// @brief Get the formatted citation, in default format.
void
Citation::get_default_formatted_citation(
	std::ostream & outstream
) const {
	bool text_exists(false);
	bool const multiple_primary_authors( primary_authors_.size() > 1 );
	platform::Size const num_authors( primary_authors_.size() + co_authors_.size() + senior_authors_.size() );
	if ( num_authors > 0 ) {
		text_exists = true;
		int pass(0);
		platform::Size author_count(0);
		utility::vector1< utility::vector1< AuthorNames > const * > const authorlists{ &primary_authors_, &co_authors_, &senior_authors_ };
		for ( utility::vector1< AuthorNames > const * authorlist : authorlists ) {
			++pass;
			if ( !(authorlist->empty()) ) {
				for ( AuthorNames const & author: (*authorlist) ) {
					++author_count;
					if ( multiple_primary_authors && pass == 1 ) outstream << "*";
					outstream << author.surname();
					if ( !author.initials().empty() ) {
						outstream << " " << author.initials();
					}
					if ( author_count < num_authors - 1 ) {
						outstream << ", ";
					} else if ( author_count == num_authors - 1 ) {
						outstream << ", and ";
					} else {
						outstream << ".";
					}
				}
			}
		}
	}
	if ( year_ != 0 ) {
		if ( text_exists ) outstream << "  ";
		outstream << "(" << year_ << ").";
	}
	if ( !article_title_.empty() ) {
		if ( text_exists ) outstream << "  ";
		text_exists = true;
		outstream << article_title_ << ".";
	}
	if ( !journal_title_.empty() ) {
		if ( text_exists ) outstream << "  ";
		text_exists = true;
		if ( journal_title_ != "XXX" ) {
			outstream << journal_title_;
		} else {
			outstream << "(Manuscript under review.)";
		}
	}
	if ( !volume_issue_pages_.empty() && journal_title_ != "XXX" ) {
		if ( text_exists ) outstream << " ";
		text_exists = true;
		outstream << volume_issue_pages_ << ".";
	}
	if ( !doi_.empty() && journal_title_ != "XXX" ) {
		if ( text_exists ) outstream << "  ";
		outstream << "doi: " << doi_ << ".";
		//text_exists = true; //Uncomment if I add anything more.  The cppcheck test doesn't like this.
	}
	if ( multiple_primary_authors ) {
		outstream << "  (*Co-primary authors.)";
	}
}


/// @brief Get the formatted citation, in Nature format.
void
Citation::get_nature_formatted_citation(
	std::ostream & outstream
) const {
	bool text_exists(false);
	platform::Size const num_authors( primary_authors_.size() + co_authors_.size() + senior_authors_.size() );
	if ( num_authors > 0 ) {
		text_exists = true;
		platform::Size author_count(0);
		utility::vector1< utility::vector1< AuthorNames > const * > const authorlists{ &primary_authors_, &co_authors_, &senior_authors_ };
		for ( utility::vector1< AuthorNames > const * authorlist : authorlists ) {
			if ( !(authorlist->empty()) ) {
				for ( AuthorNames const & author: (*authorlist) ) {
					++author_count;
					outstream << author.surname();
					if ( !author.initials().empty() ) {
						outstream << ", " << author.initials();
					}
					if ( author_count < num_authors - 1 ) {
						outstream << ", ";
					} else if ( author_count == num_authors - 1 ) {
						outstream << " & ";
					} else {
						outstream << ".";
					}
				}
			}
		}
	}
	if ( !article_title_.empty() ) {
		if ( text_exists ) outstream << "  ";
		text_exists = true;
		outstream << article_title_ << ".";
	}
	if ( !journal_title_.empty() ) {
		if ( text_exists ) outstream << "  ";
		text_exists = true;
		outstream << journal_title_;
	}
	if ( !volume_issue_pages_.empty() ) {
		if ( text_exists ) outstream << " ";
		text_exists = true;
		outstream << volume_issue_pages_;
	}
	if ( year_ != 0 ) {
		if ( text_exists ) outstream << " ";
		outstream << "(" << year_ << ").";
		text_exists = true;
	}
	if ( !doi_.empty() ) {
		if ( text_exists ) outstream << "  ";
		outstream << "doi: " << doi_ << ".";
		//text_exists = true; //Uncomment if anything more is added.  This is commented out because the cppcheck test doesn't like it.
	}
}

/// @brief Given a list of lines, parse out a list of author names from a block of the format:
/// [BEGIN_SECTION_NAME]
///      "Given names" "Surname" "Initials"
///      "Given names" "Surname" "Initials"
///      "Given names" "Surname" "Initials"
/// [END_SECTION_NAME]
/// @details Stores the result in destination.
void
Citation::parse_authors(
	std::string const & section_name,
	utility::vector1< AuthorNames > & destination,
	utility::vector1< std::string > const & lines
) {
	bool in_block( false );
	bool found( false );
	std::string const beginstring( "[BEGIN_" + section_name + "]" );
	std::string const endstring( "[END_" + section_name + "]" );
	for ( std::string const & line : lines ) {
		std::string const linestripped( utility::strip(line, " \t\n") );
		if ( !in_block ) {
			if ( linestripped == beginstring ) {
				runtime_assert_string_msg( !found, "Error in Citation::parse_authors():  Multiple " + beginstring + "..." + endstring + " blocks were found in a single citation!" );
				in_block = true;
				found = true;
				continue;
			}
		} else { //In a block
			if ( linestripped == endstring ) {
				in_block = false;
				continue;
			} else {
				utility::vector1< std::string > const author_name( utility::quoted_split( linestripped ) );
				runtime_assert_string_msg( author_name.size() == 3, "Error in Citation::parse_authors():  Could not parse the following:\n" + linestripped + "\nAuthor lines must consist of: \"Given names\" \"Surname\" \"Initial(s)\"." );
				destination.emplace_back( utility::strip( author_name[1], "\"'"), utility::strip( author_name[2], "\"'"), utility::strip( author_name[3], "\"'") );
				continue;
			}
		}
	}
}

/// @brief Given a list of lines, parse out an integer from a block of the format:
/// [BEGIN_YEAR]
///      1982
/// [END_YEAR]
/// @details Stores the result in year_.
void
Citation::parse_year(
	utility::vector1< std::string > const & lines
) {
	bool in_block( false );
	bool found( false ), foundstring(false);
	std::string const beginstring( "[BEGIN_YEAR]" );
	std::string const endstring( "[END_YEAR]" );
	for ( std::string const & line : lines ) {
		std::string const linestripped( utility::strip(line, " \t\n") );
		if ( !in_block ) {
			if ( linestripped == beginstring ) {
				runtime_assert_string_msg( !found, "Error in Citation::parse_year():  Multiple " + beginstring + "..." + endstring + " blocks were found in a single citation!" );
				in_block = true;
				found = true;
				continue;
			}
		} else { //In a block
			if ( linestripped == endstring ) {
				in_block = false;
				continue;
			} else {
				runtime_assert_string_msg( !foundstring, "Error in Citation::parse_year():  Multiple lines were found within a " + beginstring + "..." + endstring + " block!" );
				foundstring = true;
				year_ = std::atoi( linestripped.c_str() );
				continue; //Keep going through all lines to confirm no weird formatting or duplicate blocks.
			}
		}
	}
}

/// @brief Given a list of lines, parse out a string from a block of the format:
/// [BEGIN_SECTION_NAME]
///      The string to parse out
/// [END_SECTION_NAME]
/// @details Stores the result in destination.
void
Citation::parse_string_entry(
	std::string const & section_name,
	std::string & destination,
	utility::vector1< std::string > const & lines
) {
	bool in_block( false );
	bool found( false ), foundstring(false);
	std::string const beginstring( "[BEGIN_" + section_name + "]" );
	std::string const endstring( "[END_" + section_name + "]" );
	for ( std::string const & line : lines ) {
		std::string const linestripped( utility::strip(line, " \t\n") );
		if ( !in_block ) {
			if ( linestripped == beginstring ) {
				runtime_assert_string_msg( !found, "Error in Citation::parse_string_entry():  Multiple " + beginstring + "..." + endstring + " blocks were found in a single citation!" );
				in_block = true;
				found = true;
				continue;
			}
		} else { //In a block
			if ( linestripped == endstring ) {
				in_block = false;
				continue;
			} else {
				runtime_assert_string_msg( !foundstring, "Error in Citation::parse_string_entry():  Multiple lines were found within a " + beginstring + "..." + endstring + " block!" );
				foundstring = true;
				destination = linestripped;
				continue; //Keep going through all lines to confirm no weird formatting or duplicate blocks.
			}
		}
	}
}

} //basic
} //citation_manager
