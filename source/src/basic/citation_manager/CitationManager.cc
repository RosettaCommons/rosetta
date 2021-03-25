// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/CitationManager.cc
/// @brief A class that receives lists of works to cite from Rosetta modules, then returns a list of all
/// works to cite on demand.  Threadsafe.
/// @details For works with publications that have DOIs, this loads a list of references from the database,
/// indexed by DOI, once in a threadsafe manner on object creation.  This allows modules to only specify
/// the DOI.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <basic/citation_manager/CitationManager.hh>
#include <basic/citation_manager/CitationCollection.hh>
#include <basic/citation_manager/Citation.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>


// Unit headers

// Project header

// Utility headers
#include <utility/string_util.hh>
#include <utility/pointer/memory.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// C++ headers
#include <string>

// Construct tracer.
static basic::Tracer TR( "basic.citation_manager.CitationManager" );

namespace basic {
namespace citation_manager {

// Public methods /////////////////////////////////////////////////////////////
// Static constant data access

// Private methods ////////////////////////////////////////////////////////////
/// @brief Constructor triggers read from disk.
CitationManager::CitationManager()
{
	load_rosetta_citations_from_database();
}


/// @brief Clear all citations that have been collected.
/// @details Threadsafe.
void
CitationManager::clear_citations() {
#ifdef MULTI_THREADED
	std::lock_guard< std::mutex > guard( mutex_ );
#endif
	citation_list_ = CitationCollectionList{};
}

/// @brief Add citations to the list of citations that have been collected.
/// @note Only adds citations that have not been already added.
/// @details Threadsafe.
void
CitationManager::add_citations( CitationCollectionList const & input ) {
	if ( input.empty() ) return;
#ifdef MULTI_THREADED
	std::lock_guard< std::mutex > guard( mutex_ );
#endif
	citation_list_.add( input );
}

/// @brief Add a single citation to the citation manger
/// @note Will not add if a duplicate
/// @details Threadsafe
void
CitationManager::add_citation( CitationCollectionBaseCOP const & input ) {
	if ( input == nullptr ) return;
#ifdef MULTI_THREADED
	std::lock_guard< std::mutex > guard( mutex_ );
#endif
	citation_list_.add( input );
}

/// @brief Get a summary of all the citations that we've collected so far.
/// @note This is ONLY a list of citations, in a given format, on separate lines, with
/// header lines indicating the relevant module.  This does NOT include an overall header
/// explaining the context in which citations were collected.
/// @details Threadsafe.
void
CitationManager::write_collected_citations(
	std::ostream & outstream,
	utility::vector1< CitationCollectionCOP > const & published,
	CitationFormat const citation_format
) const {
	for ( auto const & collection : published ) {
		outstream << collection->module_name() << " " << collection->module_type() << "'s citation(s):\n";
		collection->get_citations_formatted(outstream, citation_format);
		outstream << "\n";
	}
}

/// @brief Write out a list of the unpublished modules.
void
CitationManager::write_unpublished_modules(
	std::ostream & outstream,
	utility::vector1< UnpublishedModuleInfoCOP > const & unpublished
) const {
	for ( auto const & module : unpublished ) {
		outstream << module->module_name() << " " << module->module_type() << "'s author(s):\n";
		module->get_author_list( outstream );
		outstream << "\n";
	}
}

/// @brief Write out all unpublished modules and citations to to the CitationManager's tracer.
void
CitationManager::write_all_citations_and_unpublished_author_info() const {
	bool anything_written(false);

	utility::vector1< CitationCollectionCOP > published;
	utility::vector1< UnpublishedModuleInfoCOP > unpublished;
	split_citations( published, unpublished );

	if ( ! published.empty() ) {
		TR << TR.bgBlue << TR.White << TR.Bold << TR.Underline;
		TR << "\nThe following Rosetta modules were used during this run of Rosetta, and should be cited:\n\n";
		TR << TR.Reset << TR.bgBlue << TR.White;
		write_collected_citations( TR, published );
		anything_written = true;
	}

	if ( ! unpublished.empty() ) {
		TR << TR.Reset << TR.bgBlue <<  TR.White << TR.Bold << TR.Underline;
		TR << "\nThe following UNPUBLISHED Rosetta modules were used during this run of Rosetta.  Their authors should be included in the author list when this work is published:\n\n";
		TR << TR.Reset << TR.bgBlue << TR.White;
		write_unpublished_modules( TR, unpublished );
		anything_written = true;
	}
	if ( anything_written ) {
		TR << TR.Reset << std::endl;
	}
}

/// @brief Write out all unpublished modules and citations in a CitationCollectionList to a given output stream.
void
CitationManager::write_all_citations_and_unpublished_author_info_from_list_to_stream(
	CitationCollectionList const & list,
	std::ostream & outstream
) const {
	utility::vector1< CitationCollectionCOP > published;
	utility::vector1< UnpublishedModuleInfoCOP > unpublished;
	split_citations( list, published, unpublished );

	if ( ! published.empty() ) {
		write_collected_citations( outstream, published );
	}

	if ( ! unpublished.empty() ) {
		write_unpublished_modules( outstream, unpublished );
	}
}

/// @brief Given a DOI string, get a Rosetta citation.
/// @details Throws if the DOI string isn't in the list of Rosetta papers in the database.
CitationCOP
CitationManager::get_citation_by_doi(
	std::string const & doi
) const {
#ifdef MULTI_THREADED
	std::lock_guard< std::mutex > guard( mutex_ );
#endif
	try {
		return doi_rosetta_citation_map_.at( doi );
	} catch( std::out_of_range const & err) {
		utility_exit_with_message( "Error in CitationManager::get_citation_by_doi(): The DOI \"" + doi + "\" was requested, but this is not associated with any paper in the Rosetta database.  Has it been added to database/citations/rosetta_citations.txt?" );
	}
}

///////////// PRIVATE FUNCTIONS /////////////

/// @brief Load Rosetta citations from the Rosetta database, and populate
/// the doi_rosetta_citation_map_.
/// @details TRIGGERS READ FROM DISK.  Threadsafe, but should only be done
/// once!
void
CitationManager::load_rosetta_citations_from_database() {
	std::string const filename( basic::database::full_name( "citations/rosetta_citations.txt" ) );
	std::string const file_contents( utility::file_contents( filename ) );
	debug_assert( !file_contents.empty() );
	TR.Debug << "Read Rosetta citations from " << filename << "." << std::endl;
	populate_doi_rosetta_citation_map( file_contents );
}

/// @brief Populate the doi_rosetta_citation_map_ from the contents of a
/// database file.
void
CitationManager::populate_doi_rosetta_citation_map(
	std::string const & database_file_contents
) {
#ifdef MULTI_THREADED
	std::lock_guard< std::mutex > guard( mutex_ );
#endif
	std::string const errmsg( "Error in CitationManager::populate_doi_rosetta_citation_map():  " );
	runtime_assert_string_msg( doi_rosetta_citation_map_.empty(), errmsg + "This function cannot be called more than once!" );

	bool in_block(false);
	std::stringstream accumulator;
	utility::vector1< std::string > const lines( utility::split_by_newlines( database_file_contents ) );
	for ( std::string const & line : lines ) {
		std::string const linestripped( utility::strip( line, " \t") );
		if ( linestripped[0] == '#' ) continue; //Skip comment lines
		if ( !in_block && linestripped == "[BEGIN_CITATION]" ) {
			in_block = true;
			accumulator.clear(); //Only clears error state.
			accumulator.str(""); //Actually clears the string.
			accumulator << linestripped << "\n";
		} else { //If in block...
			accumulator << linestripped << "\n";
			if ( linestripped == "[END_CITATION]" ) {
				in_block = false;
				CitationOP citation( utility::pointer::make_shared< Citation >( accumulator.str() ) );
				debug_assert( !citation->doi().empty() );
				runtime_assert_string_msg( doi_rosetta_citation_map_.count( citation->doi() ) == 0, errmsg + "The DOI \"" + citation->doi() + "\" appears more than once in the list of Rosetta citations." );
				doi_rosetta_citation_map_[citation->doi()] = citation;
				accumulator.clear(); //Only clears error state.
				accumulator.str(""); //Actually clears the string.
			}
		}
	}
	TR.Debug << "Parsed " << doi_rosetta_citation_map_.size() << " Rosetta citations from database file contents." << std::endl;
}

/// @brief Split the citations into published & unpublished vectors.
/// @details Return by reference.
/// Not threadsafe; list must be protected by a lock guard before calling this if there is the
/// possibility that another thread could be acting on it.
void
CitationManager::split_citations(
	CitationCollectionList const & list,
	utility::vector1< CitationCollectionCOP > & published,
	utility::vector1< UnpublishedModuleInfoCOP > & unpublished
) const {
	if ( list.empty() ) { return; }
	for ( CitationCollectionBaseCOP const & entry : list.citations() ) {
		CitationCollectionCOP cc = utility::pointer::dynamic_pointer_cast< CitationCollection const >( entry );
		if ( cc ) {
			published.push_back( cc );
			continue;
		}
		UnpublishedModuleInfoCOP um = utility::pointer::dynamic_pointer_cast< UnpublishedModuleInfo const >( entry );
		if ( um ) {
			unpublished.push_back( um );
			continue;
		}
		//TR.Error << "Unrecognized citation collection type. Typeid: " << typeid(*entry).name() << std::endl;
		utility_exit_with_message("Error in CitationManager::split_citations():  Unable to determine the type of a citation collection.");
	}
}

/// @brief Split the citations into published & unpublished vectors
/// @details Return by reference.
void
CitationManager::split_citations(
	utility::vector1< CitationCollectionCOP > & published,
	utility::vector1< UnpublishedModuleInfoCOP > & unpublished
) const {
#ifdef MULTI_THREADED
	std::lock_guard< std::mutex > guard( mutex_ );
#endif
	if ( citation_list_.empty() ) { return; }
	for ( CitationCollectionBaseCOP const & entry : citation_list_.citations() ) {
		CitationCollectionCOP cc = utility::pointer::dynamic_pointer_cast< CitationCollection const >( entry );
		if ( cc ) {
			published.push_back( cc );
			continue;
		}
		UnpublishedModuleInfoCOP um = utility::pointer::dynamic_pointer_cast< UnpublishedModuleInfo const >( entry );
		if ( um ) {
			unpublished.push_back( um );
			continue;
		}
		//TR.Error << "Unrecognized citation collection type. Typeid: " << typeid(*entry).name() << std::endl;
		utility_exit_with_message("CitationManager::split_citations():  Unable to determine the type of a citation collection.");
	}
}

} //citation_manager
} //basic
