// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/citation_manager/CitationManager
/// @brief A class that receives lists of works to cite from Rosetta modules, then returns a list of all
/// works to cite on demand.  Threadsafe.
/// @details For works with publications that have DOIs, this loads a list of references from the database,
/// indexed by DOI, once in a threadsafe manner on object creation.  This allows modules to only specify
/// the DOI.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_basic_citation_manager_CitationManager_hh
#define INCLUDED_basic_citation_manager_CitationManager_hh

// Unit headers
#include <basic/citation_manager/CitationManager.fwd.hh>
#include <basic/citation_manager/CitationCollection.fwd.hh>
#include <basic/citation_manager/Citation.fwd.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.fwd.hh>

// Utility header
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>

// Multithreading headers
#ifdef MULTI_THREADED
#include <mutex>
#endif

namespace basic {
namespace citation_manager {

/// @brief A class that receives lists of works to cite from Rosetta modules, then returns a list of all works to cite on demand.  Threadsafe.
class CitationManager : public utility::SingletonBase< CitationManager > {
	friend class utility::SingletonBase< CitationManager >;

private:  // Private methods //////////////////////////////////////////////////
	/// @brief Constructor triggers read from disk.
	CitationManager();
	CitationManager(CitationManager const & ) = delete;
	CitationManager operator=(CitationManager const & ) = delete;

public:  // Accessors //////////////////////////////////////////////////////////

	/// @brief Clear all citations that have been collected.
	/// @details Does not clear the map of DOI->Rosetta citation that was loaded from
	/// the Rosetta database.
	/// @details Threadsafe.
	void clear_citations();

	/// @brief Has the CitationModule accumulated at least one module to cite?
	/// @details Threadsafe.
	bool has_modules_to_cite() const;

	/// @brief Has the CitationModule accumulated at least one unpublished module?
	/// @details Threadsafe.
	bool has_unpublished_modules() const;

	/// @brief Add citations to the list of citations that have been collected.
	/// @note Checks module names and types, and only adds citations that have not
	/// been already added.
	/// @details Threadsafe.
	void add_citations( utility::vector1< CitationCollectionCOP > const & input );

	/// @brief Add information for one or more unpublished modules to the list that has been collected.
	/// @note Checks that these modules haven't been added, and only adds new ones.
	/// @details Threadsafe.
	void add_unpublished_modules( utility::vector1< UnpublishedModuleInfoCOP > const & input );

	/// @brief Get a summary of all the citations that we've collected so far.
	/// @note This is ONLY a list of citations, in a given format, on separate lines, with
	/// header lines indicating the relevant module.  This does NOT include an overall header
	/// explaining the context in which citations were collected.
	/// @details Threadsafe.
	void write_collected_citations( std::ostream & outstream, CitationFormat const citation_format = CitationFormat::DefaultStyle ) const;

	/// @brief Write out a list of the unpublished modules.
	void write_unpublished_modules( std::ostream & outstream ) const;

	/// @brief Write out all unpublished modules and citations to the CitationManager's tracer.
	void write_all_citations_and_unpublished_author_info() const;

	/// @brief Given a DOI string, get a Rosetta citation.
	/// @details Throws if the DOI string isn't in the list of Rosetta papers in the database.
	CitationCOP get_citation_by_doi( std::string const & doi ) const;

private: // Private fxns /////////////////////////////////////////////////////

	/// @brief Load Rosetta citations from the Rosetta database, and populate
	/// the doi_rosetta_citation_map_.
	/// @details TRIGGERS READ FROM DISK.  Threadsafe, but should only be done
	/// once!
	void load_rosetta_citations_from_database();

	/// @brief Populate the doi_rosetta_citation_map_ from the contents of a
	/// database file.
	void populate_doi_rosetta_citation_map( std::string const & database_file_contents );

private: // Private data /////////////////////////////////////////////////////

#ifdef MULTI_THREADED
	/// @brief A mutex, for thread-safety.
	mutable std::mutex mutex_;
#endif

	/// @brief A map of DOI->Rosetta citations.
	std::map< std::string, CitationCOP > doi_rosetta_citation_map_;

	/// @brief A list of everything that has been cited so far.
	utility::vector1< CitationCollectionCOP > citation_collections_;

	/// @brief A list of the modules that are unpublished that have been used so far.
	utility::vector1< UnpublishedModuleInfoCOP > unpublished_;

};

} //citation_manager
} //basic

#endif //INCLUDED_basic/citation_manager_CitationManager_fwd_hh



