// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/mhc_epitope_energy/MHCEpitopePredictorPreLoaded.cc
/// @brief MHC epitope predictor with precomputed values loaded into a peptide->score map
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <basic/database/open.hh>
#include <iostream>
#include <string>
#include <utility/string_util.hh>
#include <cppdb/frontend.h>
#include <utility/sql_database/DatabaseSessionManager.hh>

#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorPreLoaded.hh>
#include <sys/stat.h> // For getting file size with stat().  If we move to C++17, should switch to std::filesystem::file_size() function

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/map.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace mhc_epitope_energy {

using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static basic::Tracer TR("core.scoring.mhc_epitope_energy.MHCEpitopePredictorPreLoaded");

using namespace utility::sql_database;
using namespace cppdb;

/// @brief Check for the size of a file and print a warning if appropriate.
void MHCEpitopePredictorPreLoaded::check_file_size( std::string const &filename, core::Size warn_threshold ) const {
	core::Real filesize = 0; // Will be used to store the file's size, in bytes.

	// If we switch to C++17, re-write from here to replace stat() function with std::filesystem::file_size()
	// Get rid of the c-string and the platform-specific stuff.
	// Need to get filename as a C-string for use with stat(): filecstring
	char *filecstring = new char[ filename.length() + 1 ];
	std::strcpy( filecstring, filename.c_str() );

	// The following is using the stat() function to get the file's size.  If we switch to C++17, it should be switched to std::filesystem::file_size() to make it platform agnostic.
#ifdef WIN32
	struct _stat file_info; // stat object contains all file info.
	_stat( filecstring, &file_info ); // Get filename's info and store in file_info buffer object.
#else
	struct stat file_info; // stat object contains all file info.
	stat( filecstring, &file_info ); // Get filename's info and store in file_info buffer object.
#endif

	delete[] filecstring; // Deallocate the c-string.

	filesize = file_info.st_size; // Put the filesize (in bytes) into filesize.

	// End of C++17 replacement section.

	if ( filesize > warn_threshold ) {
		TR.Warning << "You are attempting to load " << filename << " into memory as a PreLoaded mhc_epitope database." << std::endl;
		TR.Warning << "This file is " << filesize / (1024*1024*1024) << "GB." << std::endl;
		TR.Warning << "Ensure that you have enough memory to support it, or use an External database stored on disk." << std::endl;
		TR.flush();
	}
}

/// @brief Initializes from an sqlite database
void MHCEpitopePredictorPreLoaded::load_database(std::string const &filename)
{
	filetype_ = LOAD_DATABASE;

	std::string db_filename = ""; //Actual location of the db filename
	// Look for a copy of the file
	utility::io::izstream infile;
	infile.open( filename );
	// If the filename is found, store in db_filename.  Otherwise, look in the database.
	if ( infile.good() ) {
		db_filename = filename;
	} else {
		db_filename = basic::database::full_name("scoring/score_functions/mhc_epitope/"+filename);
	}
	TR << "Connecting to " << db_filename << " external database." << std::endl;
	// Output a warning if db_filename is bigger than 1 GB (1073741824 bytes).
	check_file_size( db_filename, 1073741824 );

	filename_ = db_filename;

	// Basic SQL stuff following test/utility/sql_database/DatabaseSessionManagerTests.cxxtest.hh

	try {
		// Establish the connection
		DatabaseSessionManager * scm( DatabaseSessionManager::get_instance() );

		utility::sql_database::sessionOP session(scm->get_session_sqlite3( db_filename, utility::sql_database::TransactionMode::standard, 0, true, -1 )); // readonly set to true here, everything else is default settings.

		// Fetch the metadata
		// TODO: anything else important here?
		cppdb::result meta = (*session) << "SELECT name,value from meta";
		while ( meta.next() ) {
			std::string name, val;
			meta >> name >> val;
			TR.Debug << name << ":" << val << std::endl;
			if ( name == "peptide_length" ) set_peptide_length(utility::string2int(val));
		}

		// Make sure we know how long peptides are supposed to be
		if ( get_peptide_length() == 0 ) utility_exit_with_message("Database didn't specify peptide length");

		// Now load the map
		cppdb::result pep_scores = (*session) << "SELECT peptide,score from epitopes";
		while ( pep_scores.next() ) {
			std::string peptide;
			core::Real score;
			pep_scores >> peptide >> score;
			TR.Debug << peptide << ":" << score << std::endl;
			scores_[peptide] = score;
		}
	}
catch (std::exception const &e) {
	utility_exit_with_message("Unable to open valid database " + db_filename + ": " + e.what());
}
}

void MHCEpitopePredictorPreLoaded::load_csv(std::string const &filename)
{
	using namespace utility::io;

	filetype_ = LOAD_CSV;

	std::string db_filename = ""; //Actual location of the db filename
	// Look for a copy of the file
	izstream infile;
	infile.open( filename );
	// If the filename is found, store in db_filename.  Otherwise, look in the database.
	if ( infile.good() ) {
		db_filename = filename;
	} else {
		db_filename = basic::database::full_name("scoring/score_functions/mhc_epitope/"+filename);
	}
	infile.close();
	TR << "Reading epitopes from file " << db_filename << std::endl;
	// Output a warning if filename is bigger than 1 GB (1073741824 bytes).
	check_file_size( db_filename, 1073741824 );

	filename_ = db_filename;

	// TODO: very simple csv parsing (but that's presumably all we need). if is there a robust library, wouldn't hurt to use it instead
	//izstream infile;
	infile.open( db_filename );
	if ( !infile.good() ) utility_exit_with_message("ERROR: Unable to open file " + db_filename);

	std::string curline(""); //Buffer for current line.
	bool got_header = false;
	while ( getline(infile, curline) ) {
		if ( curline.size() < 1 || curline[0]=='#' ) continue; //Ignore blank lines and comments
		// First real row should be a header; ignore it too
		if ( ! got_header ) {
			got_header = true;
			continue;
		}
		utility::vector1<std::string> cols = utility::string_split(curline, ',');
		if ( cols.size() < 2 ) {
			utility_exit_with_message("ERROR: bad row " + curline + " -- not enough columns");
		}
		std::string peptide = cols[1];
		if ( get_peptide_length() == 0 ) {
			// determine peptide length from first peptide
			set_peptide_length(peptide.size());
		} else if ( peptide.size() != get_peptide_length() ) {
			// check to make sure they're all of same length
			utility_exit_with_message("all peptides should be of the same length");
		}

		core::Real score = utility::string2Real(cols[2]);
		// TODO: check for / handle duplicates? maybe have an option to warn/error
		scores_[peptide] = score;
	}
	infile.close();

	if ( get_peptide_length() == 0 ) {
		utility_exit_with_message("ERROR: no peptides in the file");
		// if epitope length remains 0, bad things happen in extracting and scoring peptides, so bail
	}
}

bool MHCEpitopePredictorPreLoaded::operator==(MHCEpitopePredictor const &other)
{
	MHCEpitopePredictorPreLoaded const *o = dynamic_cast<MHCEpitopePredictorPreLoaded const *>(&other);
	if ( !o ) return false;

	if ( o->filename_ != filename_ ) return false; // could be overly strict -- different filenames could have same data, but let's not go there
	if ( o->unseen_score_ != unseen_score_ ) return false;

	return true;
}

std::string MHCEpitopePredictorPreLoaded::report() const
{
	std::stringstream output("");

	// TODO: Any other metadata? If so, hold on to it in connect()
	output << "Map predictor from file '" << filename_ << " containing " << scores_.size() << " peptides; peptide length " << get_peptide_length() << "; unseen score " << unseen_score_;

	return output.str();
}

core::Real MHCEpitopePredictorPreLoaded::score(std::string const &pep)
{
	if ( pep.size() != get_peptide_length() ) {
		TR.Error << "Scoring peptide of size " << pep.size() << " with map expecting peptides of size " << get_peptide_length() << std::endl;
		utility_exit_with_message("MHCEpitopePredicotrMap is trying to score a peptide of the incorrect size!");
	}

	if ( scores_.count(pep)==0 ) return unseen_score_;
	return scores_[pep];
}

}//ns mhc_epitope_energy
}//ns scoring
}//ns core

#ifdef    SERIALIZATION

template< class Archive >
void
core::scoring::mhc_epitope_energy::MHCEpitopePredictorPreLoaded::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::mhc_epitope_energy::MHCEpitopePredictor >( this ) );
	arc( CEREAL_NVP( scores_ ) ); // std::map< std::string, core::Real >
	arc( CEREAL_NVP( filename_ ) ); // std::string
	arc( CEREAL_NVP( filetype_ ) ); // enum
	arc( CEREAL_NVP( unseen_score_ ) ); // core::Size
}

template< class Archive >
void
core::scoring::mhc_epitope_energy::MHCEpitopePredictorPreLoaded::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::mhc_epitope_energy::MHCEpitopePredictor >( this ) );
	arc( scores_ ); // std::map< std::string, core::Real >
	arc( filename_ ); // std::string
	arc( filetype_ ); // enum
	arc( unseen_score_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::mhc_epitope_energy::MHCEpitopePredictorPreLoaded );
CEREAL_REGISTER_TYPE( core::scoring::mhc_epitope_energy::MHCEpitopePredictorPreLoaded )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_mhc_epitope_energy_MHCEpitopePredictorPreLoaded )
#endif // SERIALIZATION
