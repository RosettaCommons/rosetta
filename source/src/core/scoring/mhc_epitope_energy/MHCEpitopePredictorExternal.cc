// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/mhc_epitope_energy/MHCEpitopePredictorExternal.cc
/// @brief MHC epitope predictor using calculations from an external program cached in an sqlite database
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

#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorExternal.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace mhc_epitope_energy {

using basic::t_info;
using basic::t_debug;
using basic::t_trace;
static basic::Tracer TR("core.scoring.mhc_epitope_energy.MHCEpitopePredictorExternal");

using namespace utility::sql_database;
using namespace cppdb;

/// @brief Initializes with a connection to an sqlite database
MHCEpitopePredictorExternal::MHCEpitopePredictorExternal(std::string const &db_filename)
{
	connect(db_filename);
}

bool MHCEpitopePredictorExternal::operator==(MHCEpitopePredictor const &other)
{
	MHCEpitopePredictorExternal const *o = dynamic_cast<MHCEpitopePredictorExternal const *>(&other);
	if ( !o ) return false;

	if ( o->filename_ != filename_ ) return false; // could be overly strict -- different filenames could have same data, but let's not go there
	if ( o->unseen_score_ != unseen_score_ ) return false;

	return true;
}

// Basic SQL stuff following test/utility/sql_database/DatabaseSessionManagerTests.cxxtest.hh

std::string MHCEpitopePredictorExternal::report() const
{
	std::stringstream output("");

	// TODO: Any other metadata? If so, hold on to it in connect()
	output << "External predictor database from file '" << filename_ << "' based on " << pred_name_ << "; peptide length " << get_peptide_length() << "; unseen score " << unseen_score_;

	return output.str();
}

core::Real MHCEpitopePredictorExternal::score(std::string const &pep)
{
	if ( pep.size() != get_peptide_length() ) {
		TR.Error << "Scoring peptide of size " << pep.size() << " with an external database expecting peptides of size " << get_peptide_length() << std::endl;
		utility_exit_with_message("MHCEpitopePredicotrExternal is trying to score a peptide of the incorrect size!");
	}

	// Simply gets the score from the database, returning unseen_score_ if it isn't there

	cppdb::result fetch = ((*session_) << "SELECT score FROM epitopes WHERE peptide=?" << pep).row(); // fetch single row if it exists
	if ( fetch.empty() ) { // if the peptide doesn't exist in the database, return the unseen_score_
		return unseen_score_;
	}
	core::Real score;
	fetch >> score;
	return score;
}

void MHCEpitopePredictorExternal::connect(std::string const & filename)
{
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
	filename_ = db_filename;

	try {
		// Establish the connection
		DatabaseSessionManager * scm( DatabaseSessionManager::get_instance() );

		session_ = sessionOP( scm->get_session_sqlite3( db_filename, utility::sql_database::TransactionMode::standard, 0, true, -1 ) ); // readonly set to true here, everything else is default settings.

		// Fetch the metadata
		// TODO: anything else important here?
		cppdb::result meta = (*session_) << "SELECT name,value from meta";
		while ( meta.next() ) {
			std::string name, val;
			meta >> name >> val;
			TR.Debug << name << ":" << val << std::endl;
			if ( name == "peptide_length" ) set_peptide_length(utility::string2int(val));
			else if ( name == "predictor" ) pred_name_ = val;
		}

		// Make sure we know how long peptides are supposed to be
		if ( get_peptide_length() == 0 ) utility_exit_with_message("Database didn't specify peptide length");
	}
catch (std::exception const &e) {
	utility_exit_with_message("Unable to open valid database " + db_filename + ": " + e.what());
}
}

}//ns mhc_epitope_energy
}//ns scoring
}//ns core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::mhc_epitope_energy::MHCEpitopePredictorExternal::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::mhc_epitope_energy::MHCEpitopePredictor >( this ) );
	// Don't serialize session_; instead re-open the session during load()
	// EXEMPT session_
	arc( CEREAL_NVP( filename_ ) ); // std::string
	arc( CEREAL_NVP( pred_name_ ) ); // std::string
	arc( CEREAL_NVP( unseen_score_ ) ); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::mhc_epitope_energy::MHCEpitopePredictorExternal::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::mhc_epitope_energy::MHCEpitopePredictor >( this ) );
	// Don't serialize session_; instead re-open the session during load()
	// EXEMPT session_
	arc( filename_ ); // std::string
	arc( pred_name_ ); // std::string
	arc( unseen_score_ ); // core::Size

	// Re-connect to the session
	connect( filename_ );
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::mhc_epitope_energy::MHCEpitopePredictorExternal );
CEREAL_REGISTER_TYPE( core::scoring::mhc_epitope_energy::MHCEpitopePredictorExternal )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_mhc_epitope_energy_MHCEpitopePredictorExternal )
#endif // SERIALIZATION
