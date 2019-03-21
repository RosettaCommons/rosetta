// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/mhc_epitope_energy/MHCEpitopePredictorExternal.hh
/// @brief MHC epitope predictor using calculations from an external program cached in an sqlite database
/// @details The peptides that might be encountered during packing should have scores stored in the "epitopes" table, with column "peptide" given the string and "score" the score.
/// For peptides not stored in the database, the current default method is to simply return the value "unseen_score" (positive is a penalty). The class could be augmented to invoke an external program on the fly and cache the computed result, but since that is slow, we started with the method of precomputing and caching scores of desireable peptides. See the tools directory for scripts to help with that; additional constraints (via e.g., FavorSequenceProfile) can be used to focus the amino acid selections on desirable choices, as demonstrated in the examples there.
/// The database should also include a "meta" table with information about the predictor stored as a map with columns "name" and "value". Currently, the only thing that's needed from there is "peptide_length" (for MHCEpitopePredictor) and the specified length.
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

#ifndef INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopePredictorExternal_hh
#define INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopePredictorExternal_hh

#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictor.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorExternal.fwd.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/vector1.hh>

#include <utility/sql_database/DatabaseSessionManager.hh>

#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace mhc_epitope_energy {

class MHCEpitopePredictorExternal: public MHCEpitopePredictor {

public:
	MHCEpitopePredictorExternal() {}

	/// @brief Initializes with a connection to an sqlite database
	MHCEpitopePredictorExternal(std::string const &db_filename);

	virtual ~MHCEpitopePredictorExternal() {}

	virtual bool operator==(MHCEpitopePredictor const & /* other */);

	std::string report() const;

	virtual core::Real score(std::string const &pep);

	/// @brief Establishes the connection to the sqlite database and gets the meta information about the predictor
	void connect(std::string const &db_filename);

	/// @brief Return the sqlite database connection pointer
	utility::sql_database::sessionOP get_database() { return session_; }

	/// @brief Sets the score for a peptide not in the database
	void set_unseen_score(core::Size u) { unseen_score_ = u; }

	/// @brief Accessor for the unseen_score_
	core::Size get_unseen_score() { return unseen_score_; }

private:
	/// @brief The name of the database filename.
	std::string filename_="";

	/// @brief The connection to the sqlite database, initialized via connect()
	utility::sql_database::sessionOP session_ = nullptr;

	/// @brief The underlying predictor that generated the database, read from database
	std::string pred_name_ = "";

	/// @brief The score for a peptide not in the database
	core::Size unseen_score_ = 100;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class MHCEpitopePredictorExternal

} // mhc_epitope_energy
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_mhc_epitope_energy_MHCEpitopePredictorExternal )
#endif // SERIALIZATION


#endif
