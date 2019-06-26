// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/mhc_epitope_energy/MHCEpitopePredictorExternal.hh
/// @brief MHC epitope predictor with precomputed values loaded into a peptide->score map
/// @details Different initializers load the map from different format files
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

#ifndef INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopePredictorPreLoaded_hh
#define INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopePredictorPreLoaded_hh

#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictor.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorPreLoaded.fwd.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/vector1.hh>
#include <map>

#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace mhc_epitope_energy {

class MHCEpitopePredictorPreLoaded: public MHCEpitopePredictor {

public:
	MHCEpitopePredictorPreLoaded() {}

	/// @brief Initializes from an sqlite database
	void load_database(std::string const &filename);

	/// @brief Initializes from a csv file with header, such that the first column is the peptide and the second column the score; other columns ignored
	void load_csv(std::string const &filename);

	virtual ~MHCEpitopePredictorPreLoaded() {}

	virtual bool operator==(MHCEpitopePredictor const & /* other */);

	std::string report() const;

	virtual core::Real score(std::string const &pep);

	/// @brief Sets the score for a peptide not in the map
	void set_unseen_score(core::Size us) { unseen_score_ = us; }

	/// @brief Accessor for the unseen_score_
	core::Size get_unseen_score() const { return unseen_score_; }

	/// @brief Accessor for the peptide->score map object
	std::map< std::string, core::Real > get_scores_map() const { return scores_; }

private:
	/// @brief Check for the size of a file and print a warning if appropriate.
	void check_file_size ( std::string const &filename, core::Size warn_threshold ) const;

private:
	/// @brief peptide->score
	std::map< std::string, core::Real > scores_;

	/// @brief The filename from which the map was loaded.
	std::string filename_="";

	/// @brief The type of file, so can reconstruct when needed (serialization)
	enum FileType { NOT_LOADED, LOAD_DATABASE, LOAD_CSV };
	FileType filetype_ = NOT_LOADED;

	/// @brief The score for a peptide not in the map
	core::Size unseen_score_ = 100;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class MHCEpitopePredictorPreLoaded

} // mhc_epitope_energy
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_mhc_epitope_energy_MHCEpitopePredictorPreLoaded )
#endif // SERIALIZATION


#endif
