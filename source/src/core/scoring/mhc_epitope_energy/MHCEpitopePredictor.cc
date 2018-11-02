
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.

/// @file   core/scoring/mhc_epitope_energy/MHCEpitopePredictor.cc
/// @brief Base class for an MHC epitope predictor, which takes a peptide (string) and returns a score predicting its risk of MHC binding (lower is lower risk, with 0 being none)
/// @details The peptides are of fixed length, specified for the predictor (e.g., for class II MHC: 9 for Propred, 15 for NetMHCII and IEDB)
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

// Unit headers
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictor.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::mhc_epitope_energy::MHCEpitopePredictor::save( Archive & arc ) const {
	arc( CEREAL_NVP( peptide_length_ ) ); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::mhc_epitope_energy::MHCEpitopePredictor::load( Archive & arc ) {
	arc( peptide_length_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::mhc_epitope_energy::MHCEpitopePredictor );
CEREAL_REGISTER_TYPE( core::scoring::mhc_epitope_energy::MHCEpitopePredictor )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_mhc_epitope_energy_MHCEpitopePredictor )
#endif // SERIALIZATION
