// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/mhc_epitope_energy/MHCEpitopePredictorMatrix.hh
/// @brief MHC epitope predictor using a position weight matrix, targeted to Propred though in theory generalizable to others
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

#ifndef INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopePredictorMatrix_hh
#define INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopePredictorMatrix_hh

#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictor.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorMatrix.fwd.hh>
#include <core/chemical/AA.hh>
#include <utility/file/FileName.fwd.hh>
#include <utility/vector1.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <map>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace mhc_epitope_energy {

/// @brief The scoring matrix for one specific allele
/// @details A position weight matrix (PWM) specifies the value for each amino acid type. A set of threshold specify the minimum score required to be in the top percentile of predicted binders.
class AlleleMatrix: public utility::pointer::ReferenceCount {
public:
	typedef std::map< char, core::Real > Weights;
	typedef utility::vector1< Weights > PWM;

	AlleleMatrix();
	AlleleMatrix(std::string name, utility::vector1< Real > threshes, PWM profile);
	~AlleleMatrix();

	bool operator==(AlleleMatrix const & /* other */);

	/// @brief Predicts whether the peptide will bind the MHC allele, with respect to the threshold for binding strength / likelihood
	bool is_hit(std::string const &pep, Real thresh);

private:
	/// @brief Arbitrary name for the allele
	std::string name_;

	/// @brief Hit thresholds: element i is the minimum score to be in the predicted top i% of binders
	utility::vector1< Real > threshes_;

	/// @brief The weights: a vector over peptide positions of a map from character (amino acid) to weight
	PWM profile_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION
};


/// @brief An predictor based on a set of allele-specific profiles of which peptides are likely to be binders
/// @details Each allele's binding is predicted separately, thresholded with respect to the top % of binders, and the score is then the number of predicted binding events
class MHCEpitopePredictorMatrix: public MHCEpitopePredictor {

public:
	MHCEpitopePredictorMatrix();
	/// @brief Loads the matrices from the file
	MHCEpitopePredictorMatrix(std::string const &fn);
	virtual ~MHCEpitopePredictorMatrix();

	virtual bool operator==(MHCEpitopePredictor const & /* other */);

	std::string report() const;

	virtual core::Real score( std::string const &pep);

	/// @brief Loads the matrices from the file
	void load_matrix(std::string const &filename);

	/// @brief Sets the threshold for what is considered to be an epitope -- top thresh% of peptides in this implementation
	void set_thresh(core::Real thresh);

	/// @brief Gets the threshold for what is considered to be an epitope -- top thresh% of peptides in this implementation
	core::Real get_thresh() { return thresh_; }

private:
	/// @brief The name of the filename containing the matrices.
	std::string filename_="";

	/// @brief The matrices for the alleles
	utility::vector1< AlleleMatrix > alleles_;

	/// @brief The threshold for what is considered to be an epitope -- top thresh% of peptides in this implementation
	core::Real thresh_ = 5;

	/// @brief A flag to indicate that the matrix is a propred matrix
	bool propred_ = false;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

}; // class MHCEpitopePredictorMatrix

} // mhc_epitope_energy
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_mhc_epitope_energy_MHCEpitopePredictorMatrix )
#endif // SERIALIZATION


#endif
