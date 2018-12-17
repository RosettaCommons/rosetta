// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/mhc_epitope_energy/MHCEpitopeEnergySetupSetup.hh
/// @brief Headers for a helper for the MHCEpitopeEnergy
/// @details Follows analogous file for Vikram K. Mulligan's NetChargeEnergy
/// The helper stores the specification of how to do the epitope scoring for a particular score term or constraint.
/// This is read from a specially formatted ".mhc" file; a default is available in database/scoring/score_function/mhc_epitope/propred8_5.mhc
/// The file specifies what predictor to use and any of its arguments (e.g., propred8 with a 5% threshold for propred, or a precomputed database)
/// It also specifies how to use the computed epitope scores in computing energies. The default is simply to minimize the total.
/// Alternatively, minimization down to an offset value with no returns after that can be used.  This will constrain the score for a peptide not be
/// be much worse than the corresponding wild-type peptide (with options to offset from the exact wild-type value) or an arbitrary threshold.
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu


#ifndef INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopeEnergySetupSetup_hh
#define INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopeEnergySetupSetup_hh

// Unit headers
#include <core/scoring/mhc_epitope_energy/MHCEpitopeEnergySetup.fwd.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictor.hh>

// Package headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

// Project headers
#include <core/types.hh>
#include <map>
#include <string>
#include <utility/vector1.hh>
#include <math.h>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {
namespace mhc_epitope_energy {

/// @brief MHCEpitopeEnergySetup, a helper class for the MHCEpitopeEnergy energy method
/// that stores all of its setup data.
class MHCEpitopeEnergySetup : public utility::pointer::ReferenceCount {
public:

	/// @brief Default constructor for MHCEpitopeEnergySetup.
	///
	MHCEpitopeEnergySetup();

	/// @brief Default destructor for MHCEpitopeEnergySetup.
	///
	virtual ~MHCEpitopeEnergySetup();

	/// @brief Clone: create a copy of this object, and return an owning pointer
	/// to the copy.
	virtual MHCEpitopeEnergySetupOP clone() const;

	virtual
	bool operator == ( MHCEpitopeEnergySetup const & /*other*/ ) const;

public:

	/**********************
	Public setup functions:
	**********************/

	/// @brief Reset all data in this data storage object.
	///
	void reset();

	/// @brief Initialize from a .mhc file.
	///
	void initialize_from_file( std::string const &filename );

	/// @brief Initialize from a string in the format of a .mhc file.
	/// @details Allows external code to initialize object without having it read
	/// directly from disk. Note however that the specs in the file may require
	/// subsequent disk reads depending on the type of predictor indicated.
	void initialize_from_file_contents( std::string const &filecontents );


public:

	/*************************
	Public accessor functions:
	*************************/

	/// @brief Get a summary of the data stored in this object
	///
	std::string report() const;

	/// @brief Is is just the default, always-zero, predictor?
	bool is_default() const { return predictor_ == NULL; }

	/// @brief How long are the peptides to be predicted?
	core::Size get_peptide_length() const { return predictor_->get_peptide_length(); }

	/// @brief The MHC epitope score for the peptide, as returned by the predictor
	core::Real raw_score(std::string peptide) const;

	/// @brief Transform the MHC epitope score, as specified for this helper
	core::Real xform(core::Real /* raw */, core::Real /* native */) const;

	/// @brief The xformed MHC epitope score, in a single step
	core::Real score(std::string peptide, core::Real native) const { return xform(raw_score(peptide), native); }

	/// @brief Getter for predictor_ (MHCEpitopePredictorOP)
	MHCEpitopePredictorOP get_predictor() { return predictor_; }

	/// @brief Return whether using relative or absolute scoring.
	bool get_relative() const { return relative_; }

	/// @brief Return whether to apply an offset to the score.
	bool get_apply_offset() const { return apply_offset_; }

	/// @brief Return whether the relative score is in additive (true) or multiplicative (false) mode
	bool get_relative_additive() const { return relative_additive_; }

	/// @brief Return the score_offset_ value
	core::Real get_score_offset() const { return score_offset_; }

private:

	/******************
	Private functions:
	******************/

	/// @brief Parse the data from an mhc-format file, and create the appropriate MHCEpitopePredictor object. Might trigger a follow-up read from disk to initialize the predictor.
	void parse_specs( utility::vector1 < std::string > const &lines );

private:

	/******************
	Private variables:
	******************/

	MHCEpitopePredictorOP predictor_;

	/// @brief is the contribution computed just from the current score (!relative) or compared to the native score (relative)
	bool relative_ = false;
	/// @brief offset the score?
	bool apply_offset_ = false;
	/// @brief if relative and offset, is the offset a value added to the native score (additive) or a factor times the native score (multiplicative / !additive)
	bool relative_additive_ = true;
	/// @brief if apply_offset_, then scores below this value contribute 0; above this value, linearly
	core::Real score_offset_ = 0.0;


#ifdef    SERIALIZATION
protected:
	friend class cereal::access;

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // mhc_epitope_energy
} // scoring
} // core

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_scoring_mhc_epitope_energy_MHCEpitopeEnergySetup )
#endif // SERIALIZATION

#endif // INCLUDED_core_scoring_EtableEnergy_HH
