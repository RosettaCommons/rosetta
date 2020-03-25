// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/mhc_epitope_energy/MHCEpitopePredictorSVM.hh
/// @brief Wraps NmerSVM into the mhc_epitope_energy scoreterm framework
/// @details Just a simple layer over Indigo King's NmerSVM, to meet the MHCEpitopePredictor interface, so can be used within mhc_epitope scoreterm machinery.
/// Note that in the expected usage here, Nmer takes a 15mer, consisting of a core 9mer plus fixed-length 3mer overhangs on both N- and C- termini. The MHCEpitopeEnergy constructs and pads the peptide appropriately.
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

#ifndef INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopePredictorSVM_hh
#define INCLUDED_core_scoring_mhc_epitope_energy_MHCEpitopePredictorSVM_hh

#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictor.hh>
#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorSVM.fwd.hh>
#include <core/scoring/nmer/NMerSVMEnergy.fwd.hh>

namespace core {
namespace scoring {
namespace mhc_epitope_energy {

class MHCEpitopePredictorSVM: public MHCEpitopePredictor {

public:
	/// @brief MHCEpitopePredictorSVM constructor, taking a NMerSVMEnergyOP configured with scoring options.
	MHCEpitopePredictorSVM(core::scoring::methods::NMerSVMEnergyOP);

	bool operator==(MHCEpitopePredictor const & /* other */) override;

	std::string report() const override;

	/// @brief Scores a peptide
	core::Real score(std::string const &pep) override;

	/// @brief Accessor for the svm_ (NMerSVMEnergy object) being used by the Predictor for scoring.
	core::scoring::methods::NMerSVMEnergyOP get_svm() const {return svm_;}

private:
	core::scoring::methods::NMerSVMEnergyOP svm_; // delegate scoring

}; // class MHCEpitopePredictorSVM

} // mhc_epitope_energy
} // scoring
} // core

#endif
