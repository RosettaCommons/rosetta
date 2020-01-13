// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/mhc_epitope_energy/MHCEpitopePredictorSVM.cc
/// @brief Wraps NmerSVM into the mhc_epitope_energy scoreterm framework
/// @details Just a simple layer over Indigo King's NmerSVM, to meet the MHCEpitopePredictor interface, so can be used within mhc_epitope scoreterm machinery.
/// Note that in the expected usage here, Nmer takes a 15mer, consisting of a core 9mer plus fixed-length 3mer overhangs on both N- and C- termini. The MHCEpitopeEnergy constructs and pads the peptide appropriately.
/// @author Chris Bailey-Kellogg, cbk@cs.dartmouth.edu; Brahm Yachnin, brahm.yachnin@rutgers.edu

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

#include <core/scoring/mhc_epitope_energy/MHCEpitopePredictorSVM.hh>
#include <core/scoring/nmer/NMerSVMEnergy.hh>

namespace core {
namespace scoring {
namespace mhc_epitope_energy {

static basic::Tracer TR("core.scoring.mhc_epitope_energy");

using namespace core::scoring::methods;
using utility::vector1;

MHCEpitopePredictorSVM::MHCEpitopePredictorSVM(NMerSVMEnergyOP svm)
: svm_(svm)
{
	set_peptide_length(svm->nmer_length() + 2*svm->term_length());
	set_overhang_length(svm->term_length());
}

bool MHCEpitopePredictorSVM::operator==(MHCEpitopePredictor const &other)
{
	//Check if the Predictor type is correct
	MHCEpitopePredictorSVM const *o = dynamic_cast<MHCEpitopePredictorSVM const *>(&other);
	if ( !o ) return false;

	//svm_ is the only prviate member variable.  All configuration is stored in the NMer object.
	if ( !( *(o->svm_) == *svm_ ) ) return false;
	//Also get the base class peptide_length_ and overhang_length_
	if ( o->get_peptide_length() != get_peptide_length() ) return false;
	if ( o->get_overhang_length() != get_overhang_length() ) return false;

	return true;
}

std::string MHCEpitopePredictorSVM::report() const
{
	std::stringstream output("");

	output << "SVM " << get_peptide_length() << "mer (" << get_svm()->nmer_length() << "mer core, " << get_svm()->term_length() << "mer termini)" << std::endl;
	output << "Configured with " << get_svm()->n_svms() << " SVMs and with avg_rank_as_energy set to " << get_svm()->avg_rank_as_energy();

	return output.str();
}

core::Real MHCEpitopePredictorSVM::score(std::string const &pep)
{
	core::Size p1( 1 + svm_->term_length() ); // pep starts with the flanking residues, so p1 is dowstream by however many of those there are

	// TODO: just compute one or the other of the scores? not the whole vector but just the total?
	// Indigo says: want to have everything available at the end, so would need to special case that scenario if limit it here.
	Real rsd_energy_avg( 0. );
	Real rsd_rank_avg( 0. );
	vector1< Real > rsd_svm_energies( svm_->n_svms(), Real( 0. ) );
	vector1< Real > rsd_svm_ranks( svm_->n_svms(), Real( 0. ) );

	svm_->get_residue_energy_from_sequence( pep, p1, rsd_energy_avg, rsd_rank_avg, rsd_svm_energies, rsd_svm_ranks );

	if ( svm_->avg_rank_as_energy() ) return rsd_rank_avg;
	return rsd_energy_avg;
}

}//ns mhc_epitope_energy
}//ns scoring
}//ns core
