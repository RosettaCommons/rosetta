// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ptm_prediction/PTMPredictionTensorflowProtocolBase.cc
/// @brief A pure virtual class, derived from RosettaTensorflowProtocolBase, which will serve as a base for protocols
/// for predicting post-translational modifications in proteins.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// The code for a Tensorflow based metric has been adpated from:
// protocols/cyclic_peptide/peptide_fold_propensity_metric/Peptide10merFoldPropensityTensorflowProtocol_v1
// which was authored by Vikram K. Mulligan (PR #3936)

// Project headers:
#include <protocols/ptm_prediction/PTMPredictionTensorflowProtocolBase.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.ptm_prediction.PTMPredictionTensorflowProtocolBase" );


namespace protocols {
namespace ptm_prediction {

/// @brief Destructor.
PTMPredictionTensorflowProtocolBase::~PTMPredictionTensorflowProtocolBase(){}

#ifdef USE_TENSORFLOW
/// @brief Allow derived classes to set the Tensorflow session.
/// @details Should only be called from derived class constructors!
void
PTMPredictionTensorflowProtocolBase::set_tensorflow_session(
	basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP session_in
) {
	debug_assert( session_in != nullptr );
	tensorflow_session_ = session_in;
}
#endif //USE_TENSORFLOW

/// @brief Given a selection, return the number of selected residues.
/// @details Static, protected.
utility::vector1< core::Size >
PTMPredictionTensorflowProtocolBase::get_selected_residues(
	core::select::residue_selector::ResidueSubset const & selected
) {
	utility::vector1< core::Size > returnvec;
	for ( core::Size i(1), imax(selected.size()); i<=imax; ++i ) {
		if ( selected[i] ) returnvec.push_back(i);
	}
	return returnvec;
}

} //ptm_prediction
} //protocols
