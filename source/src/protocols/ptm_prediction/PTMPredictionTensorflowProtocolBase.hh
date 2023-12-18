// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ptm_prediction/PTMPredictionTensorflowProtocolBase.hh
/// @brief A pure virtual class, derived from RosettaTensorflowProtocolBase, which will serve as a base for protocols
/// for predicting post-translational modifications in proteins.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// The code for a Tensorflow based metric has been adpated from:
// protocols/cyclic_peptide/peptide_fold_propensity_metric/Peptide10merFoldPropensityTensorflowProtocol_v1
// which was authored by Vikram K. Mulligan (PR #3936)


#ifndef INCLUDED_protocols_ptm_prediction_PTMPredictionTensorflowProtocolBase_hh
#define INCLUDED_protocols_ptm_prediction_PTMPredictionTensorflowProtocolBase_hh

#include <protocols/ptm_prediction/PTMPredictionTensorflowProtocolBase.fwd.hh>
#include <basic/tensorflow_manager/RosettaTensorflowProtocolBase.hh>

// Basic headers
#ifdef USE_TENSORFLOW
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.fwd.hh>
#endif //USE_TENSORFLOW

// Core headers
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
//#include <utility/pointer/ReferenceCount.hh>
// C++ header
#include <map>

// protocol headers
#include <protocols/ptm_prediction/PTMPredictionMetric.hh>

namespace protocols {
namespace ptm_prediction {

/// @brief A pure virtual class, derived from RosettaTensorflowProtocolBase, which will serve as a base for protocols
/// for predicting post-translational modifications in proteins.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)
class PTMPredictionTensorflowProtocolBase : public basic::tensorflow_manager::RosettaTensorflowProtocolBase {

public:

	/// @brief Default constructor.
	PTMPredictionTensorflowProtocolBase() = default;

	/// @brief Copy constructor.
	PTMPredictionTensorflowProtocolBase( PTMPredictionTensorflowProtocolBase const & ) = default;

	/// @brief Destructor.
	~PTMPredictionTensorflowProtocolBase() override;

public:
#ifdef USE_TENSORFLOW
	/// @brief Analyze a pose and return its post-translational modifcation probabilities (ranging from 0 to 1).
	/// @details Analyze a pose & return its post-translational modification probabilities (0-1)
	virtual
	std::map< core::Size, core::Real >
	compute_ptm_probability(
		core::pose::Pose const &pose,
		core::select::residue_selector::ResidueSubset const & selection,
//		basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP & session,
		PTMPredictionMetric::PTMPredictionMetricModes const ptm_mode
	) const = 0;

	/// @brief Analyze a given Asparagine in a pose, returning an deamidation (ranging from 0 to 1).
	/// @details Uses a neural net to predict deamidation based on structural/sequence features.
	virtual
	std::map< core::Size, core::Real >
	compute_deamidation_probability(
		core::pose::Pose const &pose,
		core::select::residue_selector::ResidueSubset const & selection
//		basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP & session
	) const = 0;

#endif // USE_TENSORFLOW
	//protected:

#ifdef USE_TENSORFLOW

protected:
	/// @brief Allow derived classes to access the Tensorflow session.
	inline basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP tensorflow_session() const { return tensorflow_session_; }

	/// @brief Allow derived classes to set the Tensorflow session.
	void set_tensorflow_session( basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP session_in );
#endif //USE_TENSORFLOW

	/// @brief Given a selection, return a vector of selected indices.
	static utility::vector1< core::Size > get_selected_residues( core::select::residue_selector::ResidueSubset const & selected);


private:

#ifdef USE_TENSORFLOW
	/// @brief The tensorflow session.
	basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP tensorflow_session_;
#endif //USE_TENSORFLOW

};

} //ptm_prediction
} //protocols

#endif //INCLUDED_protocols_ptm_prediction_PTMPredictionTensorflowProtocolBase_hh
