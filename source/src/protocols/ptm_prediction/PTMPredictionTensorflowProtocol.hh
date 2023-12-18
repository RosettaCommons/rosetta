// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ptm_prediction/PTMPredictionTensorflowProtocol.hh
/// @brief A class for predicting post-translational modifications using a Neural Net.
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)

// The code for a Tensorflow based metric has been adpated from:
// protocols/cyclic_peptide/peptide_fold_propensity_metric/Peptide10merFoldPropensityTensorflowProtocol_v1
// which was authored by Vikram K. Mulligan (PR #3936)
#ifndef INCLUDED_protocols_ptm_prediction_PTMPredictionTensorflowProtocol_hh
#define INCLUDED_protocols_ptm_prediction_PTMPredictionTensorflowProtocol_hh


#ifdef USE_TENSORFLOW
#include <tensorflow/c/c_api.h>
#endif

#include <protocols/ptm_prediction/PTMPredictionTensorflowProtocolBase.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
//#include <utility/pointer/ReferenceCount.hh>

// Basic headers
#include<basic/tensorflow_manager/RosettaTensorflowTensorContainer.hh>
// C++ header
#include <map>

// protocol headers
#include <protocols/ptm_prediction/PTMPredictionMetric.hh>

namespace protocols {
namespace ptm_prediction {

/// @brief Analyze a given residue or residues in a pose, returning an modifcation probability (ranging from 0 to 1).
/// @author Moritz Ertelt (moritz.ertelt@googlemail.com)
class PTMPredictionTensorflowProtocol : public PTMPredictionTensorflowProtocolBase {


public:

	/// @brief Default constructor.
	PTMPredictionTensorflowProtocol();

	/// @brief Copy constructor.
	PTMPredictionTensorflowProtocol(PTMPredictionTensorflowProtocol const & src);

	/// @brief Destructor.
	~PTMPredictionTensorflowProtocol() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	basic::tensorflow_manager::RosettaTensorflowProtocolBaseOP clone() const override;

	/// @brief Get the name of this protocol.
	std::string name() const override;

	/// @brief path to the model
	std::string path_to_model;

	/// @brief constructor with path to model.
	PTMPredictionTensorflowProtocol( std::string path );

public:

#ifdef USE_TENSORFLOW
	/// @brief Analyze a given residue or residues in a pose, returning an modifcation probability (ranging from 0 to 1).
	std::map< core::Size, core::Real > compute_deamidation_probability( core::pose::Pose const &pose,
		core::select::residue_selector::ResidueSubset const & selection ) const override;

	/// @brief Analyze a given residue or residues in a pose, returning an modifcation probability (ranging from 0 to 1).
	std::map< core::Size, core::Real > compute_ptm_probability( core::pose::Pose const &pose,
		core::select::residue_selector::ResidueSubset const & selection,
		PTMPredictionMetric::PTMPredictionMetricModes const ptm_mode ) const override;

	/// @brief Get the Tensorflow session used by this class.
	/// @details Implemented by derived classes.
	basic::tensorflow_manager::RosettaTensorflowSessionContainerCOP get_tensorflow_session( ) const;
#endif //USE_TENSORFLOW

private: //Methods

#ifdef USE_TENSORFLOW
/// @brief Given a pose and an already-allocated (but empty) input tensor, store the relevant pose data for the
/// deamidation site prediction in the tensor.
/// @details The input tensor is a 1D matrix. Following features are stored:
/// 1: SASA of the selected Asn.
/// 2: Phi angle of the selected Asn.
/// 3: Psi angle of the selected Asn.
/// 4: Chi_1 angle of the selected Asn.
/// 5: Chi_2 angle of the selected Asn.
/// 6: Attack distance between C-Beta of Asn and next residue N atom
/// 7: Deamidation rates from Robinson, N. E. et al. Journal of Peptide Research (2001)
/// 8-26: one-hot encoded amino acid
	void
	copy_deamidation_features_from_pose_to_tensor(
		core::pose::Pose const & pose,
		core::Size const position,
		basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > & input_tensor,
		std::vector< std::vector< float > > & deamidation_rates
	) const;

/// @brief Given a pose and an already-allocated (but empty) input tensor, store the relevant pose data for the
/// ptm prediction in the tensor.
/// @details The input tensor is a 1D matrix. Following features are stored:
/// Sequence window of -4/+4 residues, phi/psi angles window -1/+1, SASA, E/H/L secondary structure
  void
  copy_feat_to_tensor(
  	core::pose::Pose const & pose,
  	core::Size const & position,
  	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > & input_tensor_seq,
		basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > & input_tensor_struc
  ) const;

	/// @brief Read in the CSV table describing deamidation rates based on Robinson, N.E et al. (2001) J. Pept. Res.
	std::vector< std::vector< float > > parseCSV() const;
	/// @brief Check the type of modification and corresponding type of amino acid and output tensor
	void
	setup_and_check(
		std::string const & aa_select,
		PTMPredictionMetric::PTMPredictionMetricModes const & ptm_mode,
		basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > & output_tensor,
		size_t & index,
		utility::vector1< std::string > & input_names
	) const;

#endif // USE_TENSORFLOW
private: //Data

};

} //ptm_prediction
} //protocols

#endif //INCLUDED_protocols_ptm_prediction_PTMPredictionTensorflowProtocol_hh
