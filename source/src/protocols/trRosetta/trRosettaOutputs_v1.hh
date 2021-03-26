// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta/trRosettaOutputs_v1.hh
/// @brief A class for the outputs of trRosetta version 1.  This stores distance, phi, theta, and omega probability distributions.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_trRosetta_trRosettaOutputs_v1_hh
#define INCLUDED_protocols_trRosetta_trRosettaOutputs_v1_hh

#include <protocols/trRosetta/trRosettaOutputs_v1.fwd.hh>
#include <protocols/trRosetta/trRosettaOutputsBase.hh>

// Core headers:
#include <core/types.hh>

// Basic headers
#ifdef USE_TENSORFLOW
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.tmpl.hh>
#endif //USE_TENSORFLOW

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace trRosetta {

// Tensorflow model-specific defines:
#define trRosetta_v1_expected_dist_bins 37
#define trRosetta_v1_expected_phi_bins 13
#define trRosetta_v1_expected_theta_bins 25
#define trRosetta_v1_expected_omega_bins 25

/// @brief A class for the outputs of trRosetta version 1.  This stores distance, phi, theta, and omega probability distributions.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class trRosettaOutputs_v1 : public protocols::trRosetta::trRosettaOutputsBase {

public:

	/// @brief Default constructor, explicitly deleted in Tensorlfow build.
#ifdef USE_TENSORFLOW
	trRosettaOutputs_v1() = delete;
#else //Not USE_TENSORFLOW
	trRosettaOutputs_v1();
#endif

#ifdef USE_TENSORFLOW
	/// @brief Options constructor, used to package data into a convenient object.
	trRosettaOutputs_v1( core::Size const seqlength_in, utility::vector1< basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > > const & output_vec );
#endif

	/// @brief Destructor.
	~trRosettaOutputs_v1() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	trRosettaOutputs_v1OP clone() const;

public: //Accessors:

#ifdef USE_TENSORFLOW
	/// @brief Get the sequence length.
	inline core::Size sequence_length() const { return seqlength_; }

	/// @brief Get the number of distance bins (37).
	core::Size n_dist_bins() const;

	/// @brief Get the number of phi bins (13).
	core::Size n_phi_bins() const;

	/// @brief Get the number of theta bins (25).
	core::Size n_theta_bins() const;

	/// @brief Get the number of omega bins (25).
	core::Size n_omega_bins() const;

	/// @brief Access the kth distance bin for residues i and j.
	/// @details Access is 1-based.  The bin index, k, must be from 1 to 37.
	float dist( core::Size const res_i, core::Size const res_j, core::Size const bin_k ) const;

	/// @brief Access the kth phi bin for residues i and j.
	/// @details Access is 1-based.  The bin index, k, must be from 1 to 13.
	float phi( core::Size const res_i, core::Size const res_j, core::Size const bin_k ) const;

	/// @brief Access the kth theta bin for residues i and j.
	/// @details Access is 1-based.  The bin index, k, must be from 1 to 25.
	float theta( core::Size const res_i, core::Size const res_j, core::Size const bin_k ) const;

	/// @brief Access the kth omega bin for residues i and j.
	/// @details Access is 1-based.  The bin index, k, must be from 1 to 25.
	float omega( core::Size const res_i, core::Size const res_j, core::Size const bin_k ) const;

#endif //USE_TENSORFLOW

private: //Functions:

#ifdef USE_TENSORFLOW
	/// @brief Called by constructor to copy data from output_vec to internal tensors.
	void initialize_internal_tensors( utility::vector1< basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > > const & output_vec );
#endif //USE_TENSORFLOW

private: //Data:

#ifdef USE_TENSORFLOW
	/// @brief The sequence length:
	core::Size seqlength_ = 0;

	/// @brief The distance tensor:
	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > dist_tensor_;

	/// @brief The phi tensor:
	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > phi_tensor_;

	/// @brief The theta tensor:
	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > theta_tensor_;

	/// @brief The omega tensor:
	basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > omega_tensor_;
#endif //USE_TENSORFLOW

};

} //trRosetta
} //protocols

#endif //INCLUDED_protocols_trRosetta_trRosettaOutputs_v1_hh
