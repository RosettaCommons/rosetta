// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/trRosetta/trRosettaOutputs_v1.cc
/// @brief A class for the outputs of trRosetta version 1.  This stores distance, phi, theta, and omega probability distributions.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

// Project headers:
#include <protocols/trRosetta/trRosettaOutputs_v1.hh>

// Basic headers:
#include <basic/Tracer.hh>
#ifndef USE_TENSORFLOW
#include <basic/tensorflow_manager/util.hh>
#endif

// Utility headers:
#include <utility/pointer/memory.hh>

#ifdef USE_TENSORFLOW
#include <utility/fixedsizearray1.hh>
#endif

static basic::Tracer TR( "protocols.trRosetta.trRosettaOutputs_v1" );

namespace protocols {
namespace trRosetta {

#ifndef USE_TENSORFLOW
/// @details Default constructor produces error in non-tensorflow build.
trRosettaOutputs_v1::trRosettaOutputs_v1() {
	utility_exit_with_message( basic::tensorflow_manager::get_tensorflow_compilation_instructions( "trRosettaOutputs_v1 class" ) );
}
#endif

#ifdef USE_TENSORFLOW
/// @brief Options constructor, used to package data into a convenient object.
trRosettaOutputs_v1::trRosettaOutputs_v1(
	core::Size const seqlength_in,
	utility::vector1< basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > > const & output_vec
) :
	trRosettaOutputsBase(),
	seqlength_(seqlength_in),
	dist_tensor_(utility::vector1<int64_t>{static_cast<int64_t>(seqlength_in), static_cast<int64_t>(seqlength_in), trRosetta_v1_expected_dist_bins}),
	phi_tensor_(utility::vector1<int64_t>{static_cast<int64_t>(seqlength_in), static_cast<int64_t>(seqlength_in), trRosetta_v1_expected_phi_bins}),
	theta_tensor_(utility::vector1<int64_t>{static_cast<int64_t>(seqlength_in), static_cast<int64_t>(seqlength_in), trRosetta_v1_expected_theta_bins}),
	omega_tensor_(utility::vector1<int64_t>{static_cast<int64_t>(seqlength_in), static_cast<int64_t>(seqlength_in), trRosetta_v1_expected_omega_bins})
{
	initialize_internal_tensors( output_vec );
}
#endif //USE_TENSORFLOW

/// @brief Destructor.
trRosettaOutputs_v1::~trRosettaOutputs_v1() = default;

/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
trRosettaOutputs_v1OP
trRosettaOutputs_v1::clone() const {
	return utility::pointer::make_shared< trRosettaOutputs_v1 >( *this );
}

////////////////////////////////////////
// PUBLIC MEMBER FUNCTIONS:
////////////////////////////////////////

#ifdef USE_TENSORFLOW

/// @brief Get the number of distance bins (37).
core::Size
trRosettaOutputs_v1::n_dist_bins() const {
	return trRosetta_v1_expected_dist_bins;
}

/// @brief Get the number of phi bins (13).
core::Size
trRosettaOutputs_v1::n_phi_bins() const {
	return trRosetta_v1_expected_phi_bins;
}

/// @brief Get the number of theta bins (25).
core::Size
trRosettaOutputs_v1::n_theta_bins() const {
	return trRosetta_v1_expected_theta_bins;
}

/// @brief Get the number of omega bins (25).
core::Size
trRosettaOutputs_v1::n_omega_bins() const {
	return trRosetta_v1_expected_omega_bins;
}

/// @brief Access the kth distance bin for residues i and j.
/// @details Access is 1-based.  The bin index, k, must be from 1 to 37.
float
trRosettaOutputs_v1::dist(
	core::Size const res_i,
	core::Size const res_j,
	core::Size const bin_k
) const {
	debug_assert( res_i > 0 && res_i <= seqlength_ );
	debug_assert( res_j > 0 && res_j <= seqlength_ );
	debug_assert( bin_k > 0 && bin_k <= trRosetta_v1_expected_dist_bins );
	return dist_tensor_( res_i, res_j, bin_k );
}

/// @brief Access the kth phi bin for residues i and j.
/// @details Access is 1-based.  The bin index, k, must be from 1 to 13.
float
trRosettaOutputs_v1::phi(
	core::Size const res_i,
	core::Size const res_j,
	core::Size const bin_k
) const {
	debug_assert( res_i > 0 && res_i <= seqlength_ );
	debug_assert( res_j > 0 && res_j <= seqlength_ );
	debug_assert( bin_k > 0 && bin_k <= trRosetta_v1_expected_phi_bins );
	return phi_tensor_( res_i, res_j, bin_k );
}

/// @brief Access the kth theta bin for residues i and j.
/// @details Access is 1-based.  The bin index, k, must be from 1 to 25.
float
trRosettaOutputs_v1::theta(
	core::Size const res_i,
	core::Size const res_j,
	core::Size const bin_k
) const {
	debug_assert( res_i > 0 && res_i <= seqlength_ );
	debug_assert( res_j > 0 && res_j <= seqlength_ );
	debug_assert( bin_k > 0 && bin_k <= trRosetta_v1_expected_theta_bins );
	return theta_tensor_( res_i, res_j, bin_k );
}

/// @brief Access the kth omega bin for residues i and j.
/// @details Access is 1-based.  The bin index, k, must be from 1 to 25.
float
trRosettaOutputs_v1::omega(
	core::Size const res_i,
	core::Size const res_j,
	core::Size const bin_k
) const {
	debug_assert( res_i > 0 && res_i <= seqlength_ );
	debug_assert( res_j > 0 && res_j <= seqlength_ );
	debug_assert( bin_k > 0 && bin_k <= trRosetta_v1_expected_omega_bins );
	return omega_tensor_( res_i, res_j, bin_k );
}

#endif //USE_TENSORFLOW

////////////////////////////////////////
// PRIVATE MEMBER FUNCTIONS:
////////////////////////////////////////

#ifdef USE_TENSORFLOW

/// @brief Called by constructor to copy data from output_vec to internal tensors.
void
trRosettaOutputs_v1::initialize_internal_tensors(
	utility::vector1< basic::tensorflow_manager::RosettaTensorflowTensorContainer< float > > const & output_vec
) {
	using namespace basic::tensorflow_manager;

	//Expected bin counts:
	utility::fixedsizearray1< core::Size, 4 > const expected_bins{ trRosetta_v1_expected_dist_bins, trRosetta_v1_expected_omega_bins, trRosetta_v1_expected_phi_bins, trRosetta_v1_expected_theta_bins };

	//The tensors:
	utility::fixedsizearray1< RosettaTensorflowTensorContainer< float > *, 4 > const tensors{ &dist_tensor_, &omega_tensor_, &phi_tensor_, &theta_tensor_ };

	//Checks first:
#ifndef NDEBUG
	{
		debug_assert( output_vec.size() == 4 );
		for( core::Size i(1); i<=4; ++i ) {
			debug_assert( output_vec[i].n_dimensions() == 3 );
			debug_assert( output_vec[i].dimension(1) == seqlength_ );
			debug_assert( output_vec[i].dimension(2) == seqlength_ );
			debug_assert( output_vec[i].dimension(3) == expected_bins[i] );
		}
	}
#endif

	// Copy each tensor:
	for( core::Size i(1); i<=4; ++i ) {
		RosettaTensorflowTensorContainer< float > * curtensor( tensors[i] );
		for( core::Size j(1); j<=seqlength_; ++j ) {
			for( core::Size k(1); k<=seqlength_; ++k) {
				for( core::Size l(1); l<=expected_bins[i]; ++l ) {
					(*curtensor)(j,k,l) = output_vec[i](j,k,l);
				}
			}
		}
	}
}

#endif //USE_TENSORFLOW

} //trRosetta
} //protocols
