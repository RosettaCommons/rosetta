// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/protein_mpnn/util.hh
/// @brief Utility functions for ProteinMPNN
/// @author Brian Koepnick (koepnick@uw.edu)

#ifdef USE_TORCH

#ifndef INCLUDED_protocols_protein_mpnn_util_hh
#define INCLUDED_protocols_protein_mpnn_util_hh

// Core headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility>

#include <torch/script.h>

namespace protocols {
namespace protein_mpnn {

static std::string AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWYX";

typedef torch::List< torch::List< int64_t > > TiedPositions;

torch::Tensor
bb_coords_to_tensor( const core::conformation::Conformation & conf, core::Size batch_size );

torch::Tensor
randn_to_tensor( core::Size nres, core::Size batch_size );

torch::Tensor
seq_to_tensor( const std::string & seq, core::Size batch_size );

torch::Tensor
chain_encoding_to_tensor( const core::pose::Pose & pose, core::Size batch_size );

torch::Tensor
residue_idx_to_tensor( const core::pose::Pose & pose, core::Size batch_size );

torch::Tensor
coord_mask_to_tensor( const utility::vector1_bool & mask, core::Size batch_size );

torch::Tensor
chain_mask_to_tensor( const core::pose::Pose & pose, const utility::vector1_bool & chain_mask,  core::Size batch_size );

torch::Tensor
pos_mask_to_tensor( const utility::vector1_bool & pos_mask, core::Size batch_size );

torch::Tensor
omit_AAs_to_tensor( utility::vector1< char > omit_AAs );

torch::Tensor
omit_AAs_pos_to_tensor( const utility::vector1< utility::vector1< char > > & omit_AAs_pos, const core::Size batch_size );

torch::Tensor
bias_AAs_to_tensor( utility::vector1< core::Real > bias_AAs );

TiedPositions
convert_tied_positions( utility::vector1< utility::vector1< core::Size > > tied_positions );

utility::vector1< std::string >
seqs_from_tensor( const torch::Tensor & seq_t, core::Size num_res, core::Size batch_size=1 );

} //protein_mpnn
} //protocols

#endif //INCLUDED_protocols_protein_mpnn_util_hh

#endif //USE_TORCH
