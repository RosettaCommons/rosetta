// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/inverse_folding/util.hh
/// @brief Helper functions for using the MIF-ST model developed Yang et al.
/// @author moritzertelt (moritzertelt)

#ifndef INCLUDED_protocols_inverse_folding_util_hh
#define INCLUDED_protocols_inverse_folding_util_hh

#include <protocols/inverse_folding/MIFST.fwd.hh>

// Core headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/util.hh>
#include <core/pose/Pose.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

#ifdef USE_TORCH
// Torch header
#include <torch/script.h>
#endif //USE_TORCH

// Eigen header
#include <Eigen/Dense>

// std headers
#include <cmath>

namespace protocols {
namespace inverse_folding {

// constant variables defining the alphabet used
static std::string PROTEIN_ALPHABET = "ACDEFGHIKLMNPQRSTVWYBZXJOU-*#@!";
static char MASK = '#';

///@brief Helper function to return trRosetta style dihedrals
core::Real
get_dihedrals(Eigen::Vector3d const & a, Eigen::Vector3d const & b, Eigen::Vector3d const & c, Eigen::Vector3d const & d);

///@brief helper function to return trRosetta style angles
core::Real
get_angles(Eigen::Vector3d const & a, Eigen::Vector3d const & b, Eigen::Vector3d const & c);

///@brief get the coordinates of a specified ResidueSubset
std::vector<Eigen::Vector3d>
extract_coordinates(core::pose::Pose const & pose, core::select::residue_selector::ResidueSubset const & subset, std::string const & atom_name);

///@brief Calculate the Cb position following the implementation of Yang et al.
Eigen::Vector3d
calculate_Cb_position(Eigen::Vector3d const & N, Eigen::Vector3d const & Ca, Eigen::Vector3d const & C);

///@brief Calculate the distance matrix of C betas
Eigen::MatrixXd
compute_distance_matrix(std::vector<Eigen::Vector3d> const & coordinates);

///@brief function to tokenize sequence into MIF-ST Alphabet
std::vector<core::Size>
tokenize( std::string const & seq, core::Size mask_position);

#ifdef USE_TORCH

///@brief Process coordinates into features
std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor>
process_coords(core::pose::Pose const & pose, core::select::residue_selector::ResidueSubset const & subset);

///@brief replace NaNs in tensor with zero
torch::Tensor replace_nan(torch::Tensor E);

///@brief prepare the sequences for model input
torch::Tensor
prepare_sequences(std::vector<std::string> const & sequences, std::vector<core::Size> const & mask_positions);

///@brief get node features from trRosetta style dihedrals, has bool flag to define whether sines/cosines are provided
torch::Tensor
get_node_features(
        const torch::Tensor& omega,
        const torch::Tensor& theta,
        const torch::Tensor& phi,
        bool sc = false
);

///@brief get the k neighbors
torch::Tensor
get_k_neighbors(const torch::Tensor& dist, int64_t k);

///@brief get the edge features from trRosetta style dihedrals
torch::Tensor
get_edge_features(
        const torch::Tensor& dist,
        const torch::Tensor& omega,
        const torch::Tensor& theta,
        const torch::Tensor& phi,
        const torch::Tensor& E_idx,
        bool sc = false
);

///@brief create a mask for missing values
torch::Tensor
get_mask(torch::Tensor const & E);

///@brief prepare the features for the model call
std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor>
collate_structure(
        std::vector<std::string> const & sequences,
        torch::Tensor const & dist,
        torch::Tensor const & omega,
        torch::Tensor const & theta,
        torch::Tensor const & phi,
        std::vector<core::Size> const & mask_positions,
        core::Size n_connections = 20,
        core::Size n_node_features = 10,
        core::Size n_edge_features = 11
);

#endif //USE_TORCH
} //inverse_folding
} //protocols


#endif //protocols/inverse_folding_util_hh

