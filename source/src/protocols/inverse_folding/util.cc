// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/inverse_folding/util.cc
/// @brief Helper functions for using the MIF-ST model developed Yang et al.
/// @author moritzertelt (moritzertelt)

#include <protocols/inverse_folding/util.hh>

// core headers
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>

// numeric headers
#include <numeric/xyzMatrix.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.inverse_folding.util" );



namespace protocols {
namespace inverse_folding {

core::Real
get_dihedrals(Eigen::Vector3d const & a, Eigen::Vector3d const & b, Eigen::Vector3d const & c, Eigen::Vector3d const & d) {
	Eigen::Vector3d b0 = -1.0 * (b - a);
	Eigen::Vector3d b1 = c - b;
	Eigen::Vector3d b2 = d - c;

	b1 /= b1.norm();

	Eigen::Vector3d v = b0 - b0.dot(b1) * b1;
	Eigen::Vector3d w = b2 - b2.dot(b1) * b1;

	core::Real x = v.dot(w);
	core::Real y = b1.cross(v).dot(w);

	return std::atan2(y, x);
}

core::Real
get_angles(Eigen::Vector3d const & a, Eigen::Vector3d const & b, Eigen::Vector3d const & c) {
	Eigen::Vector3d v = a - b;
	v /= v.norm();

	Eigen::Vector3d w = c - b;
	w /= w.norm();

	core::Real x = v.dot(w);

	// Make sure x is in the range [-1, 1] to avoid NaN from std::acos
	x = std::max(-1.0, std::min(1.0, x));

	return std::acos(x);
}


std::vector<Eigen::Vector3d>
extract_coordinates(core::pose::Pose const & pose, core::select::residue_selector::ResidueSubset const & subset, std::string const & atom_name) {
	std::vector<Eigen::Vector3d> coordinates;
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( subset[i] ) {
			core::conformation::Residue const & residue = pose.residue(i);
			if ( residue.has(atom_name) ) {
				core::Vector const & rosetta_coord = residue.xyz(atom_name);
				coordinates.emplace_back(rosetta_coord.x(), rosetta_coord.y(), rosetta_coord.z());
			} else {
				TR.Warning << "Missing coordinates at position " << i << std::endl;
			}
		}
	}
	return coordinates;
}

Eigen::Vector3d
calculate_Cb_position(Eigen::Vector3d const & N, Eigen::Vector3d const & Ca, Eigen::Vector3d const & C) {
	Eigen::Vector3d b = Ca - N;
	Eigen::Vector3d c = C - Ca;
	Eigen::Vector3d a = b.cross(c);
	// numbers taken from the original implementation by Yang et al.
	return Ca - 0.58273431 * a + 0.56802827 * b - 0.54067466 * c;
}

Eigen::MatrixXd
compute_distance_matrix(std::vector<Eigen::Vector3d> const & coordinates) {
	auto n = static_cast< Eigen::Index >( coordinates.size() );
	Eigen::MatrixXd dists(n, n);
	for ( Eigen::Index i = 0; i < n; ++i ) {
		for ( Eigen::Index j = 0; j < n; ++j ) {
			dists(i, j) = (coordinates[i] - coordinates[j]).norm();
		}
	}
	return dists;
}




std::vector<core::Size>
tokenize( std::string const & seq, core::Size mask_position) {
	std::vector<core::Size> tokenized_seq;

	for ( core::Size i = 0; i < static_cast<core::Size>(seq.size()); ++i ) {
		char c = ( i == mask_position ) ? MASK : seq[i];  // '#' is the mask character
		core::Size index = PROTEIN_ALPHABET.find(c);

		if ( index != std::string::npos ) {
			tokenized_seq.push_back(index);
		} else {
			TR.Warning << " Residue was not recognized at tokenization step, masking it instead!" << std::endl;
			tokenized_seq.push_back( 28 );
		}
	}
	return tokenized_seq;
}

#ifdef USE_TORCH
std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor>
process_coords(core::pose::Pose const & pose, core::select::residue_selector::ResidueSubset const & subset) {

    // Extract coordinates
    std::vector<Eigen::Vector3d> N_coords = extract_coordinates(pose, subset, "N");
    std::vector<Eigen::Vector3d> Ca_coords = extract_coordinates(pose, subset, "CA");
    std::vector<Eigen::Vector3d> C_coords = extract_coordinates(pose, subset, "C");

    auto nres = static_cast< Eigen::Index >( N_coords.size() );

    // Recreate Cb given N, Ca, C
    std::vector<Eigen::Vector3d> Cb_coords;
    Cb_coords.reserve(nres);
    for (Eigen::Index i = 0; i < nres; ++i) {
        Cb_coords.push_back(calculate_Cb_position(N_coords[i], Ca_coords[i], C_coords[i]));
    }

    // Cb-Cb distance matrix
    Eigen::MatrixXd dist = compute_distance_matrix(Cb_coords);
    for (Eigen::Index i = 0; i < nres; ++i) {
        dist(i, i) = std::nan(""); // Set diagonal to NaN
    }

    // Matrices for dihedrals and angles
    // Initialize PyTorch tensors
    auto options = torch::dtype(torch::kFloat32);
    torch::Tensor dist_tensor = torch::full({nres, nres}, std::nan(""), options);
    torch::Tensor omega_tensor = torch::full({nres, nres}, std::nan(""), options);
    torch::Tensor theta_tensor = torch::full({nres, nres}, std::nan(""), options);
    torch::Tensor phi_tensor = torch::full({nres, nres}, std::nan(""), options);

    for (Eigen::Index i = 0; i < nres; ++i) {
        for (Eigen::Index j = 0; j < nres; ++j) {
            if (i != j) {
                dist_tensor[i][j] = dist(i, j);
                omega_tensor[i][j] = get_dihedrals(Ca_coords[i], Cb_coords[i], Cb_coords[j], Ca_coords[j]);
                theta_tensor[i][j] = get_dihedrals(N_coords[i], Ca_coords[i], Cb_coords[i], Cb_coords[j]);
                phi_tensor[i][j] = get_angles(Ca_coords[i], Cb_coords[i], Cb_coords[j]);
            }
        }
    }

    return std::make_tuple(dist_tensor, omega_tensor, theta_tensor, phi_tensor);
}

torch::Tensor
replace_nan(torch::Tensor E) {
    // Create a tensor that holds boolean values indicating whether each element is NaN
    torch::Tensor isnan = torch::isnan(E);
    // Replace NaN values with zeros
    E.index_put_({isnan}, 0.0);
    return E;
}

torch::Tensor
prepare_sequences(std::vector<std::string> const & sequences, std::vector<core::Size> const & mask_positions) {

    std::vector<torch::Tensor> tokenized_sequences;

    for (core::Size i = 0; i < sequences.size(); ++i) {
        auto const & seq = sequences[i];
        auto const mask_position = mask_positions[i];

        std::vector<core::Size> tokens = tokenize(seq, mask_position);

        // Convert to std::vector<int64_t>
        std::vector<int64_t> tokens_int64(tokens.begin(), tokens.end());

        torch::Tensor token_tensor = torch::tensor(tokens_int64);
        tokenized_sequences.push_back(token_tensor);
    }

    // Since all sequences have the same length, we can directly stack them
    return torch::stack(tokenized_sequences);
}

torch::Tensor
get_node_features(
        torch::Tensor const & omega,
        torch::Tensor const & theta,
        torch::Tensor const & phi,
        bool sc
) {
    // Check if all tensors are zeros
    if (omega.sum().item<float>() == 0.0 && theta.sum().item<float>() == 0.0 && phi.sum().item<float>() == 0.0) {
        return torch::zeros({omega.size(0), 10});
    }

    auto get_features = [&](
            const torch::Tensor& omega,
            const torch::Tensor& theta,
            const torch::Tensor& phi
    ) -> torch::Tensor {
        auto options = omega.options();
        auto n1 = torch::cat({torch::tensor({0.0}, options), omega.diagonal(1)});
        auto n2 = torch::cat({theta.diagonal(1), torch::tensor({0.0}, options)});
        auto n3 = torch::cat({torch::tensor({0.0}, options), theta.diagonal(-1)});
        auto n4 = torch::cat({phi.diagonal(1), torch::tensor({0.0}, options)});
        auto n5 = torch::cat({torch::tensor({0.0}, options), phi.diagonal(-1)});
        return torch::stack({n1, n2, n3, n4, n5}, 1);
    };

    torch::Tensor s, c;
    if (!sc) {
        auto ns = get_features(omega, theta, phi);
        s = torch::sin(ns);
        c = torch::cos(ns);
    } else {
        s = get_features(omega.index({0}), theta.index({0}), phi.index({0}));
        c = get_features(omega.index({1}), theta.index({1}), phi.index({1}));
    }

    return torch::cat({s, c}, 1);
}

torch::Tensor
get_k_neighbors(const torch::Tensor& dist, int64_t k) {
    // Ensure k is within bounds
    k = std::min(k, dist.size(0) - 1);

    // Use topk to get the smallest k values and their indices
    auto topk_results = torch::topk(dist, k, -1, false);
    auto idx = std::get<1>(topk_results);

    return idx;
}

torch::Tensor
get_edge_features(
        const torch::Tensor& dist,
        const torch::Tensor& omega,
        const torch::Tensor& theta,
        const torch::Tensor& phi,
        const torch::Tensor& E_idx,
        bool sc
) {
    // Check if all tensors are zeros
    if (omega.sum().item<float>() == 0.0 && theta.sum().item<float>() == 0.0 && phi.sum().item<float>() == 0.0) {
        return torch::zeros({omega.size(0), E_idx.size(1), 6});
    }

    auto get_features = [&](
            const torch::Tensor& omega,
            const torch::Tensor& theta,
            const torch::Tensor& phi,
            const torch::Tensor& E_idx
    ) -> std::vector<torch::Tensor> {
        std::vector<torch::Tensor> omega_E, theta_E, theta_Er, phi_E, phi_Er;

        for (int64_t i = 0; i < E_idx.size(0); ++i) {
            omega_E.push_back(omega.index({i, E_idx[i]}));
            theta_E.push_back(theta.index({i, E_idx[i]}));
            theta_Er.push_back(theta.index({E_idx[i], i}));
            phi_E.push_back(phi.index({i, E_idx[i]}));
            phi_Er.push_back(phi.index({E_idx[i], i}));
        }

        return {torch::stack(omega_E), torch::stack(theta_E), torch::stack(theta_Er), torch::stack(phi_E), torch::stack(phi_Er)};
    };

    std::vector<torch::Tensor> dist_E;
    for (int64_t i = 0; i < E_idx.size(0); ++i) {
        dist_E.push_back(dist.index({i, E_idx[i]}));
    }
    auto dist_E_tensor = torch::stack(dist_E);

    std::vector<torch::Tensor> angles, s, c;
    if (!sc) {
        angles = get_features(omega, theta, phi, E_idx);
        for (const auto& a : angles) {
            s.push_back(torch::sin(a));
            c.push_back(torch::cos(a));
        }
    } else {
        s = get_features(omega.index({0}), theta.index({0}), phi.index({0}), E_idx);
        c = get_features(omega.index({1}), theta.index({1}), phi.index({1}), E_idx);
    }

    s.insert(s.begin(), dist_E_tensor);
    s.insert(s.end(), c.begin(), c.end());

    return torch::stack(s, 2);
}

torch::Tensor
get_mask(torch::Tensor const & E) {
    // Sum along the last dimension and check if the sum is finite
    auto summed_E = torch::sum(E, /*dim=*/-1);
    auto mask_E = torch::isfinite(summed_E).to(torch::kFloat);

    // Reshape the mask to have an additional dimension of size 1 at the end
    auto sizes = mask_E.sizes().vec();
    sizes.push_back(1);
    mask_E = mask_E.view(torch::IntArrayRef(sizes));

    return mask_E;
}

std::tuple<torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor, torch::Tensor>
collate_structure(
        std::vector<std::string> const & sequences,
        torch::Tensor const & dist,
        torch::Tensor const & omega,
        torch::Tensor const & theta,
        torch::Tensor const & phi,
        std::vector<core::Size> const & mask_positions,
        core::Size n_connections,
        core::Size n_node_features,
        core::Size n_edge_features
) {
    // Prepare sequences
    torch::Tensor collated_seqs = prepare_sequences(sequences, mask_positions);

    core::Size max_ell = sequences[0].size();

    core::Size n = sequences.size();

    // Calculate features once, outside the loop
    torch::Tensor V = get_node_features(omega, theta, phi);
    torch::Tensor E_idx = get_k_neighbors(dist, static_cast<int64_t>(n_connections));
    torch::Tensor E = get_edge_features(dist, omega, theta, phi, E_idx);
    torch::Tensor str_mask = get_mask(E);

    E = replace_nan(E);
    V = replace_nan(V);

    torch::Tensor nodes = torch::zeros({static_cast<int>(n), static_cast<int>(max_ell), static_cast<int>(n_node_features)});
    torch::Tensor edges = torch::zeros({static_cast<int>(n), static_cast<int>(max_ell), static_cast<int>(n_connections), static_cast<int>(n_edge_features)});
    torch::Tensor connections = torch::zeros({static_cast<int>(n), static_cast<int>(max_ell), static_cast<int>(n_connections)}, torch::dtype(torch::kInt64));
    torch::Tensor edge_mask = torch::zeros({static_cast<int>(n), static_cast<int>(max_ell), static_cast<int>(n_connections), 1});

    for (core::Size i = 0; i < n; ++i) {
        core::Size nc = std::min(max_ell - 1, n_connections);

        // Use the pre-calculated structural features for each sequence
        nodes.slice(1, 0, max_ell).slice(0, i, i + 1) = V;
        edges.slice(1, 0, max_ell).slice(2, 0, nc).slice(0, i, i + 1) = E;
        connections.slice(1, 0, max_ell).slice(2, 0, nc).slice(0, i, i + 1) = E_idx;
        edge_mask.slice(1, 0, max_ell).slice(2, 0, nc).slice(0, i, i + 1) = str_mask;
    }

    return std::make_tuple(std::move(collated_seqs), std::move(nodes), std::move(edges), std::move(connections), std::move(edge_mask));
}


#endif //USE_TORCH
} //inverse_folding
} //protocols

