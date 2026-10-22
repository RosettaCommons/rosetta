// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Sparse voxelized projections of motif partner anchor geometry around selected pose residues.

#include <protocols/motifs/MotifCloud.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/motif_utils.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <set>
#include <string>
#include <utility>

namespace protocols {
namespace motifs {

static basic::Tracer ms_tr( "protocols.motifs.MotifCloud" );

MotifCloud::~MotifCloud() = default;

MotifCloud::MotifCloud(
	core::pose::PoseOP in_pose,
	protocols::motifs::MotifCOPs motif_library,
	utility::vector1<core::Size> const & source_residue_indices
) :
	working_pose_(std::move(in_pose)),
	motif_library_(std::move(motif_library))
{
	generate_fa_standard_atom_type_set();
	project_motifs_from_residues(source_residue_indices);
}

void MotifCloud::generate_fa_standard_atom_type_set()
{
	atom_type_set_ = core::chemical::ChemicalManager::get_instance()->atom_type_set("fa_standard");
}

void MotifCloud::write_visualization_pdb(
	SparseMotifMatrix const & matrix,
	std::string const & filename
) const
{
	core::Size maximum_voxel_count = 0;
	for ( auto const & voxel_entry : matrix ) {
		core::Size voxel_count = 0;
		for ( auto const & atom_type_count : voxel_entry.second ) {
			voxel_count += atom_type_count.second;
		}
		maximum_voxel_count = std::max(maximum_voxel_count, voxel_count);
	}

	std::ofstream output(filename.c_str());
	if ( !output ) {
		utility_exit_with_message("Could not open MotifCloud visualization output file: " + filename);
	}

	output << "REMARK MotifCloud voxel density visualization\n";
	output << "REMARK Each HETATM marks one occupied 1-Angstrom voxel at its center.\n";
	output << "REMARK Motif partner anchors 1 and 3 use idealized 1.4-Angstrom, 120-degree geometry.\n";
	output << "REMARK B-factor is voxel atom count as a percentage of the maximum voxel count in this cloud.\n";
	output << "REMARK Occupancy stores raw projected anchor count, capped at 999.99 by the PDB field width.\n";
	output << "REMARK Density ranges in percent of maximum: I 0-12.5, N 12.5-25, F 25-37.5, CL 37.5-50, S 50-62.5, P 62.5-75, O 75-87.5, BR 87.5-100.\n";
	std::array<std::string, 8> const density_bin_elements = {{
		"I", "N", "F", "CL", "S", "P", "O", "BR"
	}};

	core::Size serial = 0;
	for ( auto const & voxel_entry : matrix ) {
		core::Size voxel_count = 0;
		for ( auto const & atom_type_count : voxel_entry.second ) {
			voxel_count += atom_type_count.second;
		}
		if ( voxel_count == 0 ) {
			continue;
		}
		if ( serial >= 99999 ) {
			utility_exit_with_message("MotifCloud visualization exceeds the PDB 99,999 atom serial limit.");
		}
		++serial;

		core::Real const density_percent = maximum_voxel_count == 0 ? 0.0 :
			100.0 * static_cast<core::Real>(voxel_count) / static_cast<core::Real>(maximum_voxel_count);
		core::Size const density_bin_index = std::min<core::Size>(
			static_cast<core::Size>(density_percent / 12.5),
			density_bin_elements.size() - 1
		);
		std::string const & element = density_bin_elements[density_bin_index];
		std::string const atom_name = element.size() == 1 ? " " + element + "  " : " " + element + " ";
		core::Real const occupancy_count = std::min<core::Real>(static_cast<core::Real>(voxel_count), 999.99);
		output << "HETATM" << std::setw(5) << serial << ' ' << atom_name << " VXL X"
			<< std::setw(4) << 1 << "    "
			<< std::fixed << std::setprecision(3)
			<< std::setw(8) << voxel_entry.first[0] + 0.5
			<< std::setw(8) << voxel_entry.first[1] + 0.5
			<< std::setw(8) << voxel_entry.first[2] + 0.5
			<< std::setprecision(2) << std::setw(6) << occupancy_count
			<< std::setw(6) << density_percent << "          " << std::setw(2) << element << " \n";
	}
	output << "END\n";
	output.close();
	if ( !output ) {
		utility_exit_with_message("Failed while writing MotifCloud visualization output file: " + filename);
	}
}

void MotifCloud::project_motifs_from_residues(
	utility::vector1<core::Size> const & source_residue_indices
)
{
	if ( working_pose_ == nullptr ) {
		ms_tr.Error << "Cannot project motifs without an input pose." << std::endl;
		return;
	}

	utility::vector1<core::Size> selected_residues = source_residue_indices;
	if ( selected_residues.empty() ) {
		for ( core::Size residue_index = 1; residue_index <= working_pose_->size(); ++residue_index ) {
			selected_residues.push_back(residue_index);
		}
	}

	std::set<core::Size> seen_residue_indices;
	for ( core::Size source_residue_index : selected_residues ) {
		if ( source_residue_index < 1 || source_residue_index > working_pose_->size() ) {
			ms_tr.Warning << "Skipping invalid MotifCloud source residue index " << source_residue_index
				<< " for pose size " << working_pose_->size() << "." << std::endl;
			continue;
		}
		if ( !seen_residue_indices.insert(source_residue_index).second ) {
			ms_tr.Warning << "Skipping duplicate MotifCloud source residue index " << source_residue_index << "." << std::endl;
			continue;
		}

		core::conformation::Residue const & source_residue = working_pose_->residue(source_residue_index);
		std::string const source_residue_name = source_residue.name3();
		auto motif_library_iter = motif_library_map_by_residue_.find(source_residue_name);
		if ( motif_library_iter == motif_library_map_by_residue_.end() ) {
			std::string residue_name_for_lookup = source_residue_name;
			motif_library_map_by_residue_[source_residue_name] =
				get_motif_sublibrary_by_aa(motif_library_, residue_name_for_lookup);
			motif_library_iter = motif_library_map_by_residue_.find(source_residue_name);
		}
		if ( motif_library_iter->second.empty() ) {
			ms_tr.Warning << "No motifs found for source residue " << source_residue_name
				<< " at pose index " << source_residue_index << "." << std::endl;
			continue;
		}

		SparseMotifMatrix & residue_matrix = per_residue_matrices_[source_residue_index];
		for ( protocols::motifs::MotifCOP const & motif : motif_library_iter->second ) {
			if ( motif == nullptr ) {
				continue;
			}

			if ( motif->res1_atom1_index() == 0 || motif->res1_atom1_index() > source_residue.natoms() ||
					motif->res1_atom2_index() == 0 || motif->res1_atom2_index() > source_residue.natoms() ||
					motif->res1_atom3_index() == 0 || motif->res1_atom3_index() > source_residue.natoms() ) {
				ms_tr.Warning << "Skipping motif with invalid source anchor atom index for residue "
					<< source_residue_name << " at pose index " << source_residue_index << "." << std::endl;
				continue;
			}

			std::array<core::Size, 3> partner_atom_type_indices;
			std::array<std::string, 3> const partner_atom_type_names = {{
				motif->res2_atom1_name(), motif->res2_atom2_name(), motif->res2_atom3_name()
			}};
			bool partner_atom_types_valid = true;
			for ( core::Size anchor_index = 0; anchor_index < partner_atom_type_names.size(); ++anchor_index ) {
				if ( !atom_type_set_->has_atom(partner_atom_type_names[anchor_index]) ) {
					ms_tr.Warning << "Skipping motif with invalid partner anchor atom type "
						<< partner_atom_type_names[anchor_index] << " for source residue "
						<< source_residue_name << "." << std::endl;
					partner_atom_types_valid = false;
					break;
				}
				partner_atom_type_indices[anchor_index] =
					atom_type_set_->atom_type_index(partner_atom_type_names[anchor_index]);
			}
			if ( !partner_atom_types_valid ) {
				continue;
			}

			core::kinematics::Stub const source_stub(
				source_residue.xyz(motif->res1_atom2_index()),
				source_residue.xyz(motif->res1_atom1_index()),
				source_residue.xyz(motif->res1_atom2_index()),
				source_residue.xyz(motif->res1_atom3_index())
			);
			core::kinematics::Stub projected_partner_stub;
			motif->forward_jump().make_jump(source_stub, projected_partner_stub);

			core::Real constexpr approximate_bond_length = 1.4;
			core::Real constexpr approximate_bond_angle = 2.0 * 3.14159265358979323846 / 3.0;
			numeric::xyzVector<core::Real> const & frame_x = projected_partner_stub.M.col_x();
			numeric::xyzVector<core::Real> const & frame_y = projected_partner_stub.M.col_y();
			std::array<numeric::xyzVector<core::Real>, 3> const projected_anchor_xyz = {{
				projected_partner_stub.v + approximate_bond_length * frame_x,
				projected_partner_stub.v,
				projected_partner_stub.v + approximate_bond_length *
				(std::cos(approximate_bond_angle) * frame_x + std::sin(approximate_bond_angle) * frame_y)
			}};

			for ( core::Size anchor_index = 0; anchor_index < projected_anchor_xyz.size(); ++anchor_index ) {
				numeric::xyzVector<core::Real> const & xyz = projected_anchor_xyz[anchor_index];
				VoxelIndex const voxel = {{
					static_cast<int>(std::floor(xyz.x())),
					static_cast<int>(std::floor(xyz.y())),
					static_cast<int>(std::floor(xyz.z()))
				}};
				++residue_matrix[voxel][partner_atom_type_indices[anchor_index]];
				++aggregate_motif_matrix_[voxel][partner_atom_type_indices[anchor_index]];
			}
		}
	}
}

} // namespace motifs
} // namespace protocols
