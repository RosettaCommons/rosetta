// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Build aggregate and per-residue motif projection clouds and write PDB visualizations.

#include <devel/init.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>

#include <protocols/motifs/MotifCloud.hh>
#include <protocols/motifs/MotifLibrary.hh>
#include <protocols/motifs/motif_utils.hh>

#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

#include <iostream>
#include <set>
#include <string>

int
main( int argc, char * argv [] )
{
	try {
		devel::init(argc, argv);

		basic::Tracer tracer("apps.public.motif_ligand_discovery.motif_cloud");
		using basic::options::OptionKeys::in::file::s;
		using basic::options::OptionKeys::motifs::list_motifs;
		using basic::options::OptionKeys::motifs::motif_cloud_output_prefix;
		using basic::options::OptionKeys::motifs::motif_cloud_residue_indices;
		using basic::options::OptionKeys::motifs::motif_filename;
		using basic::options::option;

		if ( !option[s].user() || option[s]().size() != 1 ) {
			utility_exit_with_message("motif_cloud requires exactly one input PDB via -in:file:s.");
		}
		if ( !option[motif_cloud_residue_indices].user() ) {
			utility_exit_with_message("motif_cloud requires -motifs:motif_cloud_residue_indices with one or more 1-based pose residue indices.");
		}
		if ( !option[motif_filename].user() && !option[list_motifs].user() ) {
			utility_exit_with_message("motif_cloud requires a motif library via -motifs:motif_filename or -motifs:list_motifs.");
		}

		core::pose::PoseOP pose(new core::pose::Pose);
		core::import_pose::pose_from_file(*pose, option[s]()[1], core::import_pose::PDB_file);

		utility::vector1<int> const requested_indices = option[motif_cloud_residue_indices]();
		utility::vector1<core::Size> residue_indices;
		std::set<int> unique_indices;
		for ( int const index : requested_indices ) {
			if ( index < 1 || static_cast<core::Size>(index) > pose->size() ) {
				utility_exit_with_message("MotifCloud residue index " + std::to_string(index) +
					" is outside the input pose range 1-" + std::to_string(pose->size()) + ".");
			}
			if ( !unique_indices.insert(index).second ) {
				utility_exit_with_message("MotifCloud residue index " + std::to_string(index) + " was specified more than once.");
			}
			residue_indices.push_back(static_cast<core::Size>(index));
		}
		if ( residue_indices.empty() ) {
			utility_exit_with_message("motif_cloud requires at least one source residue index.");
		}

		utility::vector1<std::string> ligand_atom_names;
		protocols::motifs::MotifLibrary motif_library(
			protocols::motifs::get_LigandMotifLibrary_user(true, ligand_atom_names)
		);
		protocols::motifs::MotifCOPs const motifs = motif_library.library();
		if ( motifs.empty() ) {
			utility_exit_with_message("The supplied ligand motif library contains no usable motifs.");
		}

		protocols::motifs::MotifCloud cloud(pose, motifs, residue_indices);
		std::string const output_prefix = option[motif_cloud_output_prefix]();
		cloud.write_visualization_pdb(cloud.aggregate_matrix(), output_prefix + "_aggregate.pdb");

		protocols::motifs::MotifCloud::SparseMotifMatrix const empty_matrix;
		protocols::motifs::MotifCloud::PerResidueMotifMatrices const & per_residue =
			cloud.per_residue_matrices();
		for ( core::Size const index : residue_indices ) {
			auto const matrix_iter = per_residue.find(index);
			protocols::motifs::MotifCloud::SparseMotifMatrix const & matrix =
				matrix_iter == per_residue.end() ? empty_matrix : matrix_iter->second;
			cloud.write_visualization_pdb(
				matrix,
				output_prefix + "_residue_" + std::to_string(index) + ".pdb"
			);
		}

		tracer.Info << "Wrote aggregate and " << residue_indices.size()
			<< " per-residue MotifCloud PDB visualizations using a shared pose XYZ coordinate frame."
			<< std::endl;
	} catch (utility::excn::Exception const & exception) {
		exception.display();
		return -1;
	}
	return 0;
}
