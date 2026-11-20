// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/public/motif_ligand_discovery/protein_grid_binding_pocket.cc
/// @brief  Export a ProteinGrid binding-pocket representation without ligand discovery.
///
/// Required options:
///   -s <protein.pdb>
///   -binding_pocket_center_sf <x,y,z>
///   -binding_pocket_dimensions_sf <x,y,z>
///
/// The dimensions are half-widths in angstroms. Use -out:prefix to choose the output
/// filename prefix; the default is "binding_pocket". Empty-space and whole-matrix
/// output can be controlled with -use_empty_space_in_space_fill_pdb and
/// -use_whole_matrix_in_space_fill_pdb. Use -protein_grid:output_empty_space_only
/// to output only pocket empty-space atoms (oxygen).

#include <devel/init.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/motifs.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/protein_grid.OptionKeys.gen.hh>
#include <basic/options/util.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <numeric/xyzVector.hh>

#include <protocols/protein_grid/ProteinGrid.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

#include <string>

static basic::Tracer TR( "apps.public.protein_grid_binding_pocket" );

int
main( int argc, char * argv [] )
{
	try {
		devel::init( argc, argv );

		using basic::options::option;
		using namespace basic::options::OptionKeys;

		utility::vector1<std::string> const pdb_files = basic::options::start_files();
		if ( pdb_files.size() != 1 ) {
			TR.Error << "Provide exactly one input PDB with -s." << std::endl;
			return -1;
		}

		if ( !option[ motifs::binding_pocket_center_sf ].user() ||
				!option[ motifs::binding_pocket_dimensions_sf ].user() ) {
			TR.Error << "Both -binding_pocket_center_sf and -binding_pocket_dimensions_sf are required."
				<< std::endl;
			return -1;
		}

		utility::vector1<int> const center_values = option[ motifs::binding_pocket_center_sf ]();
		utility::vector1<int> const dimension_values = option[ motifs::binding_pocket_dimensions_sf ]();
		if ( center_values.size() != 3 || dimension_values.size() != 3 ) {
			TR.Error << "The binding-pocket center and dimensions must each contain exactly three "
				<< "comma-separated x,y,z values." << std::endl;
			return -1;
		}

		utility::vector1<core::Size> dimensions( 3, 0 );
		for ( core::Size axis = 1; axis <= 3; ++axis ) {
			if ( dimension_values[axis] < 0 ) {
				TR.Error << "Binding-pocket dimensions must be non-negative half-widths in angstroms."
					<< std::endl;
				return -1;
			}
			dimensions[axis] = static_cast<core::Size>( dimension_values[axis] );
		}

		numeric::xyzVector<int> const center(
			center_values[1], center_values[2], center_values[3]
		);

		core::pose::PoseOP pose( new core::pose::Pose );
		core::import_pose::pose_from_file( *pose, pdb_files[1], core::import_pose::PDB_file );
		if ( pose->size() == 0 ) {
			TR.Error << "The input PDB did not contain any residues: " << pdb_files[1] << std::endl;
			return -1;
		}

		core::Real const resolution = option[ motifs::resolution_scale_factor_float ]();
		protocols::protein_grid::ProteinGrid grid( pose, resolution, center, dimensions );
		grid.project_lj_radii();
		grid.set_print_whole_matrix( option[ motifs::use_whole_matrix_in_space_fill_pdb ]() );
		grid.set_print_empty_space( option[ motifs::use_empty_space_in_space_fill_pdb ]() );
		grid.set_print_empty_space_only( option[ protein_grid::output_empty_space_only ]() );

		std::string output_prefix = "binding_pocket";
		if ( option[ out::prefix ].user() ) {
			output_prefix = option[ out::prefix ]();
		}
		grid.export_protein_matrix_to_pdb( output_prefix );
		TR << "Exported binding-pocket matrix for " << pdb_files[1]
			<< " with center (" << center.x() << ", " << center.y() << ", " << center.z()
			<< ") and half-dimensions (" << dimensions[1] << ", " << dimensions[2]
			<< ", " << dimensions[3] << ")." << std::endl;
	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
