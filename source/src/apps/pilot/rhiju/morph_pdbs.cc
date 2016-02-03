// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief

// libRosetta headers
#include <core/types.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/copydofs/CopyDofsInfo.hh>
#include <core/init/init.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/simple_moves/CopyDofMover.hh>
#include <protocols/viewer/viewers.hh>
#include <numeric/angle.functions.hh>
#include <utility/string_util.hh>
#include <ObjexxFCL/string.functions.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <utility/excn/Exceptions.hh>

using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace basic::options;

using utility::vector1;
typedef  numeric::xyzMatrix< Real > Matrix;

OPT_KEY( Integer, nsamples )

//////////////////////////////////////////////////////////////////
//
// Morph one pdb into others, and output intermediates for making a
//  nice movie in Pymol. Based on internal coordinates, and, frankly,
//  does not work well.
//
// I think a better solution would be to numerically optimize a trajectory
//  between start and finish that is smooth but also in some way
//  minimizes the energy along the path (incl. strong contraints on bond
//  angles and distances). There's probably a fast, canned
//  algorithm for doing this in the motion planning literature, or
//  from simple Lagrangian mechanics simulation.
//
//        -- rhiju, 2015
//
//////////////////////////////////////////////////////////////////

typedef utility::vector1< std::pair< core::id::DOF_ID, core::Real > > DofsInfo;
void
morph_by_internal_coords(
	pose::Pose & pose /* for viewing */,
	pose::PoseCOP pose_start,
	pose::PoseCOP pose_end,
	utility::vector1< pose::PoseCOP > & morph_poses,
	Size const num_samples = 24 )
{
	using namespace core::id;
	using namespace core::pose;
	using namespace protocols::simple_moves;

	//	morph_poses.push_back( pose_start );

	std::map< Size, Size > res_map;
	for ( Size n = 1; n <= pose_start->total_residue(); n++ ) res_map[ n ] = n;
	CopyDofMover copy_dof_mover_end(   *pose_start, res_map );
	CopyDofMover copy_dof_mover_start( *pose_end,   res_map );

	// grab info on DOFs.
	pose = *pose_start;
	copy_dof_mover_start.apply(   pose );
	DofsInfo copy_dofs_info_start( copy_dof_mover_start.copy_dofs_info( pose ).dofs_info() );;

	copy_dof_mover_end.apply( pose );
	DofsInfo copy_dofs_info_end  ( copy_dof_mover_end.copy_dofs_info( pose ).dofs_info() );

	for ( Size n = 1; n <= num_samples; n++ ) {

		pose::copydofs::CopyDofsInfo copy_dofs_info;
		for ( Size m = 1; m <= copy_dofs_info_start.size(); m++ ) {

			DOF_ID const dof_id = copy_dofs_info_start[ m ].first;
			runtime_assert( dof_id == copy_dofs_info_end[ m ].first );

			Real dof_change = copy_dofs_info_end[ m ].second - copy_dofs_info_start[ m ].second;
			if ( dof_id.type() == PHI ) dof_change = numeric::principal_angle( dof_change ); // radians or degrees?
			Real const new_dof    = copy_dofs_info_start[ m ].second + dof_change * ( Real( n ) / num_samples );

			copy_dofs_info.push_back( std::make_pair( dof_id, new_dof) );
		}

		CopyDofMover copy_dof_mover( *pose_end, res_map );
		copy_dof_mover.set_copy_dofs_info( *pose_end, copy_dofs_info );
		copy_dof_mover.apply( pose );
		morph_poses.push_back( pose.clone() );
	}

}

/////////////////////////////////////////////////////////////////////////////////
void
morph_pdbs_test()
{
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::id;
	using namespace core::pose;

	//////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	utility::vector1< pose::PoseOP > poses;
	core::chemical::ResidueTypeSetCOP rsd_set_op( rsd_set );
	utility::vector1< std::string > infiles = option[ in::file::s ]();
	runtime_assert( infiles.size() > 0 );

	Pose pose;
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );
	for ( Size n = 1; n <= infiles.size(); n++ ) {
		std::string infile  = infiles[ n ];
		import_pose::pose_from_file( pose, *rsd_set_op, infile , core::import_pose::PDB_file);
		poses.push_back( pose.clone() );
	}

	utility::vector1< pose::PoseCOP > morph_poses;
	morph_poses.push_back( poses[ 1 ] );
	for ( Size n = 2; n <= infiles.size(); n++ ) {
		PoseCOP pose_start = poses[ n - 1 ];
		PoseCOP pose_end   = poses[ n ];
		morph_by_internal_coords( pose /* for viewing*/, pose_start, pose_end, morph_poses, option[ nsamples ]() );
	}

	runtime_assert( infiles[ 1 ].find( ".pdb" ) != std::string::npos );
	std::string tag = utility::replace_in( infiles[ 1 ], ".pdb", "" ) ;
	if ( option[ out::output_tag ].user() ) tag = option[ out::output_tag ]();
	for ( Size n = 1; n <= morph_poses.size(); n++ ) {
		std::string const output_pdb_file = tag + ".MORPH." + ObjexxFCL::lead_zero_string_of( n, 3 ) + ".pdb";
		morph_poses[ n ]->dump_pdb( output_pdb_file );
		std::cout << "Outputted: " << output_pdb_file << std::endl;
	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	morph_pdbs_test();

	option.add_relevant( in::file::s );
	option.add_relevant( out::output );

	protocols::viewer::clear_conformation_viewers();

	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		NEW_OPT( nsamples,  "number of intermediates between each pair of poses", 24 );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
