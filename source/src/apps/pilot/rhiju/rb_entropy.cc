// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

// libRosetta headers
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/rna/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>

#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/init/init.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/basic_sys_util.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/stepwise/modeler/util.hh> //has euler angle stuff.
#include <protocols/toolbox/sample_around/util.hh> //has euler angle stuff.
#include <protocols/toolbox/rigid_body/util.hh>

#include <protocols/viewer/viewers.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/sample_around.OptionKeys.gen.hh>
#include <basic/Tracer.hh>


#include <numeric/random/random_xyz.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>

#include <utility/stream_util.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <utility/excn/Exceptions.hh>

static THREAD_LOCAL basic::Tracer TR( "rb_entropy" );

using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace basic::options;

using utility::vector1;
typedef  numeric::xyzMatrix< Real > Matrix;

// all helper functions moved to protocols/toolbox/sample_around/util.cc
using namespace protocols::toolbox::sample_around;

OPT_KEY( String, nucleobase )
OPT_KEY( Boolean, twodimensional )
OPT_KEY( Real, xyz_size )

//////////////////////////////////////////////////////////////////
//
// Can we calculate free energy of forming a base pair
//  to start a helix ("init" in the nearest neighbor rules)?
//
// Basic tests of:
//
//  1. Fast RMSD calculation
//  2. Analytical calculation of RMSD distributions at "1 M" reference state.
//  3. Movers to install in thermal_sampler.
//
//        -- rhiju, 2016
//
//////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
void
rb_entropy_test()
{
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::id;
	using namespace core::kinematics;
	using numeric::random::rg;

	//////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	// Read in pose with two methane. "Z" = ligand. Note need flag:
	//         -extra_res_fa CH4.params -s two_methane.pdb
	pose::Pose pose;
	core::chemical::ResidueTypeSetCOP rsd_set_op( rsd_set );
	if ( option[ in::file::s ].user() ) {
		std::string infile  = option[ in::file::s ][1];
		import_pose::pose_from_file( pose, *rsd_set_op, infile , core::import_pose::PDB_file);
	} else {
		std::string const sequence = option[ nucleobase ]();
		runtime_assert( sequence.size() == 1 );
		make_pose_from_sequence( pose, sequence, *rsd_set_op, false /*auto_termini*/ );
		core::pose::rna::apply_Aform_torsions( pose, 1 );
	}

	pose::add_variant_type_to_pose_residue( pose, VIRTUAL_PHOSPHATE, 1 );
	pose::add_variant_type_to_pose_residue( pose, VIRTUAL_RIBOSE, 1 );

	rotate_into_nucleobase_frame( pose );
	std::string const out_prefix = option[ out::prefix]();
	pose.dump_pdb( std::string(option[ out::path::path]()) + "/" + out_prefix +  "a_rotated.pdb" );

	add_virtual_res( pose );

	kinematics::FoldTree f ( pose.fold_tree() );
	std::cout << pose.annotated_sequence() << std::endl;
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 800, 800 );


	Size cycles( 1000000 );
	core::pose::Pose start_pose = pose;
	Vector start_base_centroid = core::chemical::rna::get_rna_base_centroid( start_pose.residue( 1 ) );
	Size const probe_jump_num( 1 );
	Stub const upstream_stub = pose.conformation().upstream_jump_stub( probe_jump_num );
	kinematics::Jump start_jump( start_pose.jump( probe_jump_num ) );
	TR << "original centroid: " << start_base_centroid.x() << " " << start_base_centroid.y() << " " << start_base_centroid.z() << std::endl;

	// to compute moments of inertia in MATLAB
	Residue const & rsd( pose.residue(1) );
	utility::io::ozstream out_xyz;
	std::string filename_xyz( std::string( option[ out::path::path]() ) + "/xyz.txt" );
	out_xyz.open( filename_xyz );
	for ( Size i = rsd.first_sidechain_atom() + 1; i <= rsd.nheavyatoms(); ++i ) { //rsd.first_sidechain_atom()+1 to not include the O2prime oxygen.
		TR << rsd.atom_name( i ) << std::endl;
		if ( rsd.is_virtual( i ) ) continue;
		out_xyz << rsd.xyz(i).x() << " " << rsd.xyz(i).y() << " " << rsd.xyz(i).z() << std::endl;
	}
	out_xyz.close();
	TR << "Outputted xyz values into " << filename_xyz << std::endl;

	Matrix M;
	Real const & box_size = option[ xyz_size ]();

	utility::io::ozstream out;
	std::string filename( std::string( option[ out::path::path]() ) + "/rmsd.txt" );
	out.open( filename );
	for ( Size n = 1; n <= cycles; n++ ) {

		kinematics::Jump jump( start_jump );

		// rotation
		if ( option[ twodimensional]() ) {
			// 2D rotations, just for testing...
			M = rotation_matrix( Vector( 0.0, 0.0, 1.0 ), 2 * numeric::constants::d::pi * rg().uniform() );
		} else {
		// 3D rotation -- thank you will sheffler & quaternions
			M = numeric::random::random_rotation();
		}
		jump.rotation_by_matrix( upstream_stub, start_base_centroid, M );

		// translation
		Vector random_translation( 2.0 * rg().uniform() - 1.0,
															 2.0 * rg().uniform() - 1.0,
															 2.0 * rg().uniform() - 1.0 );
		random_translation *= box_size;
		jump.set_translation( jump.get_translation() +  random_translation );

		pose.set_jump( probe_jump_num, jump );

		Residue const & rsd1( pose.residue(1) ), rsd2( start_pose.residue(1) );
		Real rmsd( 0.0 ); Size numatoms( 0 );
		for ( Size i = rsd1.first_sidechain_atom() + 1; i <= rsd1.nheavyatoms(); ++i ) { //rsd.first_sidechain_atom()+1 to not include the O2prime oxygen.
			if ( rsd1.is_virtual( i ) ) continue;
			Vector dist = ( rsd1.xyz(i) - rsd2.xyz(i) );
			rmsd += dist.length_squared();
			numatoms++;
		}
		rmsd = sqrt( rmsd/numatoms );
		out << rmsd << std::endl;

		// tests that base_centroid is in same place:
		//	base_centroid = get_rna_base_centroid( pose.residue( 1 ) );
		//	TR << "new centroid:     " << base_centroid.x() << " " << base_centroid.y() << " " << base_centroid.z() << std::endl;

		if ( n <= 10 ) pose.dump_pdb( std::string( option[ out::path::path]() ) + "/test"+ObjexxFCL::lead_zero_string_of( n, 4 )+".pdb" );

	}
	out.close();
	TR << "Outputted " << cycles << " RMSD values into " << filename << std::endl;

}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	rb_entropy_test();
	protocols::viewer::clear_conformation_viewers();

	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		NEW_OPT( nucleobase, "nucleobase to sample around", "a" );
		NEW_OPT( twodimensional, "choose random rotation about Z-axis (2D)", false );
		NEW_OPT( xyz_size, "box half-diameter (max in x,y,z)", 1.0 );

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
