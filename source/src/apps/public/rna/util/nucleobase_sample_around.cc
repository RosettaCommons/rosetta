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
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/BinarySilentStruct.hh>

#include <core/pose/variant_util.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <devel/init.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/sys_util.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/stepwise/modeler/util.hh> //has euler angle stuff.
#include <protocols/toolbox/sample_around/util.hh> //has euler angle stuff.

#include <protocols/viewer/viewers.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/sample_around.OptionKeys.gen.hh>

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
using Matrix = numeric::xyzMatrix<Real>;

// all helper functions moved to protocols/toolbox/sample_around/util.cc
using namespace protocols::toolbox::sample_around;

OPT_KEY( Boolean, sample_water )
OPT_KEY( Boolean, sample_phosphate )
OPT_KEY( Boolean, center_on_OP2 )
OPT_KEY( Boolean, sample_another_adenosine )
OPT_KEY( Boolean, just_z )
OPT_KEY( Boolean, just_xy )
OPT_KEY( Boolean, just_xz )
OPT_KEY( Boolean, just_yz )
OPT_KEY( Boolean, quick_score )
OPT_KEY( String, copy_adenosine_adenosine_file )
OPT_KEY( String, nucleobase )
OPT_KEY( Real, xyz_increment )
OPT_KEY( Real, xyz_size )

//////////////////////////////////////////////////////////////////
//
// To sample a 'carbon' probe atom:
// nucleobase_sample_around   [-s a_RNA.pdb]
//
// To sample a water
// nucleobase_sample_around   [-s a_RNA.pdb]  -sample_water
//
// To sample an adenosine
// nucleobase_sample_around   [-s a_RNA.pdb]  -sample_another_adenosine   -copy_adenosine_adenosine_file double_A_ready_set.pdb
//
// To sample an adenosine, reading in a starting adenosine-adenosine pairing conformation.
// nucleobase_sample_around   [-s a_RNA.pdb]  -sample_another_adenosine   -copy_adenosine_adenosine_file double_A_ready_set.pdb
//
// For water runs
//
//   -extra_res ~/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/water/TP3.params
//
// since included in rna types.
//
//
//
//////////////////////////////////////////////////////////////////
//
// Scanning best phosphate/base geometries
//  according to the Rosetta potential -- looking for
//  problems in H-bond computation that might
//  be reducing strength of, e.g., G-phosphate interactions
//  in sarcin-ricin loop.
//
//        -- rhiju, 2014
//
//////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////////////
void
nucleobase_probe_score_test()
{
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::id;

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
	pose.dump_pdb(out_prefix +  "a_rotated.pdb" );

	add_virtual_res(pose);
	core::chemical::ResidueTypeSet const & residue_set = *pose.residue_type_set_for_pose();

	core::conformation::ResidueOP new_res;

	bool const sample_water_ = option[ sample_water ]();
	bool const sample_phosphate_ = option[ sample_phosphate ]();
	bool const sample_another_adenosine_ = option[ sample_another_adenosine ]();

	if ( sample_another_adenosine_ || sample_phosphate_ ) add_another_virtual_res(pose); // this is the coordinate system for the next base.

	if ( sample_water_ ) {
		new_res = ( core::conformation::ResidueFactory::create_residue ( *(residue_set.get_representative_type_name3( "HOH" )) ) );
	} else if ( sample_another_adenosine_ ) {
		new_res = pose.residue(1).clone();
	} else if ( sample_phosphate_ ) {
		new_res = pose.residue(1).clone(); // tricky -- will add a packable phosphate and some other shenanigans later.
	} else {
		new_res = ( core::conformation::ResidueFactory::create_residue ( *(residue_set.get_representative_type_name3( " CZ" )) ) );
	}

	// doing some checks for hbonds & geom_sol are being calculated -- or not calculated.
	if ( sample_water_ || sample_another_adenosine_ ) {

		for ( Size n = 1; n <= new_res->natoms(); n++ ) {

			std::cout << "PROBE ATOM:  " << n << ' ' << new_res->atom_name( n ); // << ' ' << new_res->atom_base( n ) << ' ' << new_res->abase2( n ) << std::endl;
			if ( new_res->atom_base( n ) > 0 ) std::cout << "   base: " << new_res->atom_name( new_res->atom_base( n ) );
			if ( new_res->abase2(n) > 0 ) std::cout << "   base2: " <<  new_res->atom_name( new_res->abase2( n ) );
			std::cout << std::endl;

		}
	}

	pose.append_residue_by_jump ( *new_res ,  pose.size() );


	//////////////////////////////////////////////////
	// Set up fold tree -- "chain break" between two ligands, right?
	kinematics::FoldTree f ( pose.fold_tree() );
	std::string probe_atom_name = " C1 ";
	if ( sample_water_ ) probe_atom_name = " O  ";
	if ( sample_another_adenosine_ || sample_phosphate_ ) probe_atom_name = "ORIG"; // tricky -- the jump is from one coordinate system to the next one!
	f.set_jump_atoms( 2,"ORIG",probe_atom_name);
	if ( sample_phosphate_ ) f.set_jump_atoms( 3, "ORIG", " P  " ); // hope this works.
	pose.fold_tree( f );

	if ( sample_phosphate_ ) { // "ghost phosphate" at 5' end.
		Size n = 4;
		pose::remove_variant_type_from_pose_residue( pose, VIRTUAL_PHOSPHATE, n );
		pose::add_variant_type_to_pose_residue( pose, VIRTUAL_BASE,  n );
		pose::add_variant_type_to_pose_residue( pose, FIVE_PRIME_PHOSPHATE, n );
		rotate_into_phosphate_frame( pose, n, option[ center_on_OP2 ]()  );
	}

	std::cout << pose.annotated_sequence() << std::endl;

	if ( option[ copy_adenosine_adenosine_file ].user() ) {
		pose::Pose pose_reference;
		import_pose::pose_from_file( pose_reference, *rsd_set_op, option[copy_adenosine_adenosine_file]() , core::import_pose::PDB_file);
		rotate_into_nucleobase_frame( pose_reference );

		// copy over coordinates.
		Residue const & rsd_ref = pose_reference.residue( 2 );
		Size pos2( 4 ); // the sequence of the working pose is... nucleobase-virtual-virtual-nucleobase

		for ( Size i_ref = 1; i_ref <= rsd_ref.natoms(); i_ref++ ) {
			Size i = pose.residue( pos2 ).atom_index(   rsd_ref.atom_name( i_ref ) );
			pose.set_xyz( AtomID( i, pos2 ),  rsd_ref.xyz( i_ref ) );
		}
	}

	pose.dump_pdb(out_prefix +  "START.pdb" );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 800, 800 );

	//////////////////////////////////////////////////////////////////
	// displace in z by 2.0 A... just checking coordinate system
	//This jump should code for no translation or rotation -- two_benzenes.pdb
	// has the two benzenes perfectly superimposed.
	Size const probe_jump_num( 2 );
	kinematics::Jump jump( pose.jump( probe_jump_num ) );

	jump.set_translation( Vector( 5.0, 0.0, 0.0 ) );
	pose.set_jump( probe_jump_num, jump );
	pose.dump_pdb( out_prefix + "shift_x.pdb" );

	jump.set_translation( Vector( 0.0, 5.0, 0.0 ) );
	pose.set_jump( probe_jump_num, jump );
	pose.dump_pdb( out_prefix + "shift_y.pdb" );

	jump.set_translation( Vector( 0.0, 0.0, 5.0 ) );
	pose.set_jump( probe_jump_num, jump );
	pose.dump_pdb( out_prefix + "shift_z.pdb" );

	/// This is a code snippet to test if we are sampling water rotations properly -- could make this a little class,
	//  and then include within translation scan.
	// if ( sample_water_ )  sample_all_rotations_at_jump( pose, probe_jump_num );
	bool const sample_rotations = sample_water_ || sample_phosphate_;

	//////////////////////////////////////////////////////////////////
	// OK, how about a score function?
	ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ) {
		scorefxn = scoring::get_score_function();
	} else {
		scorefxn = ScoreFunctionFactory::create_score_function( "rna/denovo/rna_hires" );
		scorefxn->set_weight( rna_sugar_close, 0.0 ); //still computed with virtual sugar? weird.
	}

	// jump.set_translation( Vector( 3.75, 1.75, 1.5 ) );
	// pose.set_jump( probe_jump_num, jump );
	// sample_all_rotations_at_jump( pose, probe_jump_num, scorefxn );
	// pose.dump_pdb( "test.pdb" );
	// exit(0);
	//  (*scorefxn)( pose );
	// scorefxn->show( std::cout, pose );

	//////////////////////////////////////////////////////////////////
	// compute scores on a plane for now.

	Real const box_size = option[ xyz_size ]();
	Real const translation_increment = option[ xyz_increment ]();
	auto box_bins = int( box_size/translation_increment );
	bool const do_all = !option[ just_xy ]() && !option[ just_yz ]() && !option[ just_z ]() && !option[ just_xz ]();

	using namespace core::io::silent;
	SilentFileOptions opts; // initialized from the command line
	SilentFileData silent_file_data( opts );
	utility::io::ozstream out;

	if ( option[ just_z ]() || do_all ) {
		//////////////////////////////////////////////
		std::cout << "Doing Z scan..." << std::endl;
		out.open( out_prefix + "score_z.table" );
		Size count( 0 );
		std::string const silent_file( out_prefix + option[ out::file::silent ]() );
		for ( int i = -box_bins; i <= box_bins; ++i ) {
			Real const x = 0.0;
			Real const y = 0.0;
			Real const z = i * translation_increment;
			jump.set_translation( Vector( x, y, z ) ) ;
			pose.set_jump( probe_jump_num, jump );
			out << z << ' ' << centroid_dist( pose, option[ sample_another_adenosine ] )  << ' ' << do_scoring( pose, scorefxn, sample_rotations, probe_jump_num ) << std::endl;

			std::string const out_file_tag = "S_"+ObjexxFCL::lead_zero_string_of( 6, count );
			BinarySilentStruct s( opts, pose, out_file_tag );
			s.add_energy( "z", z );
			silent_file_data.write_silent_struct( s, silent_file, false /*score_only*/ );
		}
		out.close();
	}


	if ( option[ just_xy ]() || do_all ) {
		//////////////////////////////////////////////
		std::cout << "Doing XY scan... Z =  0.0" << std::endl;
		do_xy_scan( pose, scorefxn, out_prefix + "score_xy_0.table", 0.0, probe_jump_num, box_bins, translation_increment, sample_rotations );
		pose.dump_pdb( out_prefix + "best_xy.pdb");
		scorefxn->show( std::cout, pose );
	}


	if ( do_all ) {
		//////////////////////////////////////////////
		std::cout << "Doing XY scan... Z = +1.5" << std::endl;
		do_xy_scan( pose, scorefxn, out_prefix + "score_xy_1.5.table", 1.5, probe_jump_num, box_bins, translation_increment, sample_rotations );

		std::cout << "Doing XY scan... Z = +3.0" << std::endl;
		do_xy_scan( pose, scorefxn, out_prefix + "score_xy_3.table", 3.0, probe_jump_num, box_bins, translation_increment, sample_rotations );

		std::cout << "Doing XY scan... Z = +4.0" << std::endl;
		do_xy_scan( pose, scorefxn, out_prefix + "score_xy_4.table", 4.0, probe_jump_num, box_bins, translation_increment, sample_rotations );
		// Following are exactly the same as +1.0 and +3.0 when modeling nucleobase.
		//std::cout << "Doing XY scan... Z = -1.0" << std::endl;
		// do_xy_scan( pose, scorefxn, "score_para_0_table", 1.0, probe_jump_num, box_bins, translation_increment, sample_rotations );

		// std::cout << "Doing XY scan... Z = -3.0" << std::endl;
		// do_xy_scan( pose, scorefxn, "score_para_0_table", 3.0, probe_jump_num, box_bins, translation_increment, sample_rotations );
	}

	if ( option[ just_xz ]() || do_all ) {
		//////////////////////////////////////////////
		std::cout << "Doing XZ scan..." << std::endl;
		out.open( out_prefix + "score_xz.table" );
		for ( int i = -box_bins; i <= box_bins; ++i ) {
			for ( int j = -box_bins; j <= box_bins; ++j ) {
				Real const x = j * translation_increment;
				Real const z = i * translation_increment;
				Real const y = 0.0;
				jump.set_translation( Vector( x, y, z ) ) ;
				pose.set_jump( probe_jump_num, jump );
				out << do_scoring( pose, scorefxn, sample_rotations, probe_jump_num ) << ' ' ;
			}
			out << std::endl;
		}
		out.close();
	}

	if ( option[ just_yz ]() || do_all ) {
		//////////////////////////////////////////////
		std::cout << "Doing YZ scan..." << std::endl;
		out.open( out_prefix + "score_yz.table" );
		for ( int i = -box_bins; i <= box_bins; ++i ) {
			for ( int j = -box_bins; j <= box_bins; ++j ) {
				Real const y = j * translation_increment;
				Real const z = i * translation_increment;
				Real const x = 0.0;
				jump.set_translation( Vector( x, y, z ) ) ;
				pose.set_jump( probe_jump_num, jump );
				out << do_scoring( pose, scorefxn, sample_rotations, probe_jump_num ) << ' ' ;
			}
			out << std::endl;
		}
		out.close();
	}

}

///////////////////////////////////////////////////////////////
void
quick_score_test(){

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::id;

	//////////////////////////////////////////////////
	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	// Read in pose with two methane. "Z" = ligand. Note need flag:
	//         -extra_res_fa CH4.params -s two_methane.pdb
	pose::Pose pose;
	std::string infile  = option[ in ::file::s ][1];
	import_pose::pose_from_file( pose, *rsd_set, infile , core::import_pose::PDB_file);

	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( pose.residue(i).is_RNA() ) {
			pose::add_variant_type_to_pose_residue( pose, VIRTUAL_PHOSPHATE, i );
			pose::add_variant_type_to_pose_residue( pose, VIRTUAL_RIBOSE, i );
		}
	}

	ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ) {
		scorefxn = scoring::get_score_function();
	} else {
		scorefxn = ScoreFunctionFactory::create_score_function( "rna/denovo/rna_hires" );
		scorefxn->set_weight( rna_sugar_close, 0.0 ); //still computed with virtual sugar? weird.
	}

	(*scorefxn)( pose );
	scorefxn->show( std::cout, pose );

}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	if ( option[quick_score]() ) {
		quick_score_test();
	} else {
		nucleobase_probe_score_test();
	}

	protocols::viewer::clear_conformation_viewers();

	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		NEW_OPT( sample_water, "use a water probe instead of carbon", false );
		NEW_OPT( sample_another_adenosine, "sample another adenosine as the 'probe'", false );
		NEW_OPT( sample_phosphate, "use a phosphate (mimicking an RNA terminal phosphate) as the 'probe'", false );
		NEW_OPT( center_on_OP2, "Define coordinate frame on phosphate to be centered on OP2, not on P", false );
		NEW_OPT( just_z, "Just scan z at x,y=0", false );
		NEW_OPT( just_xy, "Just scan x, y at z=0", false );
		NEW_OPT( just_xz, "Just scan x, z at y=0", false );
		NEW_OPT( just_yz, "Just scan y, z at x=0", false );
		NEW_OPT( quick_score, "alternative mode for checking geom_sol, etc.", false );
		NEW_OPT( copy_adenosine_adenosine_file, "get rigid body relation between two adenosines from file", "" );
		NEW_OPT( xyz_increment, "input parameter", 0.2 );
		NEW_OPT( xyz_size, "input parameter", 10.0 );
		NEW_OPT( nucleobase, "nucleobase to sample around", "a" );
		option.add_relevant( sample_around::alpha_increment );
		option.add_relevant( sample_around::cosbeta_increment );
		option.add_relevant( sample_around::gamma_increment );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
