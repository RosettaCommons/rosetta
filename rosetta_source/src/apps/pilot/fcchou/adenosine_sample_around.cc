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
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>


#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>

#include <core/scoring/rna/RNA_Util.hh>

#include <protocols/idealize/idealize.hh>

#include <core/pose/Pose.hh>
#include <core/init.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.string.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/EnergyMap.hh> //for EnergyMap
#include <core/scoring/EnergyMap.fwd.hh> //for EnergyMap
#include <core/import_pose/import_pose.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace basic::options;

using utility::vector1;

typedef  numeric::xyzMatrix< Real > Matrix;


/////////////////////////////////////////////////////////////////////////////
//FCC: Adding Virtual res
void
add_virtual_res ( core::pose::Pose & pose, bool set_res_as_root = true ) {
	int nres = pose.total_residue();

	// if already rooted on virtual residue , return
	if ( pose.residue ( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) {
		std::cout << "addVirtualResAsRoot() called but pose is already rooted on a VRT residue ... continuing." << std::endl;
		return;
	}

	// attach virt res there
	bool fullatom = pose.is_fullatom();
	core::chemical::ResidueTypeSet const & residue_set = pose.residue_type ( 1 ).residue_type_set();
	core::chemical::ResidueTypeCAPs const & rsd_type_list ( residue_set.name3_map ( "VRT" ) );
	core::conformation::ResidueOP new_res ( core::conformation::ResidueFactory::create_residue ( *rsd_type_list[1] ) );
	pose.append_residue_by_jump ( *new_res , 1 );

	// make the virt atom the root
	if ( set_res_as_root ) {
		kinematics::FoldTree newF ( pose.fold_tree() );
		newF.reorder ( nres + 1 );
		pose.fold_tree ( newF );
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
// Rhiju -- rotate to my favorite frame. Base centroid is now at origin.
//         X points to N1 atom. Z points normal to base. Y is orthonormal and points towards Hoogsteen edge, I think.
void
rotate_into_nucleobase_frame( core::pose::Pose & pose ){

	using namespace core::conformation;
	using namespace core::scoring::rna;
	using namespace core::id;
	using namespace core::kinematics;

	// assuming pose has an RNA at residue 1 -- will rotate just that residue.
	Size const base_pos( 1 );
	Residue const & rsd = pose.residue( base_pos );

	Vector centroid = get_rna_base_centroid( rsd, true /*verbose*/ );
	Matrix M = get_rna_base_coordinate_system( rsd, centroid );
	Stub stub( M, centroid );

	for (Size i = 1; i <= rsd.natoms(); i++ ){
		//		Vector xyz_new = M.transpose() * ( rsd.xyz( i ) - centroid ); // it is either M or its transpose
		Vector xyz_new = stub.global2local( rsd.xyz(i) );
		pose.set_xyz( AtomID( i, base_pos ), xyz_new );
	}

	//	std::cout << pose.xyz( NamedAtomID( " N1 ", base_pos ) ) << std::endl;

}

/////////////////////////////////////////////////////////////////////////////////
void
methane_pair_score_test()
{
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;

	//////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	// Read in pose with two methane. "Z" = ligand. Note need flag:
	//         -extra_res_fa CH4.params -s two_methane.pdb
	pose::Pose pose;
	std::string infile  = option[ in ::file::s ][1];
	import_pose::pose_from_pdb( pose, *rsd_set, infile );

	rotate_into_nucleobase_frame( pose );
	pose.dump_pdb( "a_rotated.pdb" );

	add_virtual_res(pose);
	core::chemical::ResidueTypeSet const & residue_set = pose.residue_type ( 1 ).residue_type_set();
	core::chemical::ResidueTypeCAPs const & rsd_type_list ( residue_set.name3_map ( "CCC" ) );
	core::conformation::ResidueOP new_res ( core::conformation::ResidueFactory::create_residue ( *rsd_type_list[1] ) );
	pose.append_residue_by_jump ( *new_res , 2 );

	//////////////////////////////////////////////////
	// Set up fold tree -- "chain break" between two ligands, right?
	kinematics::FoldTree f( pose.total_residue() );
	Size start( 1 ), end( 2 );
	f.new_jump( start, end, start );
	// Uh, why not?
	f.set_jump_atoms( 1," C1 "," C1 ");
	pose.fold_tree( f );

	//////////////////////////////////////////////////////////////////
	// displace in z by 2.0 A... just checking coordinate system
	//This jump should code for no translation or rotation -- two_benzenes.pdb
	// has the two benzenes perfectly superimposed.
	kinematics::Jump jump( pose.jump( 1 ) );

	jump.set_translation( Vector( 5.0, 0.0, 0.0 ) );
	pose.set_jump( 1, jump );
	pose.dump_pdb( "shift_z.pdb" );


	//////////////////////////////////////////////////////////////////
	// OK, how about a score function?
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "vdw_test" );
	//scorefxn->set_weight( hack_elec, 1.0 );

	jump.set_translation( Vector( 1.3, 0.0, 0.0 ) );
	pose.set_jump( 1, jump );
	pose.dump_pdb( "test.pdb" );
	(*scorefxn)( pose );
	scorefxn->show( std::cout, pose );

	//////////////////////////////////////////////////////////////////
	// Now iterate over a bunch of translations and rotations.
	using namespace core::io::silent;
	SilentFileData silent_file_data;
	std::string const silent_file = option[ out::file::silent  ]();

	Size n( 0 );
	for (int i = 1; i <= 600; i++ ) {
			Real const x = i * 0.01;
			Real const y = 0.0;
			Real const z = 0.0;

			jump.set_translation( Vector( x, y, z ) ) ;
			pose.set_jump( 1, jump );

			(*scorefxn)( pose );

			n++;
			std::string out_tag( "S_"+lead_zero_string_of(n,4) );
			BinaryProteinSilentStruct s( pose,  out_tag);
			s.add_energy( "x", x );
			silent_file_data.write_silent_struct( s, silent_file, true /*write score only*/ );
	}
}

///////////////////////////////////////////////////////////////
void
my_main()
{

	methane_pair_score_test();

	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	core::init(argc, argv);
	my_main();
}
