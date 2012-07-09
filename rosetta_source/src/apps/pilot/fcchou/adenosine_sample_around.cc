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
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/scoring/rna/RNA_Util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>

#include <protocols/idealize/idealize.hh>
#include <protocols/swa/StepWiseUtil.cc>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/init.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
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

OPT_KEY( Boolean, sample_water )

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

	// assuming pose has an RNA at residue 1 -- will rotate just that residue.
	Size const base_pos( 1 );
	Residue const & rsd = pose.residue( base_pos );

	Vector centroid = get_rna_base_centroid( rsd, true /*verbose*/ );
	Matrix M = get_rna_base_coordinate_system( rsd, centroid );
	kinematics::Stub stub( M, centroid );

	for (Size i = 1; i <= rsd.natoms(); i++ ){
		Vector xyz_new = stub.global2local( rsd.xyz( i ) ); // it is either this or M-inverse.
		pose.set_xyz( AtomID( i, base_pos ), xyz_new );
	}

}


///////////////////////////////////////////////////////////////////////////////////////////
// This is imported from protocols/swa/RigidBodySampler.cco
void
sample_all_rotations_at_jump( pose::Pose & pose, Size const num_jump ){

	Real alpha_, alpha_min_( 0 ), alpha_max_( 180.0 ), alpha_increment_( 60.0 );
	Real beta_, cosbeta_min_( -1.0 ), cosbeta_max_( 1.0 ), cosbeta_increment_( 0.5 );
	Real gamma_, gamma_min_( 0 ), gamma_max_( 180.0 ), gamma_increment_( 60.0 );

	Matrix M;
	Vector axis1( 1.0, 0.0, 0.0 ), axis2( 0.0, 1.0, 0.0 ), axis3( 0.0, 0.0, 1.0 );


	Size  count( 0 );
	for ( alpha_ = alpha_min_; alpha_ <= alpha_max_;  alpha_ += alpha_increment_ ){

		//std::cout << i++ << " out of " << N_SAMPLE_ALPHA << ". Current count: " << count_total_ <<
		//			". num poses that pass cuts: " << count_good_ << std::endl;

		for ( Real cosbeta = cosbeta_min_; cosbeta <= cosbeta_max_;  cosbeta += cosbeta_increment_ ){
			if ( cosbeta < -1.0 ){
				beta_ = -1.0 * degrees( std::acos( -2.0 - cosbeta ) );
			} else if ( cosbeta > 1.0 ){
				beta_ = -1.0 * degrees( std::acos( 2.0 - cosbeta ) );
			} else {
				beta_ = degrees( std::acos( cosbeta ) );
			}

			std::cout << "BETA: " << beta_ << std::endl;

			// Try to avoid singularity at pole.
			Real gamma_min_local = gamma_min_;
			Real gamma_max_local = gamma_max_;
			Real gamma_increment_local = gamma_increment_;
			if ( (beta_<-179.999 || beta_>179.999) ){
				gamma_min_local = 0.0;
				gamma_max_local = 0.0;
				gamma_increment_local = 1.0;
			}

			for ( gamma_ = gamma_min_local; gamma_ <= gamma_max_local;  gamma_ += gamma_increment_local ){

				protocols::swa::create_euler_rotation( M, alpha_, beta_, gamma_, axis1, axis2, axis3 );

				kinematics::Jump jump = pose.jump( num_jump );
				jump.set_rotation( M );
				pose.set_jump( num_jump, jump );

				pose.dump_pdb( "S_" + ObjexxFCL::string_of( count++ ) + ".pdb" );

			} // gamma
		} // beta
	}// alpha


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

	core::conformation::ResidueOP new_res;

	bool const sample_water_ = option[ sample_water ]();

	if ( sample_water_ ) {
		core::chemical::ResidueTypeCAPs const & rsd_type_list ( residue_set.name3_map ( "TP3" ) );
		new_res = ( core::conformation::ResidueFactory::create_residue ( *rsd_type_list[1] ) );
	} else {
		core::chemical::ResidueTypeCAPs const & rsd_type_list ( residue_set.name3_map ( "CCC" ) );
		new_res = ( core::conformation::ResidueFactory::create_residue ( *rsd_type_list[1] ) );
	}
	pose.append_residue_by_jump ( *new_res , 2 );

	pose.dump_pdb( "START.pdb" );


	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );


	//////////////////////////////////////////////////
	// Set up fold tree -- "chain break" between two ligands, right?
	// Uh, why not?
	kinematics::FoldTree f ( pose.fold_tree() );
	std::string probe_atom_name = sample_water_ ? " O  " : " C1 ";
	f.set_jump_atoms( 2,"ORIG",probe_atom_name);
	pose.fold_tree( f );

	//////////////////////////////////////////////////////////////////
	// displace in z by 2.0 A... just checking coordinate system
	//This jump should code for no translation or rotation -- two_benzenes.pdb
	// has the two benzenes perfectly superimposed.
	kinematics::Jump jump( pose.jump( 2 ) );

	jump.set_translation( Vector( 5.0, 0.0, 0.0 ) );
	pose.set_jump( 2, jump );
	pose.dump_pdb( "shift_x.pdb" );

	jump.set_translation( Vector( 0.0, 5.0, 0.0 ) );
	pose.set_jump( 2, jump );
	pose.dump_pdb( "shift_y.pdb" );

	jump.set_translation( Vector( 0.0, 0.0, 5.0 ) );
	pose.set_jump( 2, jump );
	pose.dump_pdb( "shift_z.pdb" );

	/// This is a code snippet to test if we are sampling water rotations properly -- could make this a little class,
	//  and then include within translation scan.
	if ( sample_water_ ) {
		sample_all_rotations_at_jump( pose, 2 );
	}

	//////////////////////////////////////////////////////////////////
	// OK, how about a score function?
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "rna_hires" );
	//scorefxn->set_weight( hack_elec, 1.0 );

	(*scorefxn)( pose );
	scorefxn->show( std::cout, pose );

	//////////////////////////////////////////////////////////////////
	// compute scores on a plane for now.
	using namespace core::io::silent;
	SilentFileData silent_file_data;
	std::string const silent_file = option[ out::file::silent  ]();

	Size n( 0 );

	for (Size i = -200; i <= 200; ++i) {
		for (Size j = -200; j <= 200; ++j) {
			Real const x = i * 0.1;
			Real const y = j * 0.1;
			Real const z = 0.0;
			jump.set_translation( Vector( x, y, z ) ) ;
			pose.set_jump( 2, jump );
			(*scorefxn)( pose );

			n++;
			std::string out_tag( "S_" + lead_zero_string_of(n,10) );
			BinaryProteinSilentStruct s( pose,  out_tag);
			s.add_energy( "x", x );
			s.add_energy( "y", y );
			silent_file_data.write_silent_struct( s, silent_file, true /*write score only*/ );
		}
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

	NEW_OPT( sample_water, "use a water probe instead of carbon", false );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	core::init(argc, argv);
	my_main();
}
