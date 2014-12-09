// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <core/scoring/ScoringManager.hh>
#include <core/chemical/rna/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <utility/io/ozstream.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/init/init.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/import_pose/import_pose.hh>


#include <protocols/stepwise/sampling/util.hh> //has euler angle stuff.

#include <protocols/viewer/viewers.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

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
using namespace numeric::conversions;
using namespace protocols::swa;

OPT_KEY( Boolean, sample_water )
OPT_KEY( Boolean, sample_another_adenosine )
OPT_KEY( Boolean, just_xy )
OPT_KEY( Boolean, quick_score )
OPT_KEY( String, copy_adenosine_adenosine_file )
OPT_KEY( Real, alpha_increment )
OPT_KEY( Real, cosbeta_increment )
OPT_KEY( Real, gamma_increment )
OPT_KEY( Real, xyz_increment )
OPT_KEY( Real, xyz_size )

//
// To sample a 'carbon' probe atom:
// adenosine_sample_around.macosgccrelease  -database ~/rosetta_database/ -s a_RNA.pdb
//
// To sample a water
// adenosine_sample_around.macosgccrelease  -database ~/rosetta_database/ -s a_RNA.pdb  -sample_water
//
// To sample an adenosine
// adenosine_sample_around.macosgccrelease  -database ~/rosetta_database/ -s a_RNA.pdb  -sample_another_adenosine   -copy_adenosine_adenosine_file double_A_ready_set.pdb
//
// To sample an adenosine, reading in a starting adenosine-adenosine pairing conformation.
// adenosine_sample_around.macosgccrelease  -database ~/rosetta_database/ -s a_RNA.pdb  -sample_another_adenosine   -copy_adenosine_adenosine_file double_A_ready_set.pdb
//
//
//  NOTE! May need to add flags:
//   -extra_res ~/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/water/TP3.params    /Users/rhiju/rosetta_database/chemical/residue_type_sets/rna/residue_types/extra/C.params
//


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
	core::chemical::ResidueTypeCOPs const & rsd_type_list ( residue_set.name3_map ( "VRT" ) );
	core::conformation::ResidueOP new_res ( core::conformation::ResidueFactory::create_residue ( *rsd_type_list[1] ) );
	pose.append_residue_by_jump ( *new_res , 1 );

	// make the virt atom the root
	if ( set_res_as_root ) {
		kinematics::FoldTree newF ( pose.fold_tree() );
		newF.reorder ( nres + 1 );
		pose.fold_tree ( newF );
	}
}

void
add_another_virtual_res ( core::pose::Pose & pose ) {
	int nres = pose.total_residue();
	// attach virt res there
	bool fullatom = pose.is_fullatom();
	core::chemical::ResidueTypeSet const & residue_set = pose.residue_type ( 1 ).residue_type_set();
	core::chemical::ResidueTypeCOPs const & rsd_type_list ( residue_set.name3_map ( "VRT" ) );
	core::conformation::ResidueOP new_res ( core::conformation::ResidueFactory::create_residue ( *rsd_type_list[1] ) );
	pose.append_residue_by_jump ( *new_res , pose.total_residue() );
}


/////////////////////////////////////////////////////////////////////////////////////////////
// Rhiju -- rotate to my favorite frame. Base centroid is now at origin.
//         X points to N1 atom. Z points normal to base. Y is orthonormal and points towards Hoogsteen edge, I think.
void
rotate_into_nucleobase_frame( core::pose::Pose & pose ){

	using namespace core::conformation;
	using namespace core::chemical::rna;
	using namespace core::id;

	// assuming pose has an RNA at residue 1 -- will rotate just that residue.
	Size const base_pos( 1 );
	Residue const & rsd = pose.residue( base_pos );

	Vector centroid = get_rna_base_centroid( rsd, true /*verbose*/ );
	Matrix M = get_rna_base_coordinate_system( rsd, centroid );
	kinematics::Stub stub( M, centroid );

	for ( Size n = 1; n <= pose.total_residue(); n++ ){
		for (Size i = 1; i <= pose.residue(n).natoms(); i++ ){
			Vector xyz_new = stub.global2local( pose.residue(n).xyz( i ) ); // it is either this or M-inverse.
			pose.set_xyz( AtomID( i, n ), xyz_new );
		}
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////
//Measure the centroid distance between two bases
Real
centroid_dist( core::pose::Pose & pose ){

	using namespace core::conformation;
	using namespace core::chemical::rna;
	using namespace core::id;


	bool const sample_another_adenosine_ = option[ sample_another_adenosine ]();
	Residue const & rsd1 = pose.residue( 1 );
	Residue const & rsd2 = sample_another_adenosine_ ? pose.residue( 4 ) : pose.residue( 3 );

	Vector centroid1, centroid2;
	if ( rsd1.is_RNA() ) {
		centroid1 = get_rna_base_centroid( rsd1, false /*verbose*/ );
	} else {
		centroid1 = rsd1.nbr_atom_xyz(); //Just use the nbr atom if not RNA
	}

	if ( rsd2.is_RNA() ) {
		centroid2 = get_rna_base_centroid( rsd2, false /*verbose*/ );
	} else {
		centroid2 = rsd2.nbr_atom_xyz(); //Just use the nbr atom if not RNA
	}

	return (centroid1 - centroid2).length();
}


///////////////////////////////////////////////////////////////////////////////////////////
// This is imported from protocols/stepwise/RigidBodySampler.cco
Real
sample_all_rotations_at_jump( pose::Pose & pose, Size const num_jump, scoring::ScoreFunctionOP scorefxn = 0 ){

	Real alpha_, alpha_min_( 0 ), alpha_max_( 180.0 ), alpha_increment_( option[ alpha_increment ]() );
	Real beta_, cosbeta_min_( -1.0 ), cosbeta_max_( 1.0 ), cosbeta_increment_( option[ cosbeta_increment ]()  );
	Real gamma_, gamma_min_( 0 ), gamma_max_( 180.0 ), gamma_increment_( option[ gamma_increment ]() );

	Matrix M;
	Vector axis1( 1.0, 0.0, 0.0 ), axis2( 0.0, 1.0, 0.0 ), axis3( 0.0, 0.0, 1.0 );

	Real const kT = 0.5;
	Real partition_function = 0;
	Size  count( 0 );
	Real  score_min( 0.0 );
	kinematics::Jump  best_jump;

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

			//std::cout << "BETA: " << beta_ << std::endl;

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

				protocols::stepwise::create_euler_rotation( M, alpha_, beta_, gamma_, axis1, axis2, axis3 );

				kinematics::Jump jump = pose.jump( num_jump );
				jump.set_rotation( M );
				pose.set_jump( num_jump, jump );

				if ( scorefxn ) {
					Real const score = (*scorefxn)( pose );
					partition_function += exp( - score / kT );
					if ( score < score_min || count == 0 ) {
						score_min = score;
						best_jump = jump;
					}
				} else {
					// this is a test
					pose.dump_pdb( "S_" + ObjexxFCL::string_of( count ) + ".pdb" );
				}

				count++;

			} // gamma
		} // beta
	}// alpha

	pose.set_jump( num_jump, best_jump );

	Real const free_E = - log( partition_function / count );

//	std::cout << "Energies: " << free_E << ' ' << score_min << std::endl;
//	return score_min;
	return free_E;

}


/////////////////////////////////////////////////////////////////////////////////
Real
do_scoring( pose::Pose & pose,
						scoring::ScoreFunctionOP scorefxn,
						bool const & sample_water_,
						Size const probe_jump_num ){

	if ( sample_water_ ){
		return sample_all_rotations_at_jump( pose, probe_jump_num, scorefxn );
	}

	return (*scorefxn)( pose );

}

//////////////////////////////////////////////////////////
void
do_xy_scan( pose::Pose & pose,
						scoring::ScoreFunctionOP scorefxn,
						std::string const outfile,
						Real const z,
						Size const probe_jump_num,
						Real const box_bins,
						Real const translation_increment,
						bool const sample_water_ ){

	kinematics::Jump jump = pose.jump( probe_jump_num );

	utility::io::ozstream out;
	out.open( outfile );

	Size count( 0 );
	Real best_score( 0.0 );
	Vector best_translation( 0.0, 0.0, 0.0 );

	for (int i = -box_bins; i <= box_bins; ++i) {
		for (int j = -box_bins; j <= box_bins; ++j) {
			Real const x = j * translation_increment;
			Real const y = i * translation_increment;
			jump.set_translation( Vector( x, y, z ) ) ;
			pose.set_jump( probe_jump_num, jump );
			Real score = do_scoring( pose, scorefxn, sample_water_, probe_jump_num );
			out << score << ' ' ;

			if (score < best_score || count++ == 0 ){
				best_translation = Vector( x, y, z );
				best_score = score;
			}

		}
		out << std::endl;
	}
	out.close();

	// return pose in lowest energy configuration that was found in scan.
	jump.set_translation( best_translation );
	pose.set_jump( probe_jump_num, jump );
	do_scoring( pose, scorefxn, sample_water_, probe_jump_num );
}

/////////////////////////////////////////////////////////////////////////////////
void
adenine_probe_score_test()
{
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::id;

	//////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	// Read in pose with two methane. "Z" = ligand. Note need flag:
	//         -extra_res_fa CH4.params -s two_methane.pdb
	pose::Pose pose;
	std::string infile  = option[ in ::file::s ][1];
	import_pose::pose_from_pdb( pose, *rsd_set, infile );

	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );
	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RIBOSE", 1 );

	rotate_into_nucleobase_frame( pose );
	pose.dump_pdb( "a_rotated.pdb" );

	add_virtual_res(pose);
	core::chemical::ResidueTypeSet const & residue_set = pose.residue_type ( 1 ).residue_type_set();

	core::conformation::ResidueOP new_res;

	bool const sample_water_ = option[ sample_water ]();
	bool const sample_another_adenosine_ = option[ sample_another_adenosine ]();

	if ( sample_another_adenosine_ ) add_another_virtual_res(pose); // this is the coordinate system for the next base.

	if ( sample_water_ ) {
		core::chemical::ResidueTypeCOPs const & rsd_type_list ( residue_set.name3_map ( "TP3" ) );
		new_res = ( core::conformation::ResidueFactory::create_residue ( *rsd_type_list[1] ) );
	} else if ( sample_another_adenosine_ ){
		new_res = pose.residue(1).clone();
	} else {
		core::chemical::ResidueTypeCOPs const & rsd_type_list ( residue_set.name3_map ( " DC" ) ); // just a carbon atom.
		new_res = ( core::conformation::ResidueFactory::create_residue ( *rsd_type_list[1] ) );
	}

	// doing some checks for hbonds & geom_sol are being calculated -- or not calculated.
	if ( sample_water_ || sample_another_adenosine_ ){

		for (Size n = 1; n <= new_res->natoms(); n++ ){

			std::cout << "PROBE ATOM:  " << n << ' ' << new_res->atom_name( n ); // << ' ' << new_res->atom_base( n ) << ' ' << new_res->abase2( n ) << std::endl;
			if ( new_res->atom_base( n ) > 0 ) std::cout << "   base: " << new_res->atom_name( new_res->atom_base( n ) );
			if ( new_res->abase2(n) > 0 ) std::cout << "   base2: " <<  new_res->atom_name( new_res->abase2( n ) );
			std::cout << std::endl;

		}
	}

	pose.append_residue_by_jump ( *new_res ,  pose.total_residue() );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 800, 800 );


	//////////////////////////////////////////////////
	// Set up fold tree -- "chain break" between two ligands, right?
	kinematics::FoldTree f ( pose.fold_tree() );
	std::string probe_atom_name = " C1 ";
	if (sample_water_) probe_atom_name = " O  ";
	if (sample_another_adenosine_) probe_atom_name = "ORIG"; // tricky -- the jump is from one coordinate system to the next one!
	f.set_jump_atoms( 2,"ORIG",probe_atom_name);
	pose.fold_tree( f );

	std::cout << pose.annotated_sequence() << std::endl;

	if ( option[ copy_adenosine_adenosine_file ].user() ){
		pose::Pose pose_reference;
		import_pose::pose_from_pdb( pose_reference, *rsd_set, option[copy_adenosine_adenosine_file]() );
		rotate_into_nucleobase_frame( pose_reference );

		// copy over coordinates.
		Residue const & rsd_ref = pose_reference.residue( 2 );
		Size pos2( 4 ); // the sequence of the working pose is... adenine-virtual-virtual-adenine

		for( Size i_ref = 1; i_ref <= rsd_ref.natoms(); i_ref++ ){
			Size i = pose.residue( pos2 ).atom_index(   rsd_ref.atom_name( i_ref ) );
			pose.set_xyz( AtomID( i, pos2 ),  rsd_ref.xyz( i_ref ) );
		}
	}
	pose.dump_pdb( "START.pdb" );


	//////////////////////////////////////////////////////////////////
	// displace in z by 2.0 A... just checking coordinate system
	//This jump should code for no translation or rotation -- two_benzenes.pdb
	// has the two benzenes perfectly superimposed.
	Size const probe_jump_num( 2 );
	kinematics::Jump jump( pose.jump( probe_jump_num ) );

	jump.set_translation( Vector( 5.0, 0.0, 0.0 ) );
	pose.set_jump( probe_jump_num, jump );
	pose.dump_pdb( "shift_x.pdb" );

	jump.set_translation( Vector( 0.0, 5.0, 0.0 ) );
	pose.set_jump( probe_jump_num, jump );
	pose.dump_pdb( "shift_y.pdb" );

	jump.set_translation( Vector( 0.0, 0.0, 5.0 ) );
	pose.set_jump( probe_jump_num, jump );
	pose.dump_pdb( "shift_z.pdb" );

	/// This is a code snippet to test if we are sampling water rotations properly -- could make this a little class,
	//  and then include within translation scan.
	//	if ( sample_water_ )		sample_all_rotations_at_jump( pose, probe_jump_num );

	//////////////////////////////////////////////////////////////////
	// OK, how about a score function?
	ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ){
		scorefxn = scoring::get_score_function();
	} else {
		scorefxn = ScoreFunctionFactory::create_score_function( "farna/rna_hires" );
		scorefxn->set_weight( rna_sugar_close, 0.0 ); //still computed with virtual sugar? weird.
	}

//	jump.set_translation( Vector( 3.75, 1.75, 1.5 ) );
//	pose.set_jump( probe_jump_num, jump );
//	sample_all_rotations_at_jump( pose, probe_jump_num, scorefxn );
//	pose.dump_pdb( "test.pdb" );
//	exit(0);
 	(*scorefxn)( pose );
	scorefxn->show( std::cout, pose );

	//////////////////////////////////////////////////////////////////
	// compute scores on a plane for now.

	Real const box_size = option[ xyz_size ]();
	Real const translation_increment = option[ xyz_increment ]();
	int box_bins = int( box_size/translation_increment );

	using namespace core::io::silent;
	SilentFileData silent_file_data;
	utility::io::ozstream out;

	//////////////////////////////////////////////
	std::cout << "Doing Z scan..." << std::endl;
	out.open( "score_z.table" );
	for (int i = -box_bins; i <= box_bins; ++i) {
		Real const x = 0.0;
		Real const y = 0.0;
		Real const z = i * translation_increment;
		jump.set_translation( Vector( x, y, z ) ) ;
		pose.set_jump( probe_jump_num, jump );
		out << z << ' ' << centroid_dist( pose )  << ' ' << do_scoring( pose, scorefxn, sample_water_, probe_jump_num ) << std::endl;
	}
	out.close();


	//////////////////////////////////////////////
	std::cout << "Doing XY scan... Z =  0.0" << std::endl;
	do_xy_scan( pose, scorefxn, "score_xy_0.table", 0.0, probe_jump_num, box_bins, translation_increment, sample_water_ );

	pose.dump_pdb( "best_xy.pdb");
	if ( option[ just_xy ]() ) return;


	//////////////////////////////////////////////
	std::cout << "Doing XY scan... Z = +1.5" << std::endl;
	do_xy_scan( pose, scorefxn, "score_xy_1.5.table", 1.5, probe_jump_num, box_bins, translation_increment, sample_water_ );

	std::cout << "Doing XY scan... Z = +3.0" << std::endl;
	do_xy_scan( pose, scorefxn, "score_xy_3.table", 3.0, probe_jump_num, box_bins, translation_increment, sample_water_ );

	std::cout << "Doing XY scan... Z = +4.0" << std::endl;
	do_xy_scan( pose, scorefxn, "score_xy_4.table", 4.0, probe_jump_num, box_bins, translation_increment, sample_water_ );
	// Following are exactly the same as +1.0 and +3.0 when modeling nucleobase.
	//std::cout << "Doing XY scan... Z = -1.0" << std::endl;
	//	do_xy_scan( pose, scorefxn, "score_para_0_table", 1.0, probe_jump_num, box_bins, translation_increment, sample_water_ );

	//	std::cout << "Doing XY scan... Z = -3.0" << std::endl;
	//	do_xy_scan( pose, scorefxn, "score_para_0_table", 3.0, probe_jump_num, box_bins, translation_increment, sample_water_ );

	//////////////////////////////////////////////
	std::cout << "Doing XZ scan..." << std::endl;
	out.open( "score_xz.table" );
	for (int i = -box_bins; i <= box_bins; ++i) {
		for (int j = -box_bins; j <= box_bins; ++j) {
			Real const x = j * translation_increment;
			Real const z = i * translation_increment;
			Real const y = 0.0;
			jump.set_translation( Vector( x, y, z ) ) ;
			pose.set_jump( probe_jump_num, jump );
			out << do_scoring( pose, scorefxn, sample_water_, probe_jump_num ) << ' ' ;
		}
		out << std::endl;
	}
	out.close();

	//////////////////////////////////////////////
	std::cout << "Doing YZ scan..." << std::endl;
	out.open( "score_yz.table" );
	for (int i = -box_bins; i <= box_bins; ++i) {
		for (int j = -box_bins; j <= box_bins; ++j) {
			Real const y = j * translation_increment;
			Real const z = i * translation_increment;
			Real const x = 0.0;
			jump.set_translation( Vector( x, y, z ) ) ;
			pose.set_jump( probe_jump_num, jump );
			out << do_scoring( pose, scorefxn, sample_water_, probe_jump_num ) << ' ' ;
		}
		out << std::endl;
	}
	out.close();


}

///////////////////////////////////////////////////////////////
void
quick_score_test(){

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::id;

	//////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	// Read in pose with two methane. "Z" = ligand. Note need flag:
	//         -extra_res_fa CH4.params -s two_methane.pdb
	pose::Pose pose;
	std::string infile  = option[ in ::file::s ][1];
	import_pose::pose_from_pdb( pose, *rsd_set, infile );

	for ( Size i = 1; i <= pose.total_residue(); i++ ){
		if (pose.residue(i).is_RNA() ){
			pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", i );
			pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RIBOSE", i );
		}
	}

	ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ){
		scorefxn = scoring::get_score_function();
	} else {
		scorefxn = ScoreFunctionFactory::create_score_function( "farna/rna_hires" );
		scorefxn->set_weight( rna_sugar_close, 0.0 ); //still computed with virtual sugar? weird.
	}

	(*scorefxn)( pose );
	scorefxn->show( std::cout, pose );

}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	if ( option[quick_score]() ){
		quick_score_test();
	} else {
		adenine_probe_score_test();
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
	NEW_OPT( just_xy, "Just scan x, y at z=0", false );
	NEW_OPT( quick_score, "alternative mode for checking geom_sol, etc.", false );
	NEW_OPT( copy_adenosine_adenosine_file, "get rigid body relation between two adenosines from file", "" );
	NEW_OPT( alpha_increment, "input parameter", 40.0 );
	NEW_OPT( cosbeta_increment, "input parameter", 0.25 );
	NEW_OPT( gamma_increment, "input parameter", 40.0 );
	NEW_OPT( xyz_increment, "input parameter", 0.2 );
	NEW_OPT( xyz_size, "input parameter", 10.0 );

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
