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
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Map.Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/util.hh>
#include <core/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/StepWiseLegacyClusterer.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/file/file_sys_util.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

#include <core/util/Tracer.hh>
using core::util::T;

// option key includes

#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/cluster.OptionKeys.gen.hh>


using core::util::Error;
using core::util::Warning;

using namespace core;
using namespace protocols;
using namespace core::options::OptionKeys;

using utility::vector1;

using ObjexxFCL::format::A;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;

using numeric::conversions::radians;
using numeric::conversions::degrees;



typedef  numeric::xyzMatrix< Real > Matrix;

OPT_KEY( Real, xyz_sample )
OPT_KEY( Real, box_radius )
OPT_KEY( Real, temperature )
OPT_KEY( Integer, n_sample )
OPT_KEY( Integer, n_sample_beta )
OPT_KEY( Real, score_cutoff )
OPT_KEY( Real, contact_cutoff )
OPT_KEY( Real, steric_dist_cutoff )
OPT_KEY( Integer, min_contacts )
OPT_KEY( Boolean, only_positive_Z )

///////////////////////////////////////////////////////////////////////
void
setup_heavy_atoms( pose::Pose const & pose,
									 utility::vector1< Vector > & pose_atoms,
									 utility::vector1< Size > const & subset_res ){

	pose_atoms.clear();

	for ( Size n = 1; n <= subset_res.size(); n++ ) {
		Size const i = subset_res[ n ];

		for ( Size j = 1; j <= pose.residue_type( i ).nheavyatoms(); j++ ){
			if ( pose.residue(i).is_virtual( j ) ) continue;
			pose_atoms.push_back( pose.xyz( core::id::AtomID(j,i) ) );
		}

	}

}


///////////////////////////////////////////////////////////////////////
bool
check_contact( Vector const & translation,
							 utility::vector1< Vector > const & moving_atoms,
							 utility::vector1< Vector > const & partner_atoms
							 ){

	using namespace core::options;
	using namespace core::options::OptionKeys;

	static Distance const DIST_CUTOFF = option[ contact_cutoff ]();
	static Real const DIST_CUTOFF_squared = DIST_CUTOFF * DIST_CUTOFF;
	static Size const MIN_CONTACTS = option[ min_contacts ]();

	Size num_contacts( 0 );

	for ( Size i = 1; i <= moving_atoms.size(); i++ ) {
		Vector const test_cbeta = moving_atoms[ i ] + translation;

		for ( Size j = 1; j <= partner_atoms.size(); j++ ) {
			Vector const & partner_cbeta = partner_atoms[ j ];

			if ( ( test_cbeta - partner_cbeta ).length_squared() < DIST_CUTOFF_squared ) {
				num_contacts ++;
				if ( num_contacts >= MIN_CONTACTS ) return true;
			}

		}
	}

	return false;

}


///////////////////////////////////////////////////////////////////////
bool
check_steric_overlap( Vector const & translation,
											utility::vector1< Vector > const & moving_atoms,
											utility::vector1< Vector > const & partner_atoms
											){

	using namespace core::options;
	using namespace core::options::OptionKeys;

	static Distance const DIST_CUTOFF = option[ steric_dist_cutoff ]();
	static Real const DIST_CUTOFF_squared = DIST_CUTOFF * DIST_CUTOFF;

	for ( Size i = 1; i <= moving_atoms.size(); i++ ) {
		Vector const test_atom = moving_atoms[ i ] + translation;

		for ( Size j = 1; j <= partner_atoms.size(); j++ ) {
			Vector const & partner_atom = partner_atoms[ j ];
			if ( ( test_atom - partner_atom ).length_squared() < DIST_CUTOFF_squared ) return false;
		}
	}

	return true;

}

///////////////////////////////////////////////////////////////
void
search_translations( pose::Pose & pose,
										 pose::Pose const & pose_to_translate,
										 utility::vector1< Size > const & moving_res,
										 utility::vector1< Size > const & partner_res,
										 Size & count,
										 Size & positive_Z_count,
										 Real & best_energy,
										 Real & Z,
										 core::io::silent::SilentFileDataOP sfd = 0 ){

	using namespace core::io::silent;
	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::scoring;
	using namespace protocols::swa;

	Real const xyz_increment = option[ xyz_sample ]();
	Real const tether_radius = option[ box_radius ]();
	Real energy_cutoff = option[ score_cutoff ]();
	Size const N_SAMPLE_TRANSLATE  = 2 * static_cast< Size >( (tether_radius/xyz_increment) + 0.5 ) + 1;
	Real const tether_radius_squared = tether_radius * tether_radius;
	Real const beta = 1.0/ option[temperature]();

	bool const positive_Z =  option[ only_positive_Z ]();

	static ScoreFunctionOP scorefxn = get_score_function();

	utility::vector1< Vector > moving_atoms, partner_atoms;
	setup_heavy_atoms( pose_to_translate, moving_atoms , moving_res );
	setup_heavy_atoms( pose_to_translate, partner_atoms, partner_res );


	for ( Size i = 1; i <= N_SAMPLE_TRANSLATE; i++ ) {
		Real const delx = ( static_cast<Real>(i) - 1 ) * xyz_increment - tether_radius;

		for ( Size j = 1; j <= N_SAMPLE_TRANSLATE; j++ ) {
			Real const dely = ( static_cast<Real>(j) - 1) * xyz_increment - tether_radius;

			for ( Size k = 1; k <= N_SAMPLE_TRANSLATE; k++ ) {
				Real const delz = ( static_cast<Real>(k) - 1) * xyz_increment - tether_radius;

				Vector shift_from_tether( delx, dely, delz );
				if ( shift_from_tether.length_squared() > tether_radius_squared ) continue;

				Vector translation = shift_from_tether;
				count++;

				std::string const tag = "S_" + ObjexxFCL::lead_zero_string_of( count, 6 );

				// These can be sped up using grid indexing.
				if ( count > 1 && !check_contact(        translation, moving_atoms, partner_atoms ) ) continue;
				if ( count > 1 && !check_steric_overlap( translation, moving_atoms, partner_atoms ) ) continue;

				translate( pose, translation, pose_to_translate, moving_res );
				Real const energy = (*scorefxn)( pose );

				// Assume first energy computed will work as a "reference" energy.
				static Real const energy_reference = energy;
				Real const Z_increment = exp( -1.0 * beta * ( energy - energy_reference ) ) - 1.0;
				//std::cout << Z_increment << std::endl;
				if ( Z_increment < 0.0 ) {
					positive_Z_count++;
					if ( positive_Z ) continue;
				}
				Z += Z_increment;

				if ( energy < (best_energy + energy_cutoff) && sfd ){
					BinarySilentStruct s( pose, tag );
					sfd->add_structure( s );
				}
				if ( energy < best_energy ) best_energy = energy;


			}
		}
	}

}


///////////////////////////////////////////////////////////////
void
define_states_test(){

	using namespace core::options;
	using namespace core::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::pose;
	using namespace core::io::silent;
	using namespace protocols::swa;

	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	// Create C-G base pair.
	Pose pose;
	make_pose_from_sequence( pose, "cg",	*rsd_set );
	FoldTree f( 2 );
	f.new_jump( 1, 2, 1);
	f.set_jump_atoms( 1,
										core::chemical::rna::chi1_torsion_atom( pose.residue( 1) ),
										core::chemical::rna::chi1_torsion_atom( pose.residue( 2) )   );
	pose.dump_pdb( "start.pdb" );

	// Virtualize backbone
	add_variant_type_to_pose_residue( pose, "VIRTUAL_BACKBONE_EXCEPT_C1PRIME", 1 );
	add_variant_type_to_pose_residue( pose, "VIRTUAL_BACKBONE_EXCEPT_C1PRIME", 2 );
	pose.dump_pdb( "virtualize.pdb" );

	// Align second residue to first. Dump it to check.
	RNA_CentroidInfo rna_centroid_info;

	Vector centroid1 = rna_centroid_info.get_base_centroid( pose.residue(1) );
	Stub s1 = rna_centroid_info.get_base_coordinate_system( pose.residue(1), centroid1 );

	Vector centroid2 = rna_centroid_info.get_base_centroid( pose.residue(2) );
	Stub s2 = rna_centroid_info.get_base_coordinate_system( pose.residue(2), centroid2 );

	utility::vector1< Size > moving_res1, moving_res2;
	moving_res1.push_back ( 1 );
	moving_res2.push_back ( 2 );

	translate( pose, -centroid1, pose, moving_res1);
	translate( pose, -centroid2, pose, moving_res2);
	pose.dump_pdb( "translated.pdb" );

	rotate( pose, s1.M.transposed(), pose, moving_res1);
	rotate( pose, s2.M.transposed(), pose, moving_res2);
	pose.dump_pdb( "rotated.pdb" );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	Vector axis1( 1.0, 0.0, 0.0);
	Vector axis2( 0.0, 1.0, 0.0);
	Vector axis3( 0.0, 0.0, 1.0);


	// Rigid body sample. Keep track of total number of states so that we can extract a Kd
	// Use input parameters to define fineness of modeler -- will look for convergence.
	// Save lowest energy states.

	/////////////////////////////////////////////////////////////
	// Sample Euler angles.
	//  This probably should become a class or something.
	//  I used this in helix sampler (stepwise_helix_test), and I
	//  think Parin has independently set this up.
	////////////////////////////////////////////////////////////
	Size const N_SAMPLE = option[ n_sample ]();
	Matrix M;
	Real alpha_increment = ( 360.0 / static_cast< Real >( N_SAMPLE ) );
	Real gamma_increment = alpha_increment;
	Real cos_beta_increment = radians( alpha_increment );
	Size N_SAMPLE_COSBETA = 2.0 / cos_beta_increment;
	if ( option[ n_sample_beta ].user() ){ //override default
		N_SAMPLE_COSBETA = option[ n_sample_beta ]();
		cos_beta_increment = 2.0 / N_SAMPLE_COSBETA;
	}
	std::cout << "N_SAMPLE_COSBETA " << N_SAMPLE_COSBETA << std::endl;

	clock_t const time_start( clock() );
	Pose pose_start = pose;

  SilentFileDataOP sfd = new SilentFileData;
	Real best_energy( 999.9 );
	Size count( 0 );
	Size positive_Z_count( 0 );
	Real Z( 0.0 );

	for ( Size i = 1; i <= N_SAMPLE; i++ ){

		Real const alpha = static_cast<Real>( i ) * alpha_increment;
		std::cout << i << " out of " << N_SAMPLE << ". Current count: " << count << ". num poses with positive Z: " << positive_Z_count << std::endl;

		for ( Size j = 1; j <= N_SAMPLE_COSBETA; j++ ){
			Real const cos_beta = -1.0 + ( j - 0.5 ) * cos_beta_increment;
			Real const beta = degrees( std::acos( cos_beta ) );

			for ( Size k = 1; k <= N_SAMPLE; k++ ){
				Real const gamma = static_cast< Real >( k ) * gamma_increment + 0.01;

				create_euler_rotation( M, alpha, beta, gamma, axis1, axis2, axis3 );

				rotate( pose, M, pose_start, moving_res2 );

				Pose pose_to_translate = pose;

				search_translations( pose,
														 pose_to_translate,
														 moving_res2,
														 moving_res1,
														 count,
														 positive_Z_count,
														 best_energy,
														 Z,
														 sfd );

			}
		}
	}


	// Cluster lowest energy states and output.
	protocols::stepwise::StepWiseLegacyClusterer stepwise_clusterer(  sfd );
	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer.set_max_decoys( max_decoys );
	stepwise_clusterer.set_cluster_by_all_atom_rmsd( true );
	stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );
	stepwise_clusterer.set_calc_rms_res( moving_res2 );
	Real cluster_radius( 0.25 );
	if ( option[ OptionKeys::cluster::radius ].user() ) cluster_radius = option[ OptionKeys::cluster::radius ]();
	stepwise_clusterer.set_cluster_radius( cluster_radius	);
	stepwise_clusterer.cluster();
	std::string const silent_file = option[ out::file::silent]();
	stepwise_clusterer.output_silent_file( silent_file );

	// Print out Kd.
	Real const xyz_increment = option[ xyz_sample ]();
	Real const concentration = (1.0 / 6.022e-4);
	std::cout << "modeler concentration " << concentration << std::endl;

	Real const additional_probability = Z * xyz_increment * xyz_increment * xyz_increment * alpha_increment / ( N_SAMPLE * N_SAMPLE_COSBETA * N_SAMPLE );
	std::cout << "additional_probability: " <<  additional_probability << std::endl;

	Real const Kd = concentration * (1.0 / additional_probability);
	std::cout << "Kd " << Kd << std::endl;

}

///////////////////////////////////////////////////////////////
void
turner_rules_test(){

	// Create C-G base pair based on state n
	// For now assume delta, chi are at ideal values.

	// Add on next C base
	// Sample its conformation [epsilon,zeta,alpha,beta,gamma]
	// Disallow steric clashes.
  // Save torsion angles, centroid, and base coordinate system.

	// Add on next G base -- same thing as above.

	// Loop over possible base pairing states of newly modeled bases.
	// Find torsion angle combinations that might qualify. Build those poses
	// and save energies.


}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace core::options;

	define_states_test();

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace core::options;


	NEW_OPT( n_sample, "number of samples per torsion angle", 18 );
	NEW_OPT( n_sample_beta, "number of samples in tilt angle beta", 18 );
	NEW_OPT( xyz_sample, "spacing in xyz search, in Angstroms", 1.0 );
	NEW_OPT( box_radius, "spacing in xyz search, in Angstroms", 10.0 );
	NEW_OPT( score_cutoff, "Scoring cutoff", 10.0 );
	NEW_OPT( temperature, "Temperature", 3.0 );
	NEW_OPT( contact_cutoff, "how close atoms need to be to define contact", 4.5 );
	NEW_OPT( steric_dist_cutoff, "how close heavy atoms need to be to define clash", 2.5 );
	NEW_OPT( min_contacts, "minimum number of contacts", 1 );
	NEW_OPT( only_positive_Z, "only allow positive contributions to partition function", false );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
