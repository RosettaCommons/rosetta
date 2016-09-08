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
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <protocols/viewer/viewers.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/options/option.hh>
#include <core/options/after_opts.hh>
#include <core/options/util.hh>
#include <core/options/option_macros.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/util/basic.hh>
#include <core/io/database/open.hh>
#include <devel/init.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/stepwise/modeler/util.hh>

// C++ headers
//#include <cstdlib>
#include <time.h>

//silly using/typedef

#include <core/util/Tracer.hh>
using core::util::T;

// option key includes

#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/swa.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>


using core::util::Error;
using core::util::Warning;

using namespace core;
using namespace protocols;
using namespace core::options;
using namespace core::options::OptionKeys;

using numeric::conversions::radians;
using numeric::conversions::degrees;

using utility::vector1;



typedef  numeric::xyzMatrix< Real > Matrix;
//typedef std::map< std::string, core::pose::PoseOP > PoseList;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, spinner )
OPT_KEY( Boolean, prepend )
OPT_KEY( Boolean, append )
OPT_KEY( Boolean, prepack )
OPT_KEY( Integer, n_sample )
OPT_KEY( Integer, n_sample_beta )
OPT_KEY( Integer, num_loop_res )
OPT_KEY( Integer, min_hydrophobic_contacts )
OPT_KEY( IntegerVector, sample_res )
OPT_KEY( Real, xyz_sample )
OPT_KEY( Real, hydrophobic_cbeta_dist_cutoff )
OPT_KEY( Real, steric_dist_cutoff )
OPT_KEY( Real, filter_rms )
OPT_KEY( String, base )


///////////////////////////////////////////////////////////////////////
void
get_moment_of_inertia( pose::Pose const & pose,
											 Vector & center_of_mass,
											 Matrix & moment_of_inertia,
											 utility::vector1< Size > const & moving_res ){

	center_of_mass = Vector( 0.0 );
	Size natoms( 0 );

	for ( Size n = 1; n <= moving_res.size(); n++ ) {
		Size const i = moving_res[ n ];

		for ( Size j = 1; j < pose.residue_type( i ).first_sidechain_atom(); j++ ) {
			center_of_mass += pose.residue( i ).xyz( j );
			natoms++;
		}
	}
	center_of_mass /= natoms;

	// Generate inertia tensor:
	for ( int k = 1; k <= 3; ++k ) {
		for ( int j = 1; j <= 3; ++j ) {
			Real cross_term = 0.0;

			for ( Size n = 1; n <= moving_res.size(); n++ ) {
				Size const i = moving_res[ n ];

				for ( Size m = 1; m < pose.residue_type( i ).first_sidechain_atom(); m++ ) {
					Vector pos_recenter = pose.residue( i ).xyz( m ) - center_of_mass;
					cross_term -= pos_recenter(j) * pos_recenter(k);
				}
			}

			moment_of_inertia(k,j) = cross_term;
		}
	}

	Real diagonal_term( 0.0 );

	for ( Size n = 1; n <= moving_res.size(); n++ ) {
		Size const i = moving_res[ n ];

		for ( Size m = 1; m < pose.residue_type( i ).first_sidechain_atom(); m++ ) {
			Vector pos_recenter = pose.residue( i ).xyz( m ) - center_of_mass;
			for ( Size j = 1; j <= 3; j++ ){
				diagonal_term += pos_recenter( j ) * pos_recenter( j );
			}
		}
	}

	for ( int j = 1; j <= 3; ++j ) {
		moment_of_inertia(j,j) += diagonal_term;
	}

	moment_of_inertia /= natoms;

	//	std::cout << M << std::endl;
	//	std::cout << center_of_mass << std::endl;
}


///////////////////////////////////////////////////////////////////////
void
get_euler_axes( Matrix const moment_of_inertia,
								Vector & axis1,
								Vector & axis2,
								Vector & axis3,
								Real & m1, Real & m2, Real & m3,
								Matrix & M
								)
{

	Vector xyz_w_w;
	Matrix xyz_eVec;

	Real const tol( 0.000001) ;
	xyz_w_w = numeric::eigenvector_jacobi( moment_of_inertia, tol, xyz_eVec );

	xyz_eVec.transpose();

	//	std::cout << "EIGENVALUES: " << xyz_w_w << std::endl;
	//	std::cout << "EIGENVECTORS: " << xyz_eVec << std::endl;

	utility::vector1< Vector > axes;
	axes.push_back( xyz_eVec.row_x() );
	axes.push_back( xyz_eVec.row_y() );
	axes.push_back( xyz_eVec.row_z() );

	// Sort by order of eigenvalue. Silly lists.
	std::list < std::pair< Real, Size > > moments;
	for ( Size k = 1; k <= 3; k++ )	moments.push_back(   std::make_pair( xyz_w_w(k), k ) );
	moments.sort();
	moments.reverse();

	std::pair< Real, Size > eval_pair;

	eval_pair = moments.front();
	moments.pop_front();
	m1 =  std::sqrt( eval_pair.first );
	axis1 = axes[ eval_pair.second ];

	eval_pair = moments.front();
	moments.pop_front();
	m2 =  std::sqrt( eval_pair.first );
	axis2 = axes[ eval_pair.second ];

	eval_pair = moments.front();
	moments.pop_front();
	m3 =  std::sqrt( eval_pair.first );
	axis3 = axes[ eval_pair.second ];

	// This is silly, but I can't get the matrix constructor from column vectors to work properly
	for ( Size k = 1; k <= 3; k++ ){
		M(k,1) = axis1(k);
		M(k,2) = axis2(k);
		M(k,3) = axis3(k);
	}

	std::cout << "MOMENT MAGNITUDES: " << m1 << ' ' << m2 << ' ' << m3 << std::endl;

}

///////////////////////////////////////////////////////////////////////
void
get_euler_axes( Matrix const moment_of_inertia,
								Vector & axis1,
								Vector & axis2,
								Vector & axis3,
								Real & m1, Real & m2, Real & m3 ){

	Matrix M;
	get_euler_axes( moment_of_inertia, axis1, axis2, axis3, m1, m2, m3, M );

}


///////////////////////////////////////////////////////////////////////
// Real
// get_backbone_rmsd_no_super( pose::Pose const & pose, pose::Pose const & ref_pose ){

// 	using namespace core::id;

// 	Real dev2( 0.0 );
// 	Size natoms( 0 );

// 	for ( Size i = 1; i <= pose.size(); i++ ) {
// 		for ( Size j = 1; j <= pose.residue_type( i ).first_sidechain_atom(); j++ ) {
// 			dev2 += ( pose.xyz( AtomID(j,i) ) - ref_pose.xyz( AtomID(j,i) ) ).length_squared();
// 			natoms++;
// 		}
// 	}
// 	dev2 /= natoms;
// 	return std::sqrt( dev2 );

// }

// ///////////////////////////////////////////////////////////////////////
// void
// check_rmsd_to_native_pose( pose::Pose const & pose_translate, pose::PoseCOP native_pose, Real & rmsd_min, pose::Pose & best_pose ) {

// 	if ( native_pose ){
// 		Real const rmsd = get_backbone_rmsd_no_super( pose_translate, *native_pose );
// 		if ( rmsd < rmsd_min ){
// 			rmsd_min = rmsd;
// 			best_pose = pose_translate;
// 		}
// 	}

// }

///////////////////////////////////////////////////////////////////////
// This is currently AD HOC. Need to survey helical bundles
// carefully -- what is the minimal set of hydrophobic residues?
//
// Later could also return a DISTANCE -- CB to most distance atom in sidechain.
//  That would be nice -- trps and mets could "reach farther" to make contacts.
//
bool
check_hydrophobic( chemical::AA aa ){
	using namespace chemical;
	if ( aa == aa_ala ) return true;
	//	if ( aa == aa_cys ) return true;
	if ( aa == aa_phe ) return true;
	if ( aa == aa_ile ) return true;
	if ( aa == aa_leu ) return true;
	if ( aa == aa_met ) return true;
	if ( aa == aa_val ) return true;
	if ( aa == aa_pro ) return true;

	return false;
}

///////////////////////////////////////////////////////////////////////
void
setup_hydrophobic_cbetas( pose::Pose const & pose,
													utility::vector1< Vector > & pose_hydrophobic_cbetas,
													utility::vector1< Size > const & subset_res ){

	pose_hydrophobic_cbetas.clear();

	for ( Size n = 1; n <= subset_res.size(); n++ ) {
		Size const i = subset_res[ n ];

		if ( pose.residue_type( i ).has_variant_type( core::chemical::VIRTUAL_RESIDUE_VARIANT ) ) continue;
		if ( check_hydrophobic( pose.aa( i ) ) ) pose_hydrophobic_cbetas.push_back( pose.xyz( id::NamedAtomID( " CB ", i )  ) );

	}

}

///////////////////////////////////////////////////////////////////////
void
setup_backbone_atoms( pose::Pose const & pose, 	utility::vector1< Vector > & backbone_atoms,
										utility::vector1< Size > const & subset_res ){


	backbone_atoms.clear();

	for ( Size n = 1; n <= subset_res.size(); n++ ) {
		Size const i = subset_res[ n ];

		if ( pose.residue_type( i ).has_variant_type( core::chemical::VIRTUAL_RESIDUE_VARIANT ) ) continue;

		for ( Size j = 1; j <= pose.residue_type( i ).first_sidechain_atom(); j++ ){

			if ( pose.residue( i ).is_virtual( j ) ) continue;
			backbone_atoms.push_back( pose.xyz( id::AtomID( j, i ) ) );

		}
	}

}


///////////////////////////////////////////////////////////////////////
bool
check_hydrophobic_contact( Vector const & translation,
													 utility::vector1< Vector > const & moving_pose_cbetas,
													 utility::vector1< Vector > const & partner_pose_cbetas
													 ){

	static Distance const DIST_CUTOFF = option[ hydrophobic_cbeta_dist_cutoff ]();
	static Real const DIST_CUTOFF_squared = DIST_CUTOFF * DIST_CUTOFF;
	static Size const MIN_CONTACTS = option[ min_hydrophobic_contacts ]();

	Size num_contacts( 0 );

	for ( Size i = 1; i <= moving_pose_cbetas.size(); i++ ) {
		Vector const test_cbeta = moving_pose_cbetas[ i ] + translation;

		for ( Size j = 1; j <= partner_pose_cbetas.size(); j++ ) {
			Vector const & partner_cbeta = partner_pose_cbetas[ j ];

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
											utility::vector1< Vector > const & moving_pose_backbone_atoms,
											utility::vector1< Vector > const & partner_pose_backbone_atoms
											){

	static Distance const DIST_CUTOFF = option[ steric_dist_cutoff ]();
	static Real const DIST_CUTOFF_squared = DIST_CUTOFF * DIST_CUTOFF;

	for ( Size i = 1; i <= moving_pose_backbone_atoms.size(); i++ ) {
		Vector const test_atom = moving_pose_backbone_atoms[ i ] + translation;

		for ( Size j = 1; j <= partner_pose_backbone_atoms.size(); j++ ) {
			Vector const & partner_atom = partner_pose_backbone_atoms[ j ];
			if ( ( test_atom - partner_atom ).length_squared() < DIST_CUTOFF_squared ) return false;
		}
	}

	return true;

}


///////////////////////////////////////////////////////////////////////
Real
check_rmsd( Vector const & translation,
						utility::vector1< Vector > const & moving_pose_backbone_atoms,
						utility::vector1< Vector > const & native_pose_backbone_atoms
						){

	Real dev2( 0.0 );
	Size natoms( 0 );

	for ( Size i = 1; i <= moving_pose_backbone_atoms.size(); i++ ) {
		Vector const test_atom = moving_pose_backbone_atoms[ i ] + translation;
		Vector const & native_atom = native_pose_backbone_atoms[ i ];
		dev2 += ( test_atom - native_atom ).length_squared();
		natoms++;
	}

	dev2 /= natoms;
	return std::sqrt( dev2 );

}


///////////////////////////////////////////////////////////////////////
bool
check_filter_rmsd( Vector const & translation,
									 utility::vector1< Vector > const & moving_pose_backbone_atoms,
									 utility::vector1< Vector > const & native_pose_backbone_atoms,
									 Real const & filter_rmsd
									 ){

	Real dev2( 0.0 );
	Real const max_dev2 = filter_rmsd * filter_rmsd * moving_pose_backbone_atoms.size();

	Size natoms( 0 );
	for ( Size i = 1; i <= moving_pose_backbone_atoms.size(); i++ ) {
		Vector const test_atom = moving_pose_backbone_atoms[ i ] + translation;
		Vector const & native_atom = native_pose_backbone_atoms[ i ];
		dev2 += ( test_atom - native_atom ).length_squared();
		//		std::cout << natoms++ << " " << dev2 << " " << test_atom( 1) << " " << native_atom( 1) << " " << (test_atom -native_atom).length_squared() << std::endl;
		if ( dev2 > max_dev2 ) return false;
	}

	return true;

}


///////////////////////////////////////////////////////////////////////
bool
search_translations( pose::Pose & pose,
										 pose::Pose const & pose_to_translate,
										 utility::vector1< Size > const & moving_res,
										 utility::vector1< Size > const & partner_res,
										 id::AtomID const & tether_atom_id,
										 Vector const & tether_xyz,
										 Real const & tether_radius,
										 Real const filter_rmsd,
										 Size & count,
										 bool const do_steric_check,
										 std::string const & silent_file,
										 pose::PoseCOP native_pose,
										 Real & rmsd_min,
										 pose::Pose & best_pose ){

	using namespace core::io::silent;
	using namespace protocols::swa;

	Vector translation_to_tether = tether_xyz - pose_to_translate.xyz( tether_atom_id );

	// short circuit for debugging
	//translate( pose, translation_to_tether, pose_to_translate );
	//	check_rmsd_to_native_pose( pose, native_pose, rmsd_min, best_pose );
	//	return;

	Real const xyz_increment = option[ xyz_sample ]();
	Size const N_SAMPLE_TRANSLATE  = 2 * static_cast< Size >( (tether_radius/xyz_increment) + 0.5 ) + 1;
	Real const tether_radius_squared = tether_radius * tether_radius;
	static Size local_count( 0 );

	////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Vector > moving_pose_hydrophobic_cbetas, partner_pose_hydrophobic_cbetas;
	setup_hydrophobic_cbetas( pose_to_translate, moving_pose_hydrophobic_cbetas , moving_res );
	setup_hydrophobic_cbetas( pose_to_translate, partner_pose_hydrophobic_cbetas, partner_res );

	//	std::cout << "NUM CBETAs " << partner_pose_hydrophobic_cbetas.size() << " " << moving_pose_hydrophobic_cbetas.size() << std::endl;
	static SilentFileData sfd;

	////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< Vector > moving_pose_backbone_atoms, partner_pose_backbone_atoms, native_pose_backbone_atoms;
	setup_backbone_atoms( pose_to_translate, moving_pose_backbone_atoms,  moving_res  );
	setup_backbone_atoms( pose_to_translate, partner_pose_backbone_atoms, partner_res );
	setup_backbone_atoms( *native_pose,      native_pose_backbone_atoms,  moving_res  );
	//	std::cout << "NUM STERIC ATOMS " << partner_pose_backbone_atoms.size() << " " << moving_pose_backbone_atoms.size() << " " << native_pose_backbone_atoms.size() << std::endl;

	bool found_better_orientation( false );

	for ( Size i = 1; i <= N_SAMPLE_TRANSLATE; i++ ) {
		Real const delx = ( static_cast<Real>(i) - 1 ) * xyz_increment - tether_radius;

		for ( Size j = 1; j <= N_SAMPLE_TRANSLATE; j++ ) {
			Real const dely = ( static_cast<Real>(j) - 1) * xyz_increment - tether_radius;

			for ( Size k = 1; k <= N_SAMPLE_TRANSLATE; k++ ) {
				Real const delz = ( static_cast<Real>(k) - 1) * xyz_increment - tether_radius;

				Vector shift_from_tether( delx, dely, delz );
				if ( shift_from_tether.length_squared() > tether_radius_squared ) continue;

				Vector translation = translation_to_tether + shift_from_tether;

				//				std::cout << "TRNASLATION " << xyz_increment() << " " << translation(1) << " " << translation_to_tether(1) << " " << shift_from_tether(1) << " " << moving_pose_backbone_atoms[1](1) << std::endl;

				if ( filter_rmsd > 0.0 && !check_filter_rmsd( translation, moving_pose_backbone_atoms, native_pose_backbone_atoms, filter_rmsd ) ) continue;
				//std::cout << "Made it past filter_rmsd" << std::endl;

				if ( !check_hydrophobic_contact( translation, moving_pose_hydrophobic_cbetas, partner_pose_hydrophobic_cbetas ) ) continue;

				// This could be sped up to O( N ) with a grid-index.
				if ( do_steric_check && !check_steric_overlap( translation, moving_pose_backbone_atoms, partner_pose_backbone_atoms ) )	continue;

				count++;

				std::string const tag = "S_" + ObjexxFCL::lead_zero_string_of( count, 6 );

				if ( silent_file.size() > 0 ){
					Real const rmsd = check_rmsd( translation, moving_pose_backbone_atoms, native_pose_backbone_atoms );
					translate( pose, translation, pose_to_translate, moving_res );
					BinarySilentStruct s( pose, tag );
					s.add_energy( "rms", rmsd );
					sfd.write_silent_struct( s, silent_file, false /*write score only*/ );
				}

				// short cut for figuring out parameters -- no pose.
				if( native_pose ){
					Real const rmsd = check_rmsd( translation, moving_pose_backbone_atoms, native_pose_backbone_atoms );

					if ( rmsd < rmsd_min ) {
						//						std::cout << rmsd << std::endl;
						translate( pose, translation, pose_to_translate, moving_res );
						//						pose.dump_pdb( tag + ".pdb" );
						best_pose = pose;
						rmsd_min = rmsd;
						found_better_orientation = true;
					}
				}


			}
		}
	}

	return found_better_orientation;

}

/////////////////////////////////////////////////////////
void
pack_it( pose::Pose & pose ){

	using namespace core::scoring;
	using namespace core::pack;
	using namespace core::pack::task;

	ScoreFunctionOP scorefxn_ = get_score_function();

	PackerTaskOP pack_task_ = pack::task::TaskFactory::create_packer_task( pose );

	pack_task_->restrict_to_repacking();
	for (Size i = 1; i <= pose.size(); i++) {
		if ( !pose.residue(i).is_protein() ) continue;
		pack_task_->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
		pack_task_->nonconst_residue_task(i).or_ex1( true );
		pack_task_->nonconst_residue_task(i).or_ex2( true );
		pack_task_->nonconst_residue_task(i).or_include_current( true );
		if ( pose.residue(i).has_variant_type( core::chemical::VIRTUAL_RESIDUE_VARIANT ) ) {
			pack_task_->nonconst_residue_task(i).prevent_repacking();
		}
	}

	pack::pack_rotamers(  pose, *scorefxn_, pack_task_ );

}


/////////////////////////////////////////////////////////////////////////////////////////////////////
void
figure_out_best_alpha_beta_gamma(
          pose::Pose const & pose1,
					pose::Pose const & pose2,
					utility::vector1< Size > const & moving_res )
{

  Vector center_of_mass1, center_of_mass2;
	Matrix moment_of_inertia1, moment_of_inertia2;

	Vector xaxis1,yaxis1,zaxis1;
	Vector xaxis2,yaxis2,zaxis2;
	Real m1, m2, m3;
	Matrix M1, M2;

	// starting pose.
	get_moment_of_inertia( pose1, center_of_mass1, moment_of_inertia1, moving_res );
	get_euler_axes( moment_of_inertia1, xaxis1, yaxis1, zaxis1, m1, m2, m3, M1 );

	get_moment_of_inertia( pose2, center_of_mass2, moment_of_inertia2, moving_res );
	get_euler_axes( moment_of_inertia2, xaxis2, yaxis2, zaxis2, m1, m2, m3, M2 );

	Real alpha, beta, gamma;
	protocols::stepwise::get_euler_angles( alpha, beta, gamma, M1, M2 );

	std::cout << "ALPHA:    " << alpha << std::endl;
	std::cout << "BETA:     " << beta << std::endl;
	std::cout << "COS_BETA: " << cos( radians(beta) ) << std::endl;
	std::cout << "GAMMA:    " << gamma << std::endl;


}

///////////////////////////////////////////////////////////////////////
void
spinner_test(){

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::kinematics;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::id;
	using namespace core::pose;
	using namespace protocols::swa;

	// Read in segment.
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	Pose pose;
	std::string infile  = option[ in::file::s ]()[1];
	io::pdb::pose_from_file( pose, *rsd_set, infile , core::import_pose::PDB_file);

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	utility::vector1< Size > const & moving_res = option[ sample_res ]();

	/////////////////////////////////////////////////////////////////
	PoseOP native_pose;
	if (option[ in::file::native ].user() ) {
		native_pose = PoseOP( new Pose );
		std::string native_pdb_file  = option[ in::file::native ];
		io::pdb::pose_from_file( *native_pose, *rsd_set, native_pdb_file , core::import_pose::PDB_file);
		native_pose->dump_pdb( "native.pdb" );

		//convert to centroid?
		Pose native_pose_centroid = *native_pose;
		core::chemical::switch_to_residue_type_set( native_pose_centroid, core::chemical::CENTROID );
		native_pose_centroid.dump_pdb( "native_centroid.pdb" );
	}


	/////////////////////////////////////////////////////////////////
	// Find center of mass, move to origin.

  Vector center_of_mass( 0.0 );
	Matrix moment_of_inertia( 0.0 );
	get_moment_of_inertia( pose, center_of_mass, moment_of_inertia, moving_res );
	translate( pose, -1.0 * center_of_mass, pose, moving_res );
	pose.dump_pdb( "recenter.pdb" );


	// sanity check
	//get_moment_of_inertia( pose, center_of_mass, moment_of_inertia );
	//	std::cout << center_of_mass( 1 ) << std::endl;

	// Find moments of inertia, x, y, z vectors.
	// How much to rotate?
	Real m1, m2, m3;
	Vector axis1,axis2,axis3;
	get_euler_axes( moment_of_inertia, axis1, axis2, axis3, m1, m2, m3 );

	///////////////////////////
	// Sample Euler angles.
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

	Real best_alpha( 0.0 ), best_gamma( 0.0 ), best_beta( 0.0 );

	Size count( 0 );
	Real rmsd_min( 999.9 );

	AtomID tether_atom_id( 0, 0);
	Vector tether_xyz( 0.0 );
	Real tether_radius( 0.0 );
	bool do_steric_check( true );

	figure_out_best_alpha_beta_gamma( pose, *native_pose, moving_res );

	Size const nloop = option[ num_loop_res ]();
	bool const prepend_ = option[ prepend ]();
	bool const append_ = option[ append ]();
	if ( prepend_ && append_ ) utility_exit_with_message( "cannot prepend and append" );
	if ( (prepend_ || append_ ) && !option[ num_loop_res ].user() ) utility_exit_with_message( "Must specify num_loop_res" );

	Size const moving_res_start = moving_res[ 1 ];
	Size const moving_res_end = moving_res[ moving_res.size() ];

	kinematics::FoldTree f( pose.size() );
	utility::vector1< Size > partner_res;

	if ( prepend_ ){

		if ( moving_res_start != 1 ) utility_exit_with_message( "If prepend, start moving_res must be 1" );

		tether_atom_id = AtomID( NamedAtomID( " C  ", moving_res_end ), pose );
		tether_xyz = pose.xyz( NamedAtomID( " N  ", moving_res_end + 1 )  );
		tether_radius = 3.2 * nloop + 1.0;
		do_steric_check = true;

		f.new_jump( moving_res_end, moving_res_end+1, moving_res_end );
		f.reorder( pose.size() );
		pose.fold_tree( f );

		for ( Size n = moving_res_end+1; n <= pose.size(); n++ ) partner_res.push_back( n );

	} else if ( append_ ) {

		if ( moving_res_end != pose.size() ) utility_exit_with_message( "If prepend, end moving_res must be end of pose" );

		tether_atom_id = AtomID( NamedAtomID( " N  ", moving_res_start ), pose );
		tether_xyz = pose.xyz( NamedAtomID( " C  ", moving_res_start - 1 ) );
		tether_radius = 3.2 * nloop + 1.0;
		do_steric_check = true;

		f.new_jump( moving_res_start-1, moving_res_start, moving_res_start );
		pose.fold_tree( f );

		for ( Size n = 1; n <= moving_res_start-1; n++ ) partner_res.push_back( n );

	} else { /*superimpose*/

		tether_atom_id = AtomID( NamedAtomID( " C  ", pose.size() ), pose );
		tether_xyz = native_pose->xyz( tether_atom_id );
		tether_radius = 4.0;
		do_steric_check = false;

	}

	//prepack? convenient for testing later steps.
	if ( option[ prepack ]() ){
		translate( pose, Vector( 1000.0, 1000.0, 1000.0),    pose, moving_res );
		pack_it( pose );
		translate( pose, Vector( -1000.0, -1000.0, -1000.0), pose, moving_res );
	}

	Pose pose_start = pose;
	Pose best_pose = pose;

	std::cout << " TETHER_RADIUS  : " << tether_radius << std::endl;
	std::cout << " DO STERIC CHECK: " << do_steric_check << std::endl;
	std::cout << " TETHER ATOM ID: "  << tether_atom_id << std::endl;
	std::cout << " TETHER XYZ: "      << tether_xyz( 1) << ' ' << tether_xyz( 2 ) << ' ' << tether_xyz( 3 ) << std::endl;

	std::string silent_file = "";
	if ( option[ out::file::silent ].user() ) silent_file = option[ out::file::silent ]();

	Real const filter_rmsd = option[ filter_rms ]();
	clock_t const time_start( clock() );

	for ( Size i = 1; i <= N_SAMPLE; i++ ){
		Real const alpha = static_cast<Real>( i ) * alpha_increment + 0.01;
		std::cout << i << " out of " << N_SAMPLE << std::endl;

		for ( Size j = 1; j <= N_SAMPLE_COSBETA; j++ ){
			Real const cos_beta = -1.0 + static_cast< Real >( j ) * cos_beta_increment - 0.01;
			Real const beta = degrees( std::acos( cos_beta ) );

			for ( Size k = 1; k <= N_SAMPLE; k++ ){
				Real const gamma = static_cast< Real >( k ) * gamma_increment + 0.01;

				protocols::stepwise::create_euler_rotation( M, alpha, beta, gamma, axis1, axis2, axis3 );

				rotate( pose, M, pose_start, moving_res );

				Pose pose_to_translate = pose;

				bool found_better_orientation =	 search_translations( pose,
																															pose_to_translate,
																															moving_res,
																															partner_res,
																															tether_atom_id,
																															tether_xyz,
																															tether_radius,
																															filter_rmsd,
																															count,
																															do_steric_check,
																															silent_file,
																															native_pose,
																															rmsd_min,
																															best_pose );

				if ( found_better_orientation ){
					std::cout << "found better orientation. alpha, cos_beta, gamma: " << alpha << " " << cos_beta << " " << gamma << std::endl;
				}

			}
		}
	}

	std::cout << "Sampled number of poses: " << count << std::endl;

	std::cout << "Minimum rmsd: " << rmsd_min << std::endl;

	//	Real const rmsd = get_backbone_rmsd_no_super( best_pose, *native_pose );
	//	std::cout << "Check rmsd: " << rmsd << std::endl;
	best_pose.dump_pdb( "best.pdb" );


	std::cout << "Total time in spinner test: " <<
		static_cast<Real>(clock() - time_start) / CLOCKS_PER_SEC << std::endl;

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace core::options;

	if ( option[ spinner ] ){
		spinner_test();
	}

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );

}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace core::options;

	utility::vector1< Size > blank_size_vector;
	utility::vector1< std::string > blank_string_vector;

	NEW_OPT( spinner, "spinner", false );
	NEW_OPT( prepend, "prepend", false );
	NEW_OPT( append, "append", false );
	NEW_OPT( prepack, "prepack", false );
	NEW_OPT( sample_res, "sample residues", blank_size_vector );
	NEW_OPT( n_sample, "number of samples per torsion angle", 18 );
	NEW_OPT( n_sample_beta, "number of samples in tilt angle beta", 18 );
	NEW_OPT( xyz_sample, "spacing in xyz search, in Angstroms", 1.0 );
	NEW_OPT( filter_rms, "rmsd cut on moving segment", 0.0 );
	NEW_OPT( hydrophobic_cbeta_dist_cutoff, "how close cbetas need to be to define contact", 6.0 );
	NEW_OPT( steric_dist_cutoff, "how close cbetas need to be to define contact", 6.0 );
	NEW_OPT( num_loop_res, "number of intervening loop residues", 0 );
	NEW_OPT( min_hydrophobic_contacts, "minimum number of contacts", 2 );
	NEW_OPT( base, "pdb file for base pose", "" );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );

	exit( 0 );

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
