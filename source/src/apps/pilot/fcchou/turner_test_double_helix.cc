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
#include <core/init/init.hh>
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <protocols/rna/RNA_SuiteAssign.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/Energies.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/io/silent/RNA_SilentStruct.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/swa/StepWiseClusterer.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <core/io/pdb/pose_io.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/gzip_util.hh>
#include <numeric/random/random.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace basic::options;
using ObjexxFCL::fmt::A;
using ObjexxFCL::fmt::I;
using ObjexxFCL::fmt::F;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using io::pdb::dump_pdb;
using utility::vector1;
using utility::tools::make_vector1;

typedef  numeric::xyzMatrix< Real > Matrix;

OPT_KEY( String, out_scores_prefix )
OPT_KEY( String, save_torsions )
OPT_KEY( String, force_field )
OPT_KEY( String, seq )
OPT_KEY( String, out_pdb )
OPT_KEY( String, out_fold )
OPT_KEY( String, out_unfold )
OPT_KEY( String, algorithm )
OPT_KEY( Real, kT )
OPT_KEY( Integer, n_cycle )
OPT_KEY( Boolean, o2star_trials )
OPT_KEY( Boolean, check_clash )
OPT_KEY( Boolean, sample_near_1a )
OPT_KEY( Integer, fixed_pair_state_number )
OPT_KEY( String,  reference_rigid_body_samples_fixed_pair )
OPT_KEY( RealVector, kT_list )
OPT_KEY( RealVector, ST_weight_list )

static const scoring::rna::RNA_FittedTorsionInfo rna_fitted_torsion_info;
static numeric::random::RandomGenerator RG(5075);  // <- Magic number, do not change it!

typedef std::pair< unsigned int, float[4]> RNA_scores; //count, total_score, hbond_sc, fa_stack, rna_torsion
typedef std::pair< core::Size, utility::vector1< utility::vector1< Real > > > Torsions;
/*
//////////Binary IO////////////////
void 
write_scores( utility::vector1 < RNA_scores > const & scores, std::string const & out_name) {
	std::ofstream out_file (out_name.c_str(), std::ios::out | std::ios::binary);
	Size const vector_size = scores.size();
	out_file.write( (char*) & vector_size, sizeof(vector_size) );
	out_file.write( (char*) & scores[1], sizeof(RNA_scores) * vector_size );
	out_file.close();
}

void
read_scores( utility::vector1 < RNA_scores > & scores, std::string const & in_name) {
	std::ifstream in_file (in_name.c_str(), std::ios::in | std::ios::binary);
	Size vector_size = 0;
	in_file.read( (char*) & vector_size, sizeof(vector_size) );
	scores.clear();
	scores.resize(vector_size);
	in_file.read( (char*) & scores[1], sizeof(RNA_scores) * vector_size );
	in_file.close();
}
*/

//////////////////////////////////
//////////////////////////////////
utility::vector1< Real >
get_suite_ideal_A_form_torsions(){

	using namespace scoring::rna;

	static utility::vector1< Real >  ideal_A_form_torsions;
	if (ideal_A_form_torsions.size() == 0) {
		ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_epsilon_north()[1].center );
		ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_minus()[1].center );
		ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_alpha()[1].center );
		ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_beta()[1].center );
		ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_gamma()[1].center );
		ideal_A_form_torsions.push_back( rna_fitted_torsion_info.ideal_delta_north() );
		ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[1].center );
		for (Size i = 1; i <= ideal_A_form_torsions.size(); ++i) {
			if (ideal_A_form_torsions[i] > 360) {
				ideal_A_form_torsions[i] -= 360;
			} else if (ideal_A_form_torsions[i] <=  0) {
				ideal_A_form_torsions[i] += 360;
			}	
		}
	}

	return ideal_A_form_torsions;
}
//////////////////////////////////
void
apply_suite_torsions( utility::vector1< Real > const & torsion_set, 
											pose::Pose & pose,
											Size const moving_suite,
											bool const sample_3prime_pucker = true){

	using namespace id;
	using namespace scoring::rna;

	pose.set_torsion( TorsionID( moving_suite, id::BB, 5 ), torsion_set[1] );   //epsilon
	pose.set_torsion( TorsionID( moving_suite, id::BB, 6 ), torsion_set[2] );   //zeta
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, 1 ), torsion_set[3] ); //alpha
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, 2 ), torsion_set[4] ); //beta
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, 3 ), torsion_set[5] ); //gamma

	Real delta, nu2, nu1;
	if (torsion_set[6] < 115) { //North pucker, [6] is delta angle (only pick one of the two states)
		delta = rna_fitted_torsion_info.ideal_delta_north();
		nu2 = rna_fitted_torsion_info.ideal_nu2_north();
		nu1 = rna_fitted_torsion_info.ideal_nu1_north();
	} else { //South pucker
		delta = rna_fitted_torsion_info.ideal_delta_south();
		nu2 = rna_fitted_torsion_info.ideal_nu2_south();
		nu1 = rna_fitted_torsion_info.ideal_nu1_south();
	}

	if (sample_3prime_pucker) {
		pose.set_torsion( TorsionID( moving_suite+1, id::BB,  4 ), delta );
		pose.set_torsion( TorsionID( moving_suite+1, id::CHI, 2 ), nu2 );
		pose.set_torsion( TorsionID( moving_suite+1, id::CHI, 3 ), nu1 );
		pose.set_torsion( TorsionID( moving_suite+1, id::CHI, 1 ), torsion_set[7] ); //chi
	} else {
		pose.set_torsion( TorsionID( moving_suite, id::BB,  4 ), delta );
		pose.set_torsion( TorsionID( moving_suite, id::CHI, 2 ), nu2 );
		pose.set_torsion( TorsionID( moving_suite, id::CHI, 3 ), nu1 );
		pose.set_torsion( TorsionID( moving_suite, id::CHI, 1 ), torsion_set[7] ); //chi
	}
}
//////////////////////////////////
void
initialize_o2star_pack( pose::Pose const & pose, 
												scoring::ScoreFunctionOP const scorefxn,
												scoring::ScoreFunctionOP o2star_pack_scorefxn,
												pack::task::PackerTaskOP o2star_pack_task ) {

	using namespace pack;
	using namespace pack::task;
	using namespace scoring;
	using namespace scoring::rna;

	for (Size i = 1; i <= pose.total_residue(); ++i) {
		o2star_pack_task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
		o2star_pack_task->nonconst_residue_task(i).or_ex4( true ); //extra rotamers?? Parin S. Jan 28, 2010
		o2star_pack_task->nonconst_residue_task(i).or_include_current( true );
	}
	// Each of the following terms have been pretty optimized for the packer (trie, etc.)
	o2star_pack_scorefxn->set_weight( fa_atr, scorefxn->get_weight( fa_atr ) );
	o2star_pack_scorefxn->set_weight( fa_rep, scorefxn->get_weight( fa_rep ) );
	o2star_pack_scorefxn->set_weight( hbond_lr_bb_sc, scorefxn->get_weight( hbond_lr_bb_sc ) );
	o2star_pack_scorefxn->set_weight( hbond_sr_bb_sc, scorefxn->get_weight( hbond_sr_bb_sc ) );
	o2star_pack_scorefxn->set_weight( hbond_sc, scorefxn->get_weight( hbond_sc ) );
	o2star_pack_scorefxn->set_energy_method_options( scorefxn->energy_method_options() );
	// note that geom_sol is not optimized well --> replace with lk_sol for now.
	o2star_pack_scorefxn->set_weight( fa_sol, scorefxn->get_weight( lk_nonpolar ) );
}
//////////////////////////////////
Real
create_random_angle_from_range(Real lower_bound = 0, Real upper_bound = 360) {
	return RG.uniform() * (upper_bound - lower_bound) + lower_bound;
}
//////////////////////////////////
Real
create_random_angle_from_range_list(utility::vector1< std::pair <Real, Real> > const & range_list) {
	Real range_sum = 0;
	utility::vector1< Real > range_edge_list;
	for (Size i = 1; i <= range_list.size(); ++i) {
		Real range_size = range_list[i].second - range_list[i].first;
		range_sum += range_size;
		range_edge_list.push_back(range_sum);
	}
	Real angle = RG.uniform() * range_sum;

	for (Size i = 1; i <= range_list.size(); ++i) {
		if (angle > range_edge_list[i]) continue;
		if (i != 1) angle -= range_edge_list[i-1];
		angle += range_list[i].first;
		break;
	}

	return angle;
}
//////////////////////////////////
void
create_random_torsions(utility::vector1< Real > & torsion_list, bool const modify_pucker = true) {
	static const Real delta_north = rna_fitted_torsion_info.ideal_delta_north();
	static const Real delta_south = rna_fitted_torsion_info.ideal_delta_south();

	//Sample all ranges
	const Real alpha   = create_random_angle_from_range();
	const Real beta    = create_random_angle_from_range();
	const Real gamma   = create_random_angle_from_range();
	const Real epsilon = create_random_angle_from_range();
	const Real zeta    = create_random_angle_from_range();
	const Real chi     = create_random_angle_from_range();
	const Real delta   = (modify_pucker && RG.uniform() < 0.5) ? delta_north : delta_south;

	torsion_list.clear();
	torsion_list.push_back(epsilon);
	torsion_list.push_back(zeta);
	torsion_list.push_back(alpha);
	torsion_list.push_back(beta);
	torsion_list.push_back(gamma);
	torsion_list.push_back(delta);
	torsion_list.push_back(chi);
}
//////////////////////////////////
void
sample_near_current_torsion(utility::vector1< Real > & torsion_list, 
														Real const stddev, 
														utility::vector1< Real > const & upper_bound,
														utility::vector1< Real > const & lower_bound,
														bool const modify_pucker = true, 
														bool const cyclic = true) {
	static const Real delta_north = rna_fitted_torsion_info.ideal_delta_north();
	static const Real delta_south = rna_fitted_torsion_info.ideal_delta_south();
	static Real new_torsion;

	for (Size i = 1; i <= 7; ++i) {
		if (i == 6) {
			if (modify_pucker && RG.uniform() < 0.2) {
				torsion_list[i]  = (RG.uniform() < 0.5) ? delta_north : delta_south;
			}

		} else {
			if (cyclic) {
				torsion_list[i] += RG.gaussian() * stddev;
				while (true) {
					if (torsion_list[i] > upper_bound[i]) {
						torsion_list[i] -= upper_bound[i];
						torsion_list[i] += lower_bound[i];
					} else if (torsion_list[i] <= lower_bound[i]) {
						torsion_list[i] += upper_bound[i];
						torsion_list[i] -= lower_bound[i];
					} else {
						break;
					}
				}

			} else {
				new_torsion = torsion_list[i] + RG.gaussian() * stddev;
				if (new_torsion <= upper_bound[i] && new_torsion > lower_bound[i]) {
					torsion_list[i] = new_torsion;
				}
			}
		}
	}
}
//////////////////////////////////
void
output_torsion_list( std::string const & outfile,
										 utility::vector1< utility::vector1< Real > > const & data_list,
										 Real const volume_factor){

	using namespace utility::io;

	ozstream out;
	out.open_append( outfile );

	for ( Size n = 1; n <= data_list.size(); n++ ){
		utility::vector1< Real > const & info = data_list[ n ];
		for ( Size i = 1; i <= info.size(); i++ ){
			out << info[i] << ' ';
		}
		out << volume_factor <<std::endl;
	}

	out.close();
}
//////////////////////////////////
void
read_rigid_body_settings( std::string const infile,
													utility::vector1< utility::vector1< Real > > & input_rigid_body_settings ){
		using namespace io::silent;
		using namespace utility::io;
		izstream input( infile );

		if ( !input ) {
			std::cerr << "No file: " << infile << std::endl;
			utility_exit_with_message( "No file" );
		}

		utility::vector1< Real >  vals;
		for ( Size i = 1; i <= 8; i++ ) vals.push_back( 0.0 );
		//Real alpha, beta, gamma, x, y, z, energy, volume;

		while ( input >> vals[1] ) {
			input >> vals[2] >> vals[3] >> vals[4] >> vals[5] >> vals[6] >> vals[7] >> vals[8] >> skip;
			input_rigid_body_settings.push_back( vals );
		}
}
//////////////////////////////////
void
apply_rigid_body_settings( pose::Pose & pose, 
													 pose::Pose const & pose_start,
													 utility::vector1< Size > moving_res,
													 Real const alpha,
													 Real const beta,
													 Real const gamma,
													 Real const x,
													 Real const y,
													 Real const z ){

	using namespace protocols::swa;
	static Matrix M;
	static const Matrix reference_axes = Matrix::identity();
	Vector const & axis1 = reference_axes.col_x();
	Vector const & axis2 = reference_axes.col_y();
	Vector const & axis3 = reference_axes.col_z();
	Vector const reference_centroid = Vector( 0.0, 0.0, 0.0 );

	create_euler_rotation( M, alpha, beta, gamma, axis1, axis2, axis3 );
	rotate( pose, M, pose_start, moving_res, reference_centroid );
	translate( pose, Vector( x,y,z), pose, moving_res );

}

//////////////////////////////////
void
setup_double_helix_pose ( pose::Pose & pose){
	using namespace chemical;
	using namespace kinematics;
	using namespace pose;
	using namespace io::silent;
	using namespace id;
	using namespace protocols::swa;
	using namespace protocols::swa::rna;
	using namespace scoring;
	using namespace scoring::rna;
	using namespace optimization;

	ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	std::string sequence = "";
	if ( option[ seq ].user() ) sequence = option[ seq ]();
	if ( sequence.size() % 2 != 0) {
		utility_exit_with_message( "Improper sequence length: should be multiples of 2" );
	}	

	pose::make_pose_from_sequence( pose, sequence, *rsd_set );
	Size n_res = pose.n_residue();
	Size chain_len = n_res/2;

	FoldTree f( n_res );
	f.new_jump( 1, n_res, chain_len );
	f.set_jump_atoms( 1,
										scoring::rna::chi1_torsion_atom( pose.residue( 1 ) ),
										scoring::rna::chi1_torsion_atom( pose.residue( n_res ) )   );
	pose.fold_tree( f );

	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );
	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", chain_len+1 );

	utility::vector1< Size > strand1_res, strand2_res;
	for (Size i = 1; i <= chain_len; ++i) {
		strand1_res.push_back(i);
		strand2_res.push_back(i + chain_len);
	}

	translate_and_rotate_residue_to_origin( pose, 1, strand1_res );
	translate_and_rotate_residue_to_origin( pose, n_res, strand2_res );

	//Need to setup starting base pair with user-input rigid body setting.
	utility::vector1< utility::vector1< Real > > reference_rigid_body_settings;
	std::string const infile_reference = option[ reference_rigid_body_samples_fixed_pair ]();
	std::cout << "Reading rigid body settings from " << infile_reference << std::endl;
	read_rigid_body_settings( infile_reference, reference_rigid_body_settings );
	utility::vector1< Real > const & rbs = reference_rigid_body_settings[ option[ fixed_pair_state_number ]() ];
	apply_rigid_body_settings( pose, pose, strand2_res, rbs[1],rbs[2],rbs[3],rbs[4],rbs[5],rbs[6] );

	// For now assume delta, chi are at ideal values.
	//Set the pucker and base-plane angle (delta and chi) to be ideal
	utility::vector1< Real > strand1_torsion_set = get_suite_ideal_A_form_torsions();
	for ( Size i = 1; i <= pose.total_residue() - 1; i++ ){
		apply_suite_torsions( strand1_torsion_set, pose, i );
	}
	apply_suite_torsions( strand1_torsion_set, pose, 1, false );

}
//////////////////////////////////
Size
score2bin(Real const score, Real const min_score, Real const max_score, Real const bin_size) {
	if (score > max_score || score < min_score) {
		utility_exit_with_message( "score2bin: score is larger than max_score!" );
	}
	return std::ceil( (score - min_score) / bin_size );
}
//////////////////////////////////
void
update_system(utility::vector1< Real > & torsion_list, Real const sampling_range) {
	static bool const sample_near_1a_ = option[ sample_near_1a ] ();
	bool random_sampl = false;
	if (sample_near_1a_) {
		if (sampling_range > 30) {
			random_sampl = true;
		}
	} else if (sampling_range > 300) {
		random_sampl = true;
	}

	static bool need_upper_lower_initiate = true;
	static utility::vector1< Real > const A_form_torsion = get_suite_ideal_A_form_torsions();
	static Size const torsion_list_size = A_form_torsion.size();	
	static Real delta_torsion;
	static utility::vector1< Real >  upper_bound, lower_bound;

	if (need_upper_lower_initiate) { 
		for ( Size i = 1; i <= torsion_list_size; ++i ) {
			if (sample_near_1a_) {
				upper_bound.push_back(A_form_torsion[i] + 40);
				lower_bound.push_back(A_form_torsion[i] - 40);
			} else {
				upper_bound.push_back(360);
				lower_bound.push_back(0);
			}
		}
		need_upper_lower_initiate = false;
	}

	if ( sample_near_1a_ ) { //Sample around 1a torsion (+-40)
		if (random_sampl) {
			for ( Size i = 1; i <= torsion_list_size; ++i ) {
				if (i != 6) {
					torsion_list[i] = create_random_angle_from_range(A_form_torsion[i] - 40, A_form_torsion[i] + 40);
				}
			}
		} else {
			sample_near_current_torsion(torsion_list, sampling_range, upper_bound, lower_bound, false, false);
		}
	} else {
		if (random_sampl) {
			create_random_torsions(torsion_list);
		} else {
			sample_near_current_torsion(torsion_list, sampling_range, upper_bound, lower_bound);
		}
	}
}
///////////////////////////////////
void
helix_minimize (core::pose::Pose & pose, scoring::ScoreFunctionOP scorefxn) {
	using namespace optimization;

	AtomTreeMinimizer minimizer;
	float const dummy_tol ( 0.00000001 );
	MinimizerOptions min_options1 ( "dfpmin", dummy_tol, false, false, false );

	kinematics::MoveMap mm;
	mm.set_bb ( true );
	mm.set_chi ( true );
	mm.set_jump ( false );

	minimizer.run ( pose, mm, *scorefxn, min_options1 );
}
//////////////////////////////////
Real
rmsd_compute (core::pose::Pose const & pose1, core::pose::Pose const & pose2) {
	Real sum_rmsd = 0;
	Size n_atom = 0;
	for (Size i = 1; i <= pose1.total_residue(); ++i) {
		conformation::Residue const & rsd1 = pose1.residue(i);
		conformation::Residue const & rsd2 = pose2.residue(i);
		for (Size j = 1; j <= rsd1.nheavyatoms(); ++j) {
			sum_rmsd += ( rsd1.xyz(j) - rsd2.xyz(j) ).length_squared();
			++n_atom;
		}
	}
	return sqrt(sum_rmsd / double(n_atom) );
}
///////////////////////////////////
bool 
is_atom_clash (core::pose::Pose const & pose, Real const dist_cutoff = 1.2) {
	Size const n_res = pose.total_residue();
	Real const dist_cutoff_sq = dist_cutoff * dist_cutoff;
	Real dist_sq, dist_x, dist_y, dist_z;
	for (Size i = 1; i <= n_res; ++i) {
		for (Size j = 1; j <= n_res; ++j) {
			if (i == j) continue;
			for (Size atom1 = 1; atom1 <= pose.residue(i).natoms(); ++atom1) {
				for (Size atom2 = 1; atom2 <= pose.residue(j).natoms(); ++atom2) {
					if ( pose.residue(i).is_virtual(atom1) || pose.residue(j).is_virtual(atom2) ) continue;
					Vector const & xyz1 = pose.residue(i).xyz(atom1);
					Vector const & xyz2 = pose.residue(j).xyz(atom2);
					dist_x = xyz1[0] - xyz2[0];
					dist_y = xyz1[1] - xyz2[1];
					dist_z = xyz1[2] - xyz2[2];
					if (dist_x > dist_cutoff) continue;
					if (dist_y > dist_cutoff) continue;
					if (dist_z > dist_cutoff) continue;

					dist_sq = dist_x * dist_x + dist_y * dist_y + dist_z * dist_z;
					if (dist_sq < dist_cutoff_sq) return true;
				}
			}
		}
	}
	return false;
}
//////////////////////////////////
Size
torsion2bin( Real const torsion, Size const total_bin) {
	Real torsion_test = torsion;
	if (torsion_test < 0) {
		torsion_test += 360;
	} else if (torsion > 360) {
		torsion_test -= 360;
	}

	return std::ceil( torsion_test / (360.0 / total_bin) );
}
//////////////////////////////////
bool
pose_list_compare( std::pair <pose::Pose, Real> const & i, std::pair <pose::Pose, Real> const & j) {
	return (i.second < j.second);
}
//////////////////////////////////
void
double_helix_test(){
	using namespace chemical;
	using namespace scoring;
	using namespace scoring::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace pack;
	using namespace pack::task;
	using namespace utility::io;
	using namespace scoring::rna;

	std::string output_pdb = option[out_pdb] ();

	//Setup scoring function
	std::string const force_field_name = option[force_field];
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( force_field_name );
	ScoreFunctionOP rep_scorefxn = new ScoreFunction;
	rep_scorefxn->set_weight( fa_rep, scorefxn -> get_weight( fa_rep ) );


	Pose pose;
	setup_double_helix_pose(pose);
	pose.dump_pdb( "ideal.pdb" );
	scorefxn -> show(pose);

	std::string const score_out_prefix = option[ out_scores_prefix ] ();
	std::string const torsion_list_out = option[ save_torsions ] ();
	std::string const outfile = option[ out::file::o ] ();
	Size const num_cycle = option[ n_cycle ]();
	Real const kT_sys = option[ kT ] ();
	bool const is_check_clash = option[ check_clash ] ();
	std::string const sequence = option[ seq ]();
	Size const n_res = pose.n_residue();

	// initialize for o2star rotamer trials. 	
	PackerTaskOP o2star_pack_task =  pack::task::TaskFactory::create_packer_task( pose );
	ScoreFunctionOP o2star_pack_scorefxn = new ScoreFunction;
	if ( option[ o2star_trials ]() ) {
		initialize_o2star_pack(pose, scorefxn, o2star_pack_scorefxn, o2star_pack_task );
	}

	utility::vector1< utility::vector1< Real > > suite_torsion, suite_torsion_new;

	utility::vector1< Real > A_form_torsion_set = get_suite_ideal_A_form_torsions();

	for (Size i = 1; i <= (n_res - 2); ++i) {
		suite_torsion.push_back(A_form_torsion_set);
	}
	suite_torsion_new = suite_torsion;

	utility::vector1< Size > hist;

	Real const max_score =  500.005;
	Real const min_score = -90.005;
	Real const bin_size  =   0.01;
	for (Real i = min_score; i <= max_score; i += bin_size) {
		hist.push_back(0);
	}

	bool const is_kT_inf = (kT_sys > 999);
	Real sampling_range = kT_sys * 4.0 / double(n_res - 2);
	Real score = (*scorefxn)( pose );	
	Real score_new, score_saved, rep_score;
	Size bin;
	Size n_accpet = 0;
	Size const n_res_half = (n_res - 2) / 2;
	Size const torsion_list_size = suite_torsion_new.size();
	Real const rep_score_cutoff = 100;
	clock_t const time_start( clock() );
	Real E_cum = 0;

	//Initialize binary score saving
	utility::vector1 <RNA_scores> scores_list;
	utility::vector1 <Torsions> torsions_list;
	bool is_save_torsions = ( torsion_list_out != "" );
	ScoreFunctionOP fa_stack_scorefxn = new ScoreFunction;
	ScoreFunctionOP hbond_sc_scorefxn = new ScoreFunction;
	ScoreFunctionOP hbond_intra_scorefxn = new ScoreFunction;
	ScoreFunctionOP rna_torsion_scorefxn = new ScoreFunction;
	fa_stack_scorefxn->set_weight( fa_stack, 1 );
	hbond_sc_scorefxn->set_weight( hbond_sc, 1 );
	hbond_intra_scorefxn->set_weight( hbond_intra, 1 );
	rna_torsion_scorefxn->set_weight( rna_torsion, 1 );

	RNA_scores current_scores;
	current_scores.first = 1;
	current_scores.second[0] = score;
	current_scores.second[1] = (*hbond_sc_scorefxn)(pose);
	current_scores.second[2] = (*fa_stack_scorefxn)(pose);
	current_scores.second[3] = (*rna_torsion_scorefxn)(pose);

	Torsions current_torsions;
	current_torsions.first = 1;
	current_torsions.second = suite_torsion;
	///////////////////////////////////

	for (Size cycle = 1; cycle <= num_cycle; cycle++) {
	
		for (Size i = 1; i <= torsion_list_size; ++i) {
			if (is_kT_inf) {
				update_system( suite_torsion_new[i], 9999 );
			} else {
				update_system( suite_torsion_new[i], sampling_range );
			}

			if (i <= n_res_half) {
				apply_suite_torsions(suite_torsion_new[i], pose, i, true);
			} else {
				apply_suite_torsions(suite_torsion_new[i], pose, i+1, false);
			}
		}

		rep_score = (*rep_scorefxn)(pose);
		if (is_check_clash && rep_score > rep_score_cutoff) {
				score_new = rep_score;
		} else {
			score_new = (*scorefxn)( pose );
		}

		if (is_kT_inf || score_new < score || RG.uniform() < exp( (score - score_new) / kT_sys )) {
			score = score_new;
			suite_torsion = suite_torsion_new;
			++n_accpet;

			scores_list.push_back(current_scores);
			if (is_save_torsions) torsions_list.push_back(current_torsions);

			current_scores.first = 1;
			current_scores.second[0] = score;
			current_scores.second[1] = (*hbond_sc_scorefxn)(pose);
			current_scores.second[2] = (*fa_stack_scorefxn)(pose);
			current_scores.second[3] = (*rna_torsion_scorefxn)(pose);

			current_torsions.first = 1;
			current_torsions.second = suite_torsion;

		} else {
			suite_torsion_new = suite_torsion;
			++current_scores.first;
			++current_torsions.first;
		}

		E_cum += score;
		score_saved = score;
		if (score_saved < min_score) {
			std::cout << "Warning: score < min_score" << std::endl;
			continue;
		} else if (score_saved > max_score) {
			score_saved = max_score;
		}
		bin = score2bin( score_saved, min_score, max_score, bin_size );
		++hist[bin];

	}

	std::cout << "Total number of rotamers applied: " << num_cycle << std::endl;
	std::cout << "accept rate:" << (1.0 * n_accpet / num_cycle) << std::endl;
	Real const time_in_test = static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC;
	std::cout << "Time in sampler: " <<  time_in_test << std::endl;

	Size first_bin = 0;
	Size last_bin  = 0;
	for (Size i = 1; i <= hist.size(); ++i) {
		if (first_bin == 0) {
			if (hist[i] == 0) continue;
			first_bin = i;
		}
		if (hist[i] != 0) last_bin = i;
	}
	
	ozstream out;
	out.open( outfile );
	out << "Score N_sample" << std::endl;
	out << "Total " << num_cycle << std::endl;
	out << "E_avg " << E_cum / num_cycle << std::endl;
	for (Size i = first_bin; i <= last_bin; ++i) {
		Real score = min_score + (static_cast<Real>(i)- 0.5) * static_cast<Real>(bin_size);
		out << score << " " << hist[i] << std::endl;
	}
	out.close();

	if (score_out_prefix != "") {
		ozstream out;
		std::ostringstream oss;
		oss << score_out_prefix << '_' << std::fixed << std::setprecision(2) << kT_sys << ".out.gz";
		out.open(oss.str());
		out << "count total hb stack torsion" << std::endl;
		for (Size i = 1; i <= scores_list.size(); ++i) {
			out << scores_list[i].first << ' ';
			for (Size j = 0; j != 4; ++j) {
				if ( scores_list[i].second[j] < 0.005 && scores_list[i].second[j] > -0.005 ) {
					out << "0 ";
				} else if ( scores_list[i].second[j] > 999 ) {
					out << "999 ";
				} else {
					out << std::fixed << std::setprecision( 2 ) << scores_list[i].second[j] << ' ';
				}
			}
			out << std::endl;
		}
		out.close();
	}

	if (is_save_torsions) {
		ozstream out;
		std::ostringstream oss;
		oss << torsion_list_out << '_' << std::fixed << std::setprecision(2) << kT_sys << ".out.gz";
		out.open(oss.str());
		for (Size i = 1; i <= torsions_list.size(); ++i) {
			out << torsions_list[i].first << ' ';
			for (Size j = 1; j <= torsions_list[i].second.size(); ++j) {
				for (Size k = 1; k <= torsions_list[i].second[j].size(); ++k) {
					out << std::fixed << std::setprecision( 2 ) << torsions_list[i].second[j][k] << ' ';
				}
			}
			out << std::endl;
		}
		out.close();
	}
}
//////////////////////////////////
void
helix_ST(){
	using namespace chemical;
	using namespace scoring;
	using namespace scoring::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace pack;
	using namespace pack::task;
	using namespace utility::io;
	using namespace scoring::rna;

	std::string output_pdb = option[out_pdb] ();

	//Setup scoring function
	std::string const force_field_name = option[force_field];
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( force_field_name );
	ScoreFunctionOP rep_scorefxn = new ScoreFunction;
	rep_scorefxn->set_weight( fa_rep, scorefxn -> get_weight( fa_rep ) );


	Pose pose;
	setup_double_helix_pose(pose);
	pose.dump_pdb( "ideal.pdb" );
	scorefxn -> show(pose);

	std::string const score_out_prefix = option[ out_scores_prefix ] ();
	std::string const outfile = option[ out::file::o ] ();
	Size const num_cycle = option[ n_cycle ]();
	utility::vector1< Real > const kT_sys_list = option[ kT_list ] ();
	utility::vector1< Real > const weight_list = option[ ST_weight_list ] ();
	bool const is_check_clash = option[ check_clash ] ();
	std::string const sequence = option[ seq ]();
	Size const n_res = pose.n_residue();

	if ( kT_sys_list.size() != weight_list.size() ) {
		utility_exit_with_message("kT_list and weight_list have different sizes!!!!!");
	}
	Size const n_temp =  kT_sys_list.size();
	std::cout << "kT list: ";
	for (Size i = 1; i <= n_temp; ++i) std::cout << kT_sys_list[i] << ' ';
	std::cout << std::endl;

	std::cout << "weight list: ";
	for (Size i = 1; i <= n_temp; ++i) std::cout << weight_list[i] << ' ';
	std::cout << std::endl;

	// initialize for o2star rotamer trials. 	
	PackerTaskOP o2star_pack_task =  pack::task::TaskFactory::create_packer_task( pose );
	ScoreFunctionOP o2star_pack_scorefxn = new ScoreFunction;
	if ( option[ o2star_trials ]() ) {
		initialize_o2star_pack(pose, scorefxn, o2star_pack_scorefxn, o2star_pack_task );
	}

	//suite torsion list
	utility::vector1< utility::vector1< Real > > suite_torsion, suite_torsion_new;
	utility::vector1< Real > const A_form_torsion_set = get_suite_ideal_A_form_torsions();
	for (Size i = 1; i <= (n_res - 2); ++i) {
		suite_torsion.push_back(A_form_torsion_set);
	}
	suite_torsion_new = suite_torsion;

	//histogram list
	utility::vector1< utility::vector1< Size > > hist;
	utility::vector1< Size > hist_1d;
	Real const max_score =  500.005;
	Real const min_score = -200.005;
	Real const bin_size  =   0.01;
	for (Real i = min_score; i <= max_score; i += bin_size) {
		hist_1d.push_back(0);
	}
	for (Size i = 1; i <= n_temp; ++i) {
		hist.push_back(hist_1d);
	}

	//score list
	ScoreFunctionOP fa_stack_scorefxn = new ScoreFunction;
	ScoreFunctionOP hbond_sc_scorefxn = new ScoreFunction;
	ScoreFunctionOP hbond_intra_scorefxn = new ScoreFunction;
	ScoreFunctionOP rna_torsion_scorefxn = new ScoreFunction;
	fa_stack_scorefxn->set_weight( fa_stack, 1 );
	hbond_sc_scorefxn->set_weight( hbond_sc, 1 );
	hbond_intra_scorefxn->set_weight( hbond_intra, 1 );
	rna_torsion_scorefxn->set_weight( rna_torsion, 1 );
	Real score = (*scorefxn)( pose );

	RNA_scores current_scores;
	current_scores.first = 1;
	current_scores.second[0] = score;
	current_scores.second[1] = (*hbond_sc_scorefxn)(pose);
	current_scores.second[2] = (*fa_stack_scorefxn)(pose);
	current_scores.second[3] = (*rna_torsion_scorefxn)(pose);

	utility::vector1< utility::vector1 < RNA_scores > > scores_list (n_temp);

	//Other params
	Size kT_id = 1;
	Size new_kT_id = 0;
	Real kT_sys = kT_sys_list[kT_id];
	bool is_kT_inf = (kT_sys > 999);
	Real sampling_range = kT_sys * 4.0 / double(n_res - 2);
	Real score_new, score_saved, new_kT, beta_old, beta_new,log_prob, rep_score;
	Size bin (0), n_accpet (0), n_exchange(0), n_accp_exch (0);
	Size const n_res_half = (n_res - 2) / 2;
	Size const torsion_list_size = suite_torsion_new.size();
	Real const rep_score_cutoff = 100;
	clock_t const time_start( clock() );
	Real lowest_score = 999;
	Pose lowest_pose = pose;
	
	//Start MC
	for (Size cycle = 1; cycle <= num_cycle; cycle++) {
		//Normal MCMC
		for (Size i = 1; i <= torsion_list_size; ++i) {
			if (is_kT_inf) {
				update_system( suite_torsion_new[i], 9999 );
			} else {
				update_system( suite_torsion_new[i], sampling_range );
			}

			if (i <= n_res_half) {
				apply_suite_torsions(suite_torsion_new[i], pose, i, true);
			} else {
				apply_suite_torsions(suite_torsion_new[i], pose, i+1, false);
			}
		}

		rep_score = (*rep_scorefxn)(pose);
		if (is_check_clash && rep_score > rep_score_cutoff) {
				score_new = rep_score;
		} else {
			score_new = (*scorefxn)( pose );
		}

		if (is_kT_inf || score_new < score || RG.uniform() < exp( (score - score_new) / kT_sys )) {
			score = score_new;
			suite_torsion = suite_torsion_new;
			++n_accpet;
			
			if (score < lowest_score) {
				lowest_score = score;
				lowest_pose = pose;
			}

			scores_list[kT_id].push_back(current_scores);
			current_scores.first = 1;
			current_scores.second[0] = score;
			current_scores.second[1] = (*hbond_sc_scorefxn)(pose);
			current_scores.second[2] = (*fa_stack_scorefxn)(pose);
			current_scores.second[3] = (*rna_torsion_scorefxn)(pose);

		} else {
			suite_torsion_new = suite_torsion;
			++current_scores.first;
		}

		score_saved = score;
		if (score_saved < min_score) {
			std::cout << "Warning: score < min_score" << std::endl;
			continue;
		} else if (score_saved > max_score) {
			score_saved = max_score;
		}

		bin = score2bin( score_saved, min_score, max_score, bin_size );
		++hist[kT_id][bin];

		//Try changing kT
		if (RG.uniform() < 0.2) {
			++n_exchange;
			(RG.uniform() < 0.5) ? new_kT_id = kT_id + 1 : new_kT_id = kT_id - 1;
			if (new_kT_id < 1 || new_kT_id > n_temp) continue;
			new_kT = kT_sys_list[new_kT_id];
			beta_old = (kT_sys > 999) ? 0 : (1.0 / kT_sys);
			beta_new = (new_kT > 999) ? 0 : (1.0 / new_kT);
			log_prob = - (beta_new - beta_old) * score + (weight_list[new_kT_id] - weight_list[kT_id]);
			if (log_prob > 0 || RG.uniform() < exp (log_prob) ) { //check if we want to exchange kT
				scores_list[kT_id].push_back(current_scores);
				current_scores.first = 1;
				++n_accp_exch;
				kT_id = new_kT_id;
				kT_sys = new_kT;
				is_kT_inf = (kT_sys > 999);
				sampling_range = kT_sys * 4.0 / double(n_res - 2);
			}
		}
	}

	std::cout << "Total number of rotamers applied: " << num_cycle << std::endl;
	std::cout << "accept rate:" << (1.0 * n_accpet / num_cycle) << std::endl;
	std::cout << "exchange rate:" << (1.0 * n_exchange / num_cycle) << std::endl;
	std::cout << "exchange accept rate:" << (1.0 * n_accp_exch / n_exchange) << std::endl;
	Real const time_in_test = static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC;
	std::cout << "Time in sampler: " <<  time_in_test << std::endl;

	for (Size id = 1; id <= hist.size(); ++id) {
		Size first_bin = 0;
		Size last_bin  = 0;
		for (Size i = 1; i <= hist[id].size(); ++i) {
			if (first_bin == 0) {
				if (hist[id][i] == 0) continue;
				first_bin = i;
			}
			if (hist[id][i] != 0) last_bin = i;
		}

		Real const kT = kT_sys_list[id];
		std::ostringstream oss;
		oss << outfile << '_' << std::fixed << std::setprecision(2) << kT << ".out";
	
		ozstream out;
		out.open( oss.str() );
		out << "Score N_sample" << std::endl;
		Size total_data = 0;
		for (Size i = first_bin; i <= last_bin; ++i) {
			Real score = min_score + (static_cast<Real>(i)- 0.5) * static_cast<Real>(bin_size);
			out << std::setprecision(10) << score << " " << hist[id][i] << std::endl;
			total_data += hist[id][i];
		}
		out << "Total " << total_data << std::endl;
		out.close();

		if (score_out_prefix != "") {
			ozstream out;
			std::ostringstream oss;
			oss << score_out_prefix << '_' << std::fixed << std::setprecision(2) << kT << ".out.gz";
			out.open(oss.str());
			out << "count total hb stack torsion" << std::endl;
			for (Size i = 1; i <= scores_list[id].size(); ++i) {
				out << scores_list[id][i].first << ' ';
				for (Size j = 0; j != 4; ++j) {
					if ( scores_list[id][i].second[j] < 0.005 && scores_list[id][i].second[j] > -0.005 ) {
						out << "0 ";
					} else if ( scores_list[id][i].second[j] > 999 ) {
						out << "999 ";
					} else {
						out << std::fixed << std::setprecision( 2 ) << scores_list[id][i].second[j] << ' ';
					}
				}
				out << std::endl;
			}
			out.close();
		}
	}

}
///////////////////////////////////
void
minimize_and_score () {
	using namespace chemical;
	using namespace scoring;
	using namespace scoring::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace pack;
	using namespace pack::task;
	using namespace utility::io;
	using namespace scoring::rna;

	//Setup scoring function
	std::string const force_field_name = option[force_field];
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( force_field_name );

	std::string output_pdb = option[out_pdb] ();

	Pose pose;
	setup_double_helix_pose(pose);

	helix_minimize(pose, scorefxn);
	scorefxn -> show ( std::cout, pose );
	pose.dump_pdb ( output_pdb );
}
///////////////////////////////////
void*
my_main( void* )
{
	std::string const algorithm_name = option[algorithm];

	if ( algorithm_name == "double_helix" ) {
		double_helix_test();
	} else if ( algorithm_name == "helix_ST" ) {
		helix_ST();
	} else if ( algorithm_name == "minimize_and_score" ) {
		minimize_and_score();
	}
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}
//////////////////////////////////
int
main( int argc, char * argv [] )
{
    try {
	using namespace core;
	utility::vector1< Size > blank_size_vector;
	utility::vector1< Real > blank_size_vector_real;

	NEW_OPT(out_scores_prefix, "", "");
	NEW_OPT(save_torsions, "", "");
	NEW_OPT(force_field, "score_file", "rna/rna_hires_fang");
	NEW_OPT( seq, "sequence to model", "" );
	NEW_OPT( out_pdb, "output pdb name", "test.pdb" );
	NEW_OPT( out_unfold, "name of output for unfolded states", "unfold.out" );
	NEW_OPT( out_fold, "name of output for folded states", "fold.out" );
	NEW_OPT( algorithm, "Specify algorithm to execute", "two_bp");
	NEW_OPT( kT, "", 9999.99 );
	NEW_OPT( n_cycle, "cycle number for Random sampling", 0 );
	NEW_OPT( o2star_trials, "in dinucleotide test, do rotamer trials", false );
	NEW_OPT( reference_rigid_body_samples_fixed_pair, "input file with alpha, beta, gamma, x, y, and z", "rbs_cluster.txt" );
	NEW_OPT( fixed_pair_state_number, "from reference file, which rigid body setting to use", 1 );
	NEW_OPT( kT_list, "list of kT for ST", blank_size_vector_real );
	NEW_OPT( ST_weight_list, "list of ST weights", blank_size_vector_real );
	NEW_OPT( check_clash, "check clahsed conformers", true );
	NEW_OPT( sample_near_1a, "sample near 1a conformation", true );

	/////////////////////////////
	// setup
	//////////////////////////////
	core::init::init ( argc, argv );
	//////////////////////////////
	// end of setup
	//////////////////////////////

	protocols::viewer::viewer_main( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
                                  }
        return 0;
    }
