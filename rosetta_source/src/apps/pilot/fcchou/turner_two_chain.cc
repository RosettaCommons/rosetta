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
#include <core/init.hh>
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
#include <protocols/rna/RNA_IdealCoord.hh>
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
#include <boost/array.hpp>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>

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

static const scoring::rna::RNA_FittedTorsionInfo rna_fitted_torsion_info;
static numeric::random::RandomGenerator RG(5075);  // <- Magic number, do not change it!

typedef std::pair< core::Size, utility::vector1<float> > Trajectory_point;
typedef boost::array<float,5> Backbone_Torsion;
typedef boost::array<float,2> Nucleoside_Torsion;

OPT_KEY( String, force_field )
OPT_KEY( String, seq1 )
OPT_KEY( String, seq2 )
OPT_KEY( String, algorithm )
OPT_KEY( Integer, n_cycle )
OPT_KEY( String,  reference_rigid_body_samples_fixed_pair )
OPT_KEY( Integer, fixed_pair_state_number )
OPT_KEY( RealVector, kT_list )
OPT_KEY( RealVector, weight_list )
OPT_KEY( RealVector, input_torsion )
OPT_KEY( String, output_prefix )
OPT_KEY( Boolean, save_torsions )
OPT_KEY( Boolean, save_score_terms )

//////////////////////////////////
//Torsion angles setup and apply to pose.
//////////////////////////////////
std::pair < Backbone_Torsion, Nucleoside_Torsion >
ideal_A_form_torsions(){

	using namespace scoring::rna;

	std::pair < Backbone_Torsion, Nucleoside_Torsion > ideal_torsions;

	ideal_torsions.first[0] = rna_fitted_torsion_info.gaussian_parameter_set_epsilon_north()[1].center;
	ideal_torsions.first[1] = rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_minus()[1].center;
	ideal_torsions.first[2] = rna_fitted_torsion_info.gaussian_parameter_set_alpha()[1].center;
	ideal_torsions.first[3] = rna_fitted_torsion_info.gaussian_parameter_set_beta()[1].center;
	ideal_torsions.first[4] = rna_fitted_torsion_info.gaussian_parameter_set_gamma()[1].center;

	ideal_torsions.second[0] = rna_fitted_torsion_info.ideal_delta_north();
	ideal_torsions.second[1] = rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[1].center;
	
	return ideal_torsions;
}

//////////////////////////////////
void
apply_backbone( Backbone_Torsion const & backbone, 
								pose::Pose & pose,
								Size const suite ){

	using namespace id;
	using namespace scoring::rna;

	pose.set_torsion( TorsionID( suite  , id::BB, 5 ), backbone[0] ); //epsilon
	pose.set_torsion( TorsionID( suite  , id::BB, 6 ), backbone[1] ); //zeta
	pose.set_torsion( TorsionID( suite+1, id::BB, 1 ), backbone[2] ); //alpha
	pose.set_torsion( TorsionID( suite+1, id::BB, 2 ), backbone[3] ); //beta
	pose.set_torsion( TorsionID( suite+1, id::BB, 3 ), backbone[4] ); //gamma

}
//////////////////////////////////
void
apply_nucleoside( Nucleoside_Torsion const & nucleoside, 
								  pose::Pose & pose,
								  Size const residue ){

	using namespace id;
	using namespace scoring::rna;

	static protocols::rna::RNA_IdealCoord const ideal_coord_rna;

	if (nucleoside[0]  < 115) { //North
		if ( pose.torsion( TorsionID( residue, id::BB, 4 ) ) > 115) {
			ideal_coord_rna.apply(pose, residue, true);
		}
	} else { //South
		if ( pose.torsion( TorsionID( residue, id::BB, 4 ) ) < 115) {
			ideal_coord_rna.apply(pose, residue, false);
		}
	}
	
	pose.set_torsion( TorsionID( residue, id::CHI, 1 ), nucleoside[1] ); //chi
}
//////////////////////////////////
void
apply_all( utility::vector1< Backbone_Torsion > const & backbones,
					 utility::vector1< Nucleoside_Torsion > const & nucleosides,
					 pose::Pose & pose )
{
	Size i;
	for (i = 1; i <= backbones.size(); ++i) {
		apply_backbone( backbones[i], pose, i );
		apply_nucleoside( nucleosides[i], pose, i );
	}
	apply_nucleoside( nucleosides[i], pose, i );
}
/////////////////////////////////
//Create the starting model.
/////////////////////////////////
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
													 utility::vector1< Size > const & moving_res,
													 utility::vector1< Real > const & rbs) {

	using namespace protocols::swa;
	static Matrix M;
	static const Matrix reference_axes = Matrix::identity();
	static Vector const & axis1 = reference_axes.col_x();
	static Vector const & axis2 = reference_axes.col_y();
	static Vector const & axis3 = reference_axes.col_z();
	static Vector const reference_centroid = Vector( 0.0, 0.0, 0.0 );

	create_euler_rotation( M, rbs[1], rbs[2], rbs[3], axis1, axis2, axis3 );
	rotate( pose, M, pose_start, moving_res, reference_centroid );
	translate( pose, Vector( rbs[4], rbs[5], rbs[6] ), pose, moving_res );

}
//////////////////////////////////
void
setup_pose ( pose::Pose & pose){
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

	std::string const sequence1 = option[ seq1 ]();
	std::string const sequence2 = option[ seq2 ]();
	std::string const total_seq = sequence1 + sequence2;
	Size const total_len = total_seq.size();
	Size const len1      = sequence1.size();
	Size const len2      = sequence2.size();

		
	if (len1 == 0) utility_exit_with_message( "User need to specify at least -seq1 !!" ); 

	ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	pose::make_pose_from_sequence( pose, total_seq, *rsd_set );

	static std::pair < Backbone_Torsion,  Nucleoside_Torsion > const ideal_torsions = ideal_A_form_torsions();
	for ( Size i = 1; i <= total_len-1; i++ ){
		apply_backbone( ideal_torsions.first, pose, i );
		apply_nucleoside( ideal_torsions.second, pose, i );
	}
	apply_nucleoside( ideal_torsions.second, pose, total_len );
	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );

	if (len2 == 0) return; // Only one chain

	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", len1 + 1 );

	FoldTree f( total_len );
	f.new_jump( 1, total_len, len1 );
	f.set_jump_atoms( 1,
										scoring::rna::chi1_torsion_atom( pose.residue( 1 ) ),
										scoring::rna::chi1_torsion_atom( pose.residue( total_len ) )   );
	pose.fold_tree( f );

	utility::vector1< Size > strand1_res, strand2_res;
	for (Size i = 1; i <= len1; ++i) strand1_res.push_back(i);
	for (Size i = len1 + 1; i <= total_len; ++i) strand2_res.push_back(i);

	translate_and_rotate_residue_to_origin( pose, 1, strand1_res );
	translate_and_rotate_residue_to_origin( pose, total_len, strand2_res );

	//Need to setup starting base pair with user-input rigid body setting.
	utility::vector1< utility::vector1< Real > > reference_rigid_body_settings;
	std::string const infile_reference = option[ reference_rigid_body_samples_fixed_pair ]();
	read_rigid_body_settings( infile_reference, reference_rigid_body_settings );

	utility::vector1< Real > const & rbs = reference_rigid_body_settings[ option[ fixed_pair_state_number ]() ];
	apply_rigid_body_settings( pose, pose, strand2_res, rbs );
}
//////////////////////////////////
//Random angle sampling functions
/////////////////////////////////
void
random_angle(float & angle, Real const lower_bound = -180, Real const upper_bound = 180) {
	angle = RG.uniform() * (upper_bound - lower_bound) + lower_bound;
	if (angle < -180) {
		angle += 360;
	} else if (angle >= 180) {
		angle -= 360;
	}
}
/////////////////////////////////
void
random_delta(float & angle) {
	static const float delta_north = rna_fitted_torsion_info.ideal_delta_north();
	static const float delta_south = rna_fitted_torsion_info.ideal_delta_south();
	angle = (RG.uniform() < 0.5) ? delta_north : delta_south;
}
/////////////////////////////////
void
gaussian_angle_move(float & angle, Real const stdev) {
	angle += RG.gaussian() * stdev;
	if (angle < -180) {
		angle += 360;
	} else if (angle >= 180) {
		angle -= 360;
	}
}
/////////////////////////////////
void
nucleoside_sampling( Nucleoside_Torsion & nucleoside, 
										 Real const sample_stdev, 
										 Real const lowerbound,
										 Real const upperbound,
										 Real const pucker_prob ) 
{
	if (sample_stdev == 0) return; //Skip sampling

	static std::pair < Backbone_Torsion,  Nucleoside_Torsion > const ideal_torsions = ideal_A_form_torsions();

	if ( pucker_prob != 0 && (pucker_prob == 1 || RG.uniform() < pucker_prob) ) {
		random_delta( nucleoside[0] );
 	}

	if (sample_stdev < 0 || sample_stdev > 180 ) { //use random sampling for chi
		random_angle( nucleoside[1], lowerbound + ideal_torsions.second[1], upperbound + ideal_torsions.second[1] );
		return;
	}

	float chi = nucleoside[1];
	gaussian_angle_move( chi, sample_stdev );
	if ( lowerbound == -180.0 && upperbound == 180.0 ) {
		nucleoside[1] = chi;
	} else {
		float diff = chi - ideal_torsions.second[1];
		if (diff < -180) {
			diff += 360;
		} else if (diff >= 180) {
			diff -= 360;
		}

		if ( diff >= lowerbound && diff < upperbound ) nucleoside[1] = chi;
	}
}
/////////////////////////////////
void
backbone_sampling( Backbone_Torsion & backbone,
									 Real const sample_stdev, 
									 Real const lowerbound,
									 Real const upperbound )
{
	if (sample_stdev == 0) return; //Skip sampling
	static std::pair < Backbone_Torsion,  Nucleoside_Torsion > const ideal_torsions = ideal_A_form_torsions();

	if  (sample_stdev < 0 || sample_stdev > 180 ) { //use random sampling 
		for (Size i = 0; i != 5; ++i) random_angle( backbone[i], lowerbound + ideal_torsions.first[i], upperbound + ideal_torsions.first[i] );
		return;
	}
		
	for (Size i = 0; i != 5; ++i) {
		float angle = backbone[i];
		gaussian_angle_move( angle, sample_stdev );
		if ( lowerbound == -180.0 && upperbound == 180.0 ) {
			backbone[i] = angle;
		} else {
			float diff = angle - ideal_torsions.first[i];
			if (diff < -180) {
				diff += 360;
			} else if (diff >= 180) {
				diff -= 360;
			}
			if ( diff >= lowerbound && diff < upperbound ) backbone[i] = angle;
		}
	}
}
//////////////////////////////////
void
update_nucleoside( utility::vector1 < Nucleoside_Torsion > & nucleosides,
									 utility::vector1 < Size > const & sample_modes,
									 Real const stdev_bp,
									 Real const stdev_free,
									 Real const lowerbound_bp,
									 Real const upperbound_bp,
									 Real const pucker_prob )
{
	for (Size i = 1; i <= nucleosides.size(); ++i) {
		if ( sample_modes[i] == 1 ) { //free mode
			nucleoside_sampling( nucleosides[i], stdev_free, -180, 180, pucker_prob );
		} else if ( sample_modes[i] == 2 ) { //bp mode
			nucleoside_sampling( nucleosides[i], stdev_bp, lowerbound_bp, upperbound_bp, 0 );
		}
	}
}

//////////////////////////////////
void
update_backbone( utility::vector1 < Backbone_Torsion > & backbones,
								 utility::vector1 < Size > const & sample_modes,
								 Real const stdev_bp,
								 Real const stdev_free,
								 Real const lowerbound_bp,
								 Real const upperbound_bp )
								 
{
	for (Size i = 1; i <= backbones.size(); ++i) {
		if ( sample_modes[i] == 1 ) { //free mode
			backbone_sampling( backbones[i], stdev_free, -180, 180 );
		} else if ( sample_modes[i] == 2 ) { //bp mode
			backbone_sampling( backbones[i], stdev_bp, lowerbound_bp, upperbound_bp );
		}
	}
}
//////////////////////////////////
//Data saving
//////////////////////////////////
void
get_torsion_list( utility::vector1< float > & data,
									utility::vector1 < Backbone_Torsion > const & backbones,
									utility::vector1 < Nucleoside_Torsion > const & nucleosides,
									Real const missing_suite )
{
	for (Size i = 1; i <= backbones.size(); ++i) {
		if (i == missing_suite) continue;
		for (Size j = 0; j != 5; ++j) { 
			data.push_back( backbones[i][j] );
		}
	}

	for (Size i = 1; i <= nucleosides.size(); ++i) {
		for (Size j = 0; j != 2; ++j) { 
			data.push_back( nucleosides[i][j] );
		}
	}
}
//////////////////////////////////
void
get_score_terms( utility::vector1< float > & data,
								 pose::Pose & pose,
								 utility::vector1 < scoring::ScoreFunctionOP > const & scorefxns )
{
	for (Size i = 1; i <= scorefxns.size(); ++i) {
		data.push_back( (*scorefxns[i])(pose) );
	}
}
//////////////////////////////////
void
MC_run () {
	using namespace chemical;
	using namespace scoring;
	using namespace scoring::rna;
	using namespace kinematics;
	using namespace pose;
	using namespace utility::io;
	using namespace scoring::rna;

	//Setup scoring function
	std::string const force_field_name = option[force_field];
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( force_field_name );
	ScoreFunctionOP rep_scorefxn = new ScoreFunction;
	rep_scorefxn->set_weight( fa_rep, scorefxn -> get_weight( fa_rep ) );

	//Scorefxn list (hbond_sc, fa_stack, rna_torsion)
	utility::vector1< ScoreFunctionOP > scorefxns;
	ScoreFunctionOP hbond_sc_fxn = new ScoreFunction;
	hbond_sc_fxn -> set_weight( hbond_sc, 1.0 );
	scorefxns.push_back( hbond_sc_fxn );
	ScoreFunctionOP fa_stack_fxn = new ScoreFunction;
	fa_stack_fxn -> set_weight( fa_stack, 1.0 );
	scorefxns.push_back( fa_stack_fxn );
	ScoreFunctionOP rna_torsion_fxn = new ScoreFunction;
	rna_torsion_fxn -> set_weight( rna_torsion, 1.0 );
	scorefxns.push_back( rna_torsion_fxn );
	//

	Pose pose;
	setup_pose(pose);
	scorefxn -> show(pose);

	bool const is_save_torsions = option[ save_torsions ]();
	bool const is_save_score_terms = option[ save_score_terms ]();
	std::string const out_prefix = option[ output_prefix ]();
	Size const num_cycle = option[ n_cycle ]();
	utility::vector1< Real > const kTs = option[ kT_list ] ();
	utility::vector1< Real > weights = option[ weight_list ] ();
	if ( weights.empty() ) weights.push_back(0);
	if ( weights.size() != kTs.size() ) 
		utility_exit_with_message("kT_list and weight_list have different sizes!!!!!");
	Size const kT_id_max = kTs.size();

	std::string const sequence1 = option[ seq1 ]();
	std::string const sequence2 = option[ seq2 ]();
	std::string const total_seq = sequence1 + sequence2;
	Size const n_res = total_seq.size();
	Size const len1  = sequence1.size();
	Size const len2  = sequence2.size();

	std::pair < Backbone_Torsion,  Nucleoside_Torsion > const ideal_torsions = ideal_A_form_torsions();

	utility::vector1< Backbone_Torsion > backbones, backbones_new;
	utility::vector1< Nucleoside_Torsion > nucleosides, nucleosides_new;

	nucleosides.push_back( ideal_torsions.second );
	for (Size i = 1; i <= n_res-1; ++i) {
		backbones.push_back( ideal_torsions.first );
		nucleosides.push_back( ideal_torsions.second );
	}


	utility::vector1< Real > stdev_bp, stdev_free, pucker_prob;
	for (Size i = 1; i <= kT_id_max; ++i) {
		Real const kT = kTs[i];
		bool const is_kT_inf = (kT < 0);
		stdev_bp.push_back( (is_kT_inf) ? -1.0 : kT * 4.0 / double(n_res - 2) );
		stdev_free.push_back( (is_kT_inf) ? -1.0 : kT * 24.0 / double(n_res) );
		pucker_prob.push_back( (is_kT_inf) ? 1.0 : 0.2 ); 
	}

	Real const upperbound_bp = 60;
	Real const lowerbound_bp = -60;

	utility::vector1< Size > backbones_mode ( backbones.size(), 2 );
	utility::vector1< Size > nucleosides_mode ( nucleosides.size(), 2 );

	Real const n_bp = std::min(len1, len2);

	if (len1 > len2) {
		for (Size i = n_bp + 1; i <= len1; ++i) {
			nucleosides_mode[i] = 1;
			if (i-1 > 0) backbones_mode[i-1] = 1;
		}
		if (len1 <= backbones_mode.size())
			backbones_mode[len1] = 0;
	} else {
		for (Size i = len1 + 1; i <= n_res - n_bp; ++i) {
			nucleosides_mode[i] = 1;
			backbones_mode[i]   = 1;
		}
		backbones_mode[len1] = 0;
	}

	apply_all( backbones, nucleosides, pose );
	Real score = (*scorefxn)( pose );	
	Real score_new, rep_score;
	Size n_accpet = 0;
	Size n_exch = 0;
	Size n_acpt_exch = 0;
	Size kT_id = 1;
	Size kT_id_new = kT_id;
	Real kT = kTs[kT_id];
	Real kT_new = kT;
	Real log_prob = 0;
	Size current_count = 1;

	//Saved data
	utility::vector1< utility::vector1< Trajectory_point > > data_list;
	utility::vector1< float > data_point;
	for (Size i = 1; i <= kT_id_max; ++i) {
		utility::vector1< Trajectory_point > empty_data;
		data_list.push_back( empty_data );
	}

	Real const rep_score_cutoff = 100;
	Real const exchange_rate = 0.2;
	clock_t const time_start( clock() );

	///////////////////////////////////

	for (Size cycle = 1; cycle <= num_cycle; cycle++) {
		backbones_new = backbones;
		nucleosides_new = nucleosides;

		update_backbone( backbones_new, backbones_mode, stdev_bp[kT_id], 
										 stdev_free[kT_id], lowerbound_bp, upperbound_bp);
		update_nucleoside( nucleosides_new, nucleosides_mode,  stdev_bp[kT_id], 
											 stdev_free[kT_id], lowerbound_bp, upperbound_bp, pucker_prob[kT_id]);

		apply_all( backbones_new, nucleosides_new, pose );
	
		rep_score = (*rep_scorefxn)(pose);
		if ( rep_score > rep_score_cutoff ) {
			score_new = rep_score;
		} else {
			score_new = (*scorefxn)( pose );
		}

		if (kT < 0 || score_new < score || RG.uniform() < exp( (score - score_new) / kT )) {
			data_point.clear();
			data_point.push_back( score );
			if (is_save_torsions) get_torsion_list(data_point, backbones, nucleosides, len1);
			if (is_save_score_terms) get_score_terms(data_point, pose, scorefxns);
			data_list[kT_id].push_back( Trajectory_point (current_count, data_point) );
			
			current_count = 1;
			score = score_new;
			backbones = backbones_new;
			nucleosides = nucleosides_new;
			++n_accpet;
		} else {
			++current_count;
		}

		if (kT_id_max == 1) continue; //Only one kT, skip the T-jumping stage.

		if (RG.uniform() < exchange_rate) {
			++n_exch;
			if ( kT_id_max == 2 ) { //only 2 kTs
				kT_id_new = (kT_id == 1) ? 2 : 1;
			} else {
				kT_id_new = (RG.uniform() < 0.5) ? kT_id + 1 : kT_id - 1;
				if (kT_id_new < 1 || kT_id_new > kT_id_max) continue;
			}
			kT_new = kTs[kT_id_new];
			log_prob = -( (kT < 0) ? 0 : (score / kT) ) - ( (kT_new < 0) ? 0 : (score / kT_new) ) + weights[kT_id_new] - weights[kT_id];
			if (log_prob > 0 || RG.uniform() < exp (log_prob) ) { //check if we want to exchange kT
				++n_acpt_exch;
				kT_id = kT_id_new;
				kT = kT_new;

				data_point.clear();
				data_point.push_back( score );
				if (is_save_torsions) get_torsion_list(data_point, backbones, nucleosides, len1);
				if (is_save_score_terms) get_score_terms(data_point, pose, scorefxns);
				data_list[kT_id].push_back( Trajectory_point (current_count, data_point) ); 
				current_count = 1;
			}
		}
	}

	std::cout << "Total number of cycles = " << num_cycle << std::endl;
	std::cout << "Accept rate:" << (1.0 * n_accpet / num_cycle) << std::endl;
	if ( kT_id_max != 1) 
		std::cout << "Exchange Accept rate:" << (1.0 * n_acpt_exch / n_exch) << std::endl;

	Real const time_in_test = static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC;
	std::cout << "Time in sampler: " <<  time_in_test << std::endl;

	//Output data
	for (Size id = 1; id <= kT_id_max; ++id) {
		Real const kT = kTs[id];
		ozstream out;
		std::ostringstream oss;
		oss << out_prefix << '_' << std::fixed << std::setprecision(2) << kT << ".out.gz";
		out.open(oss.str());
		for (Size i = 1; i <= data_list[id].size(); ++i) {
			out << data_list[id][i].first << ' ';
			for (Size j = 1; j <= data_list[id][i].second.size(); ++j) {
				Real const & value = data_list[id][i].second[j];
				if ( value < 0.00005 && value > -0.00005 ) {
					out << "0 ";
				} else if ( value > 999 ) {
					out << "999 ";
				} else {
					out << std::fixed << std::setprecision( 4 ) << value << ' ';
				}
			}
			out << std::endl;
		}
		out.close();
	}
}
//////////////////////////////////
void
torsion2pdb()
{
	utility::vector1< Real > const all_torsions = option[ input_torsion ] ();

	Pose pose;
	setup_pose(pose);

	std::string const sequence1 = option[ seq1 ]();
	std::string const sequence2 = option[ seq2 ]();
	std::string const total_seq = sequence1 + sequence2;
	Size const n_res = total_seq.size();
	Size const len1  = sequence1.size();
	Size const len2  = sequence2.size();

	Size n_torsion = 2 * n_res;
	if (len2 == 0 ) {
		n_torsion += 5 * (n_res - 1);
	} else {
		n_torsion += 5 * (n_res - 2);
	}

	if ( all_torsions.size() != n_torsion ) utility_exit_with_message("Invalid number of torsions!!");

	std::pair < Backbone_Torsion,  Nucleoside_Torsion > const ideal_torsions = ideal_A_form_torsions();

	utility::vector1< Backbone_Torsion > backbones;
	utility::vector1< Nucleoside_Torsion > nucleosides;

	Size curr_position = 1;
	for (Size i = 1; i <= n_res-1; ++i) {
		if (i == len1) {
			backbones.push_back( ideal_torsions.first );
		} else {
			Backbone_Torsion backbone;
			for (Size j = curr_position; j <= curr_position + 4; ++j) backbone[j-curr_position] = all_torsions[j];
			backbones.push_back( backbone );
			curr_position += 5;
		}
	}
	for (Size i = 1; i <= n_res; ++i) {
			Nucleoside_Torsion nucleoside;
			for (Size j = curr_position; j <= curr_position + 1; ++j) nucleoside[j-curr_position] = all_torsions[j];
			nucleosides.push_back( nucleoside );
			curr_position += 2;
	}
	apply_all( backbones, nucleosides, pose );
	pose.dump_pdb("dump.pdb");
}

//////////////////////////////////
void*
my_main( void* )
{
	std::string const algorithm_name = option[algorithm];

	if ( algorithm_name == "" || algorithm_name == "MC" ) {
		MC_run();
	} else if ( algorithm_name == "torsion2pdb" ) {
		torsion2pdb();
	}
		
	exit( 0 );
}
//////////////////////////////////
int
main( int argc, char * argv [] )
{
	using namespace core;
	utility::vector1< Size > blank_size_vector;
	utility::vector1< Real > blank_size_vector_real;

	NEW_OPT( force_field, "score_file", "rna/rna_hires_fang");
	NEW_OPT( seq1, "sequence 1 to model, 3' to 5' ", "" );
	NEW_OPT( seq2, "sequence 2 to model, 3' to 5' ", "" );
	NEW_OPT( algorithm, "Specify algorithm to execute", "two_bp");
	NEW_OPT( n_cycle, "cycle number for Random sampling", 0 );
	NEW_OPT( reference_rigid_body_samples_fixed_pair, "input file with alpha, beta, gamma, x, y, and z", "" );
	NEW_OPT( fixed_pair_state_number, "from reference file, which rigid body setting to use", 1 );
	NEW_OPT( kT_list, "list of kT for ST", blank_size_vector_real );
	NEW_OPT( weight_list, "list of ST weights", blank_size_vector_real );
	NEW_OPT( input_torsion, "list of torsions", blank_size_vector_real );
	NEW_OPT( output_prefix, "prefix for the out file", "" );
	NEW_OPT( save_torsions, "save torsion angles", false );
	NEW_OPT( save_score_terms, "save score_terms", false );
	/////////////////////////////
	// setup
	//////////////////////////////
	core::init ( argc, argv );
	//////////////////////////////
	// end of setup
	//////////////////////////////

	protocols::viewer::viewer_main( my_main );
}
