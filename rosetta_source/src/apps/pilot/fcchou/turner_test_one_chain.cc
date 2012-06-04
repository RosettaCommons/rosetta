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

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
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
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/swa/StepWiseUtil.hh>
#include <protocols/swa/StepWiseClusterer.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <core/init.hh>
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
#include <time.h>

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
using utility::vector1;
using utility::tools::make_vector1;

OPT_KEY( String, force_field )
OPT_KEY( String, seq )
OPT_KEY( String, out_decoy )
OPT_KEY( String, algorithm )
OPT_KEY( Real, out_decoy_score_cutoff )
OPT_KEY( Real, score_cutoff )
OPT_KEY( Real, sampling_stddev )
OPT_KEY( Real, kT )
OPT_KEY( Integer, n_sample )
OPT_KEY( Boolean, check_clash )
OPT_KEY( Boolean, o2star_trials )
OPT_KEY( Boolean, sample_3_prime_end )
OPT_KEY( Boolean, add_virt_res )
OPT_KEY( RealVector, torsion_list )

static const scoring::rna::RNA_FittedTorsionInfo rna_fitted_torsion_info;
static numeric::random::RandomGenerator RG(245075);  // <- Magic number, do not change it!
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
	}
	return ideal_A_form_torsions;
}
//////////////////////////////////
void
apply_nucleoside_torsion( utility::vector1< Real > const & torsion_set,
													pose::Pose & pose,
													Size const moving_res){

	using namespace id;
	using namespace scoring::rna;

	Real delta, nu2, nu1;
	if (torsion_set[1] < 115) { //North pucker, [6] is delta angle (only pick one of the two states)
		delta = rna_fitted_torsion_info.ideal_delta_north();
		nu2 = rna_fitted_torsion_info.ideal_nu2_north();
		nu1 = rna_fitted_torsion_info.ideal_nu1_north();
	} else { //South pucker
		delta = rna_fitted_torsion_info.ideal_delta_south();
		nu2 = rna_fitted_torsion_info.ideal_nu2_south();
		nu1 = rna_fitted_torsion_info.ideal_nu1_south();
	}

	pose.set_torsion( TorsionID( moving_res, id::BB,  4 ), delta );
	pose.set_torsion( TorsionID( moving_res, id::CHI, 2 ), nu2 );
	pose.set_torsion( TorsionID( moving_res, id::CHI, 3 ), nu1 );
	pose.set_torsion( TorsionID( moving_res, id::CHI, 1 ), torsion_set[2] );
}
//////////////////////////////////
void
apply_suite_torsion( utility::vector1< Real > const & torsion_set,
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
setup_one_chain_pose ( pose::Pose & pose, bool is_virtualize = true ) {

	using namespace chemical;
	using namespace kinematics;
	using namespace pose;
	using namespace io::silent;
	using namespace id;
	using namespace protocols::swa;
	using namespace protocols::swa::rna;
	using namespace scoring::rna;

	ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	bool is_add_virt_res = option[ add_virt_res ]();
	std::string sequence = option[ seq ]();
	if (is_add_virt_res) sequence.insert(0, "a");

	pose::make_pose_from_sequence( pose, sequence, *rsd_set );

	if (is_add_virt_res && is_virtualize) {
		pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESIDUE", 1 );
		add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 2 );
	} else {
		pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );
	}

	//Set the pucker and base-plane angle (delta and chi) to be ideal
	utility::vector1< Real > A_form_torsion_set = get_suite_ideal_A_form_torsions();
	utility::vector1< Real > nucleoside1_torsion_set;
	nucleoside1_torsion_set.push_back(A_form_torsion_set[6]);
	nucleoside1_torsion_set.push_back(A_form_torsion_set[7]);

	for (Size i = 1; i <= pose.n_residue() - 1; ++i) {
		apply_suite_torsion( A_form_torsion_set, pose, i );
	}
	apply_nucleoside_torsion( nucleoside1_torsion_set, pose, 1 );
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
create_random_suite_torsion(utility::vector1< Real > & torsion_list) {
	static const Real delta_north = rna_fitted_torsion_info.ideal_delta_north();
	static const Real delta_south = rna_fitted_torsion_info.ideal_delta_south();

	//Sample all ranges
	const Real alpha   = create_random_angle_from_range();
	const Real beta    = create_random_angle_from_range();
	const Real gamma   = create_random_angle_from_range();
	const Real epsilon = create_random_angle_from_range();
	const Real zeta    = create_random_angle_from_range();
	const Real chi     = create_random_angle_from_range();
	const Real delta   = (RG.uniform() < 0.5) ? delta_north : delta_south;

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
create_random_nucleoside_torsion(utility::vector1< Real > & torsion_list) {
	static const Real delta_north = rna_fitted_torsion_info.ideal_delta_north();
	static const Real delta_south = rna_fitted_torsion_info.ideal_delta_south();

	//Sample all ranges
	const Real chi     = create_random_angle_from_range();
	const Real delta   = (RG.uniform() < 0.5) ? delta_north : delta_south;

	torsion_list.clear();
	torsion_list.push_back(delta);
	torsion_list.push_back(chi);
}
//////////////////////////////////
void
sample_near_suite_torsion(utility::vector1< Real > & torsion_list, Real const stddev) {
	static const Real delta_north = rna_fitted_torsion_info.ideal_delta_north();
	static const Real delta_south = rna_fitted_torsion_info.ideal_delta_south();

	torsion_list[1] += RG.gaussian() * stddev;
	torsion_list[2] += RG.gaussian() * stddev;
	torsion_list[3] += RG.gaussian() * stddev;
	torsion_list[4] += RG.gaussian() * stddev;
	torsion_list[5] += RG.gaussian() * stddev;
	torsion_list[7] += RG.gaussian() * stddev;
	if (RG.uniform() < 0.2) {
		torsion_list[6]  = (RG.uniform() < 0.5) ? delta_north : delta_south;
	}
	for (Size i = 0; i <= torsion_list.size(); ++i) {
		if (torsion_list[i] > 360) {
			torsion_list[i] -= 360;
		} else if (torsion_list[i] <=  0) {
			torsion_list[i] += 360;
		}
	}
}
//////////////////////////////////
void
sample_near_nucleoside_torsion(utility::vector1< Real > & torsion_list, Real const stddev) {
	static const Real delta_north = rna_fitted_torsion_info.ideal_delta_north();
	static const Real delta_south = rna_fitted_torsion_info.ideal_delta_south();

	if (RG.uniform() < 0.2) {
		torsion_list[1]  = (RG.uniform() < 0.5) ? delta_south : delta_north;
	}
	torsion_list[2] += RG.gaussian() * stddev;

	if (torsion_list[2] > 360) {
		torsion_list[2] -= 360;
	} else if (torsion_list[2] <=  0) {
		torsion_list[2] += 360;
	}

}

//////////////////////////////////
Size
score2bin(Real const score, Real const min_score, Real const max_score, Real const bin_size) {
	if (score > max_score || score < min_score) {
		utility_exit_with_message( "score2bin: score is larger than max_score!" );
	}
	return std::ceil( (score - min_score) / bin_size );
}


///////////////////////////////////////////////////////////////////////////////////////
void
save_decoy_to_silent_file_data( core::io::silent::SilentFileDataOP & silent_file_data,
																std::string const tag,
																pose::Pose const & pose, pose::Pose const & pose_reference ){
	using namespace core::io::silent;
	BinaryRNASilentStruct s( pose, tag );
	s.add_energy( "rms", core::scoring::all_atom_rmsd( pose, pose_reference ) );
	silent_file_data->add_structure( s );
}


///////////////////////////////////////////////////////////////////////////////////////
void
cluster_and_output_silent_file( core::io::silent::SilentFileDataOP & silent_file_data,
																 std::string const & silent_file ){
	using namespace protocols::swa;

	StepWiseClustererOP stepwise_clusterer = new StepWiseClusterer( silent_file_data );
	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer->set_max_decoys( max_decoys );
	stepwise_clusterer->set_rename_tags( true /*option[ rename_tags ]*/ );
	Real cluster_radius( 0.7 );
	if ( option[ OptionKeys::cluster::radius ].user() ) cluster_radius = option[ OptionKeys::cluster::radius ]();
	stepwise_clusterer->cluster();
	stepwise_clusterer->output_silent_file( silent_file );

}

//////////////////////////////////
void
one_chain_MC_sampling(){
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
	using namespace io::silent;

	Pose pose, pose_full;
	setup_one_chain_pose( pose );
	setup_one_chain_pose( pose_full, false );

	pose.dump_pdb("ideal.pdb");
	pose_full.dump_pdb("ideal_full.pdb");
	Size const n_residue = pose.n_residue();
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	Pose pose_start = pose;

	//Setup scoring function
	std::string const force_field_name = option[force_field];
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( force_field_name );
	ScoreFunctionOP rep_scorefxn = new ScoreFunction;
	rep_scorefxn->set_weight( fa_rep, scorefxn -> get_weight( fa_rep ) );
	scorefxn -> show(pose);

	std::string const outfile = option[ out::file::o ] ();
	std::string const outfile_decoy = option[out_decoy] ();
	Real const decoy_cutoff = option[ out_decoy_score_cutoff ]();
	Size const num_cycle = option[ n_sample ]();
	bool const is_check_clash = option[ check_clash ] ();
	bool const is_add_virt_res = option[ add_virt_res ] ();
	Real const kT_sys = option[ kT ] ();
	bool const is_saving_decoy = !(outfile_decoy == "");
	Size first_res = 0;
	if (is_add_virt_res) {
		first_res = 2;
	} else {
		first_res = 1;
	}

	Real sample_range;
	if (option[ sampling_stddev ].user()) {
		sample_range = option[ sampling_stddev ] ();
	} else {
		sample_range = kT_sys * 24.0 / double(n_residue);
	}

	// initialize for o2star rotamer trials.
	PackerTaskOP o2star_pack_task =  pack::task::TaskFactory::create_packer_task( pose );
	ScoreFunctionOP o2star_pack_scorefxn = new ScoreFunction;
	if ( option[ o2star_trials ]() ) {
		initialize_o2star_pack(pose, scorefxn, o2star_pack_scorefxn, o2star_pack_task );
	}

	// Get ready for main loop.
	clock_t const time_start( clock() );

	utility::vector1< Size > hist;
	utility::vector1< Real > nucleoside_torsion, nucleoside_torsion_new;
	utility::vector1< utility::vector1< Real > > suite_torsion, suite_torsion_new;

	Real const max_score =  500.005;
	Real const min_score = -90.005;
	Real const bin_size  =   0.01;
	for (Real i = min_score; i <= max_score; i += bin_size) {
		hist.push_back(0);
	}

	utility::vector1< Real > A_form_torsion_set = get_suite_ideal_A_form_torsions();
	nucleoside_torsion.push_back(A_form_torsion_set[6]);
	nucleoside_torsion.push_back(A_form_torsion_set[7]);
	nucleoside_torsion_new = nucleoside_torsion;
	apply_nucleoside_torsion(nucleoside_torsion, pose, first_res);

	for (Size i = first_res; i <= n_residue - 1; ++i) {
		apply_suite_torsion( A_form_torsion_set, pose, i );
		suite_torsion.push_back( A_form_torsion_set );
		suite_torsion_new.push_back( A_form_torsion_set );
	}
	Real score = (*scorefxn)( pose );
	pose.dump_pdb("test1.pdb");

	bool const is_kT_inf = (kT_sys > 999);
	Size n_accpet = 0;

	ozstream output_decoy;
	if (is_saving_decoy) {
		output_decoy.open(outfile_decoy);
	}
	SilentFileDataOP silent_file_data;
	std::string const silent_file = option[ out::file::silent ]();
	bool const save_silent = ( silent_file.size() > 0 );
	if ( save_silent ) silent_file_data = new SilentFileData;

	std::cout << "Start MC sampling" << std::endl;
	for (Size count = 1; count <= num_cycle; count++) {

		if (is_kT_inf) {
			for (Size i = first_res; i <= n_residue - 1; ++i) {
				create_random_suite_torsion(suite_torsion_new[i + 1 - first_res]);
				apply_suite_torsion( suite_torsion_new[i + 1 - first_res], pose, i );
			}
			create_random_nucleoside_torsion(nucleoside_torsion_new);
			apply_nucleoside_torsion(nucleoside_torsion_new, pose, first_res);
		} else {
			for (Size i = first_res; i <= n_residue - 1; ++i) {
				sample_near_suite_torsion(suite_torsion_new[i + 1 - first_res], sample_range);
				apply_suite_torsion( suite_torsion_new[i + 1 - first_res], pose, i );
			}
			sample_near_nucleoside_torsion(nucleoside_torsion_new, sample_range);
			apply_nucleoside_torsion(nucleoside_torsion_new, pose, first_res);
		}


//		if ( option[ o2star_trials ]() ) pack::rotamer_trials( pose, *o2star_pack_scorefxn, o2star_pack_task );

		Real const rep_score = (*rep_scorefxn)( pose_full );
		Real score_new;
		if (is_check_clash && rep_score > 1.3 * max_score) {
			score_new = 9999.99;
		} else {
			score_new = (*scorefxn)( pose );
		}

		if (is_kT_inf || score_new <= score || RG.uniform() < exp( (score - score_new) / kT_sys )) {
			score = score_new;
			nucleoside_torsion = nucleoside_torsion_new;
			suite_torsion = suite_torsion_new;
			++n_accpet;
			if (is_saving_decoy && score < decoy_cutoff) {
				output_decoy << score << ' ' << nucleoside_torsion[1] << ' ' << nucleoside_torsion[2];
				for (Size i = first_res; i <= first_res /*n_residue - 1*/ ; ++i) {
					for (Size j = 1; j <= suite_torsion[i].size(); ++j) {
						output_decoy << ' ' << suite_torsion[i][j];
					}
				}
				output_decoy << std::endl;

				std::string const tag = "S_"+lead_zero_string_of(n_accpet,6);
				if (save_silent) save_decoy_to_silent_file_data( silent_file_data, tag, pose, pose_start /*reference*/ );
			}
		} else {
			nucleoside_torsion_new = nucleoside_torsion;
			suite_torsion_new = suite_torsion;
		}

		Real score_saved = score;
		if (score_saved < min_score) {
			std::cout << "Warning: score < min_score" << std::endl;
			continue;
		} else if (score_saved > max_score) {
			score_saved = max_score;
		}

		Size const bin = score2bin( score_saved, min_score, max_score, bin_size );
		++hist[bin];

	}

	std::cout << "Total number of rotamers applied: " << num_cycle << std::endl;
	std::cout << "sampling stdev: " << sample_range << std::endl;
	std::cout << "accept rate:" << (1.0 * n_accpet / num_cycle) << std::endl;
	Real const time_in_dinucleotide_test = static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC;
	std::cout << "Time in dinucleotide sampler: " <<  time_in_dinucleotide_test << std::endl;

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
	for (Size i = first_bin; i <= last_bin; ++i) {
		Real score = min_score + (static_cast<Real>(i)- 0.5) * static_cast<Real>(bin_size);
		out << score << " " << hist[i] << std::endl;
	}
	out.close();

	if ( save_silent ) cluster_and_output_silent_file( silent_file_data, silent_file );

}

//////////////////////////////////
void*
my_main( void* ) {
	std::string const algorithm_name = option[algorithm];

	if ( algorithm_name == "one_chain_MC" ) {
		one_chain_MC_sampling();
	}

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}
//////////////////////////////////
int
main( int argc, char * argv [] ) {

	utility::vector1< Size > blank_size_vector;
	utility::vector1< Real > blank_size_vector_real;

	NEW_OPT(force_field, "score_file", "rna_hires_07232011_with_intra_base_phosphate");
	NEW_OPT( seq, "sequence to model", "" );
	NEW_OPT( out_decoy, "out file name for decoys", "" );
	NEW_OPT( out_decoy_score_cutoff, "out file name for decoys", 999.99 );
	NEW_OPT( algorithm, "Specify algorithm to execute", "");
	NEW_OPT( score_cutoff, "Cutoff in energy score", 99999.99 );
	NEW_OPT( n_sample, "Sample number for Random sampling", 0 );
	NEW_OPT( sampling_stddev, "Sampling standard deviation in degree", 0.0 );
	NEW_OPT( kT, "kT of simulation in RU", 9999.99 );
	NEW_OPT( check_clash, "check clahsed conformers", true );
	NEW_OPT( add_virt_res, "Append a virtual residue", true );
	NEW_OPT( o2star_trials, "do rotamer trials", false );
	//////////////////////////////
	// setup
	//////////////////////////////
	core::init(argc, argv);

	//////////////////////////////
	// end of setup
	//////////////////////////////

	protocols::viewer::viewer_main( my_main );
}
