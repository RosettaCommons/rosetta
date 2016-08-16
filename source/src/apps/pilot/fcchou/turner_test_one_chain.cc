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
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/chemical/rna/util.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
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
#include <core/io/silent/BinarySilentStruct.hh>
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
#include <protocols/stepwise/sampling/util.hh>
#include <protocols/farna/RNA_SuiteAssign.hh>
#include <protocols/stepwise/StepWiseLegacyClusterer.hh>
#include <protocols/stepwise/sampling/rna/util.hh>
#include <core/init/init.hh>
#include <core/io/pdb/pdb_writer.hh>
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

#include <utility/excn/Exceptions.hh>


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
using ObjexxFCL::format::A;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using utility::vector1;
using utility::tools::make_vector1;

OPT_KEY( String, out_scores_prefix )
OPT_KEY( String, force_field )
OPT_KEY( String, seq )
OPT_KEY( String, algorithm )
OPT_KEY( Real, score_cutoff )
OPT_KEY( Real, sampling_stddev )
OPT_KEY( Real, kT )
OPT_KEY( Integer, n_sample )
//OPT_KEY( Boolean, output_binary )
//OPT_KEY( String, output_binary_prefix )
OPT_KEY( Boolean, check_clash )
OPT_KEY( Boolean, o2prime_trials )
OPT_KEY( Boolean, sample_3_prime_end )
OPT_KEY( Boolean, add_virt_res )
OPT_KEY( Boolean, rmsd_vs_energy )
OPT_KEY( RealVector, torsion_list )
OPT_KEY( RealVector, kT_list )
OPT_KEY( RealVector, ST_weight_list )
OPT_KEY( Real, exchange_rate )

static const chemical::rna::RNA_FittedTorsionInfo rna_fitted_torsion_info;


//////////Binary IO////////////////

typedef std::pair<unsigned int, float[4]> RNA_scores; //count, total_score, hbond_sc, fa_stack, rna_torsion
/*
void
write_scores( utility::vector1 < RNA_scores > const & scores, std::string const & out_name) {
	std::ofstream out_file (out_name.c_str(), std::ios::out | std::ios::binary);
	Size const vector_size = scores.size();
	Size const RNA_scores_size = sizeof(RNA_scores);
	std::cout << "RNA_scores_size = " << RNA_scores_size << std::endl;
	out_file.write( (char*) & vector_size, sizeof(vector_size) );
	out_file.write( (char*) & RNA_scores_size, sizeof(RNA_scores_size) );
	out_file.write( (char*) & scores[1], RNA_scores_size * vector_size );
	out_file.close();
}

void
read_scores( utility::vector1 < RNA_scores > & scores, std::string const & in_name) {
	std::ifstream in_file (in_name.c_str(), std::ios::in | std::ios::binary);
	Size vector_size = 0;
	Size RNA_scores_size = 0;
	in_file.read( (char*) & vector_size, sizeof(vector_size) );
	in_file.read( (char*) & RNA_scores_size, sizeof(RNA_scores_size) );
	scores.clear();
	scores.resize(vector_size);
	std::cout << "RNA_scores_size = " << RNA_scores_size << std::endl;
	in_file.read( (char*) & scores[1], RNA_scores_size * vector_size );
	in_file.close();
}
*/
//////////////////////////////////
utility::vector1< Real >
get_suite_ideal_A_form_torsions(){

	using namespace chemical::rna;

	static utility::vector1< Real >  ideal_A_form_torsions;
	if (ideal_A_form_torsions.size() == 0) {
		ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_epsilon_north()[1].center );
		ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_minus()[1].center );
		ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_alpha()[1].center );
		ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_beta()[1].center );
		ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_gamma()[1].center );
		ideal_A_form_torsions.push_back( rna_fitted_torsion_info.delta_north() );
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
	using namespace chemical::rna;

	Real delta, nu2, nu1;
	if (torsion_set[1] < 115) { //North pucker, [6] is delta angle (only pick one of the two states)
		delta = rna_fitted_torsion_info.delta_north();
		nu2 = rna_fitted_torsion_info.nu2_north();
		nu1 = rna_fitted_torsion_info.nu1_north();
	} else { //South pucker
		delta = rna_fitted_torsion_info.delta_south();
		nu2 = rna_fitted_torsion_info.nu2_south();
		nu1 = rna_fitted_torsion_info.nu1_south();
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
	using namespace chemical::rna;

	pose.set_torsion( TorsionID( moving_suite, id::BB, 5 ), torsion_set[1] );   //epsilon
	pose.set_torsion( TorsionID( moving_suite, id::BB, 6 ), torsion_set[2] );   //zeta
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, 1 ), torsion_set[3] ); //alpha
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, 2 ), torsion_set[4] ); //beta
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, 3 ), torsion_set[5] ); //gamma

	Real delta, nu2, nu1;
	if (torsion_set[6] < 115) { //North pucker, [6] is delta angle (only pick one of the two states)
		delta = rna_fitted_torsion_info.delta_north();
		nu2 = rna_fitted_torsion_info.nu2_north();
		nu1 = rna_fitted_torsion_info.nu1_north();
	} else { //South pucker
		delta = rna_fitted_torsion_info.delta_south();
		nu2 = rna_fitted_torsion_info.nu2_south();
		nu1 = rna_fitted_torsion_info.nu1_south();
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
	using namespace protocols::stepwise::sampling::rna;
	using namespace chemical::rna;

	ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

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
initialize_o2prime_pack( pose::Pose const & pose,
												scoring::ScoreFunctionOP const scorefxn,
												scoring::ScoreFunctionOP o2prime_pack_scorefxn,
												pack::task::PackerTaskOP o2prime_pack_task ) {

	using namespace pack;
	using namespace pack::task;
	using namespace scoring;
	using namespace chemical::rna;

	for (Size i = 1; i <= pose.total_residue(); ++i) {
		o2prime_pack_task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
		o2prime_pack_task->nonconst_residue_task(i).or_ex4( true ); //extra rotamers?? Parin S. Jan 28, 2010
		o2prime_pack_task->nonconst_residue_task(i).or_include_current( true );
	}
	// Each of the following terms have been pretty optimized for the packer (trie, etc.)
	o2prime_pack_scorefxn->set_weight( fa_atr, scorefxn->get_weight( fa_atr ) );
	o2prime_pack_scorefxn->set_weight( fa_rep, scorefxn->get_weight( fa_rep ) );
	o2prime_pack_scorefxn->set_weight( hbond_lr_bb_sc, scorefxn->get_weight( hbond_lr_bb_sc ) );
	o2prime_pack_scorefxn->set_weight( hbond_sr_bb_sc, scorefxn->get_weight( hbond_sr_bb_sc ) );
	o2prime_pack_scorefxn->set_weight( hbond_sc, scorefxn->get_weight( hbond_sc ) );
	o2prime_pack_scorefxn->set_energy_method_options( scorefxn->energy_method_options() );
	// note that geom_sol is not optimized well --> replace with lk_sol for now.
	o2prime_pack_scorefxn->set_weight( fa_sol, scorefxn->get_weight( lk_nonpolar ) );
}
//////////////////////////////////
Real
create_random_angle_from_range(Real lower_bound = 0, Real upper_bound = 360) {
	return numeric::random::rg().uniform() * (upper_bound - lower_bound) + lower_bound;
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
	Real angle = numeric::random::rg().uniform() * range_sum;

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
	static const Real delta_north = rna_fitted_torsion_info.delta_north();
	static const Real delta_south = rna_fitted_torsion_info.delta_south();

	//Sample all ranges
	const Real alpha   = create_random_angle_from_range();
	const Real beta    = create_random_angle_from_range();
	const Real gamma   = create_random_angle_from_range();
	const Real epsilon = create_random_angle_from_range();
	const Real zeta    = create_random_angle_from_range();
	const Real chi     = create_random_angle_from_range();
	const Real delta   = (numeric::random::rg().uniform() < 0.5) ? delta_north : delta_south;

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
	static const Real delta_north = rna_fitted_torsion_info.delta_north();
	static const Real delta_south = rna_fitted_torsion_info.delta_south();

	//Sample all ranges
	const Real chi     = create_random_angle_from_range();
	const Real delta   = (numeric::random::rg().uniform() < 0.5) ? delta_north : delta_south;

	torsion_list.clear();
	torsion_list.push_back(delta);
	torsion_list.push_back(chi);
}
//////////////////////////////////
void
sample_near_suite_torsion(utility::vector1< Real > & torsion_list, Real const stddev) {
	static const Real delta_north = rna_fitted_torsion_info.delta_north();
	static const Real delta_south = rna_fitted_torsion_info.delta_south();

	torsion_list[1] += numeric::random::rg().gaussian() * stddev;
	torsion_list[2] += numeric::random::rg().gaussian() * stddev;
	torsion_list[3] += numeric::random::rg().gaussian() * stddev;
	torsion_list[4] += numeric::random::rg().gaussian() * stddev;
	torsion_list[5] += numeric::random::rg().gaussian() * stddev;
	torsion_list[7] += numeric::random::rg().gaussian() * stddev;
	if (numeric::random::rg().uniform() < 0.2) {
		torsion_list[6]  = (numeric::random::rg().uniform() < 0.5) ? delta_north : delta_south;
	}
	for (Size i = 1; i <= torsion_list.size(); ++i) {
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
	static const Real delta_north = rna_fitted_torsion_info.delta_north();
	static const Real delta_south = rna_fitted_torsion_info.delta_south();

	if (numeric::random::rg().uniform() < 0.2) {
		torsion_list[1]  = (numeric::random::rg().uniform() < 0.5) ? delta_south : delta_north;
	}
	torsion_list[2] += numeric::random::rg().gaussian() * stddev;

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
//////////////////////////////////
Real rmsd_compute(core::pose::Pose const & pose, core::pose::Pose const & ref_pose) {
	Size n_atom = 0;
	Real sd = 0;
	for (Size i = 1; i <= pose.total_residue(); ++i) {
		for (Size j = 1; j <= pose.residue(i).nheavyatoms(); ++j) {
			++n_atom;
			sd += ( pose.residue(i).xyz(j) - ref_pose.residue(i).xyz(j) ).length_squared();
		}
	}
	return sqrt(sd / double(n_atom) );
}

////////////////////////////////////
void
one_chain_MC_sampling(){
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace pack;
	using namespace pack::task;
	using namespace utility::io;
	using namespace chemical::rna;

	Pose pose, pose_full, lowest_pose;
	Real lowest_score = 9999;
	setup_one_chain_pose( pose );
	setup_one_chain_pose( pose_full, false );

	Pose const ideal_pose = pose;

	pose.dump_pdb("ideal.pdb");
	pose_full.dump_pdb("ideal_full.pdb");
	Size const n_residue = pose.n_residue();
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	//Setup scoring function
	std::string const force_field_name = option[force_field];
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( force_field_name );
	ScoreFunctionOP rep_scorefxn = new ScoreFunction;
	rep_scorefxn->set_weight( fa_rep, scorefxn -> get_weight( fa_rep ) );
	scorefxn -> show(pose);

	std::string const score_out_prefix = option[out_scores_prefix];
	std::string const outfile = option[ out::file::o ] ();
	Size const num_cycle = option[ n_sample ]();
	bool const is_check_clash = option[ check_clash ] ();
	bool const is_add_virt_res = option[ add_virt_res ] ();
	Real const kT_sys = option[ kT ] ();
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

	// initialize for o2prime rotamer trials.
	PackerTaskOP o2prime_pack_task =  pack::task::TaskFactory::create_packer_task( pose );
	ScoreFunctionOP o2prime_pack_scorefxn = new ScoreFunction;
	if ( option[ o2prime_trials ]() ) {
		initialize_o2prime_pack(pose, scorefxn, o2prime_pack_scorefxn, o2prime_pack_task );
	}

	// Get ready for main loop.
	clock_t const time_start( clock() );

	utility::vector1< Size > hist;
	utility::vector1< Real > nucleoside_torsion, nucleoside_torsion_new;
	utility::vector1< utility::vector1< Real > > suite_torsion, suite_torsion_new;

	Real const rep_score_cutoff = 100;
	Real const max_score =  500.005;
	Real const min_score = -200.005;
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

	//Initialize binary score saving
	utility::vector1 <RNA_scores> scores_list;
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


	bool const is_kT_inf = (kT_sys > 999);
	Size n_accpet = 0;

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


//		if ( option[ o2prime_trials ]() ) pack::rotamer_trials( pose, *o2prime_pack_scorefxn, o2prime_pack_task );

		Real const rep_score = (*rep_scorefxn)( pose_full );
		Real score_new;
		if (is_check_clash && rep_score > rep_score_cutoff) {
			score_new = rep_score;
		} else {
			score_new = (*scorefxn)( pose );
		}

		if (is_kT_inf || score_new <= score || numeric::random::rg().uniform() < exp( (score - score_new) / kT_sys )) {
			score = score_new;
			nucleoside_torsion = nucleoside_torsion_new;
			suite_torsion = suite_torsion_new;
			++n_accpet;

			scores_list.push_back(current_scores);
			current_scores.first = 1;
			current_scores.second[0] = score;
			current_scores.second[1] = (*hbond_sc_scorefxn)(pose);
			current_scores.second[2] = (*fa_stack_scorefxn)(pose);
			current_scores.second[3] = (*rna_torsion_scorefxn)(pose);

			if (score < lowest_score) {
				lowest_score = score;
				lowest_pose = pose;
			}
		} else {
			nucleoside_torsion_new = nucleoside_torsion;
			suite_torsion_new = suite_torsion;
			++current_scores.first;
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

	lowest_pose.dump_pdb("lowest.pdb");

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
				} else if ( scores_list[i].second[j] > 999 ){
					out << "999 ";
				} else {
					out << std::fixed << std::setprecision( 2 ) << scores_list[i].second[j] << ' ';
				}
			}
			out << std::endl;
		}
		out.close();
	}
}

//////////////////////////////////
void
one_chain_ST_MC () {
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace pack;
	using namespace pack::task;
	using namespace utility::io;
	using namespace chemical::rna;

	Pose pose, pose_full, lowest_pose;
	Real lowest_score = 9999;
	setup_one_chain_pose( pose );
	setup_one_chain_pose( pose_full, false );

	Pose const ideal_pose = pose;

	pose.dump_pdb("ideal.pdb");
	pose_full.dump_pdb("ideal_full.pdb");
	Size const n_residue = pose.n_residue();
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	//Setup scoring function
	std::string const force_field_name = option[force_field];
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( force_field_name );
	ScoreFunctionOP rep_scorefxn = new ScoreFunction;
	rep_scorefxn->set_weight( fa_rep, scorefxn -> get_weight( fa_rep ) );
	scorefxn -> show(pose);

	std::string const outfile = option[ out::file::o ] ();
	Size const num_cycle = option[ n_sample ]();
	bool const is_check_clash = option[ check_clash ] ();
	bool const is_add_virt_res = option[ add_virt_res ] ();

	Real const exchange_rate_ = option[ exchange_rate ]();
	utility::vector1< Real > const kT_sys_list = option[ kT_list ] ();
	utility::vector1< Real > const weight_list = option[ ST_weight_list ] ();
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


	Size first_res = 0;
	if (is_add_virt_res) {
		first_res = 2;
	} else {
		first_res = 1;
	}

	// initialize for o2prime rotamer trials.
	PackerTaskOP o2prime_pack_task =  pack::task::TaskFactory::create_packer_task( pose );
	ScoreFunctionOP o2prime_pack_scorefxn = new ScoreFunction;
	if ( option[ o2prime_trials ]() ) {
		initialize_o2prime_pack(pose, scorefxn, o2prime_pack_scorefxn, o2prime_pack_task );
	}

	// Get ready for main loop.
	clock_t const time_start( clock() );

	utility::vector1< Size > hist;
	utility::vector1< Real > nucleoside_torsion, nucleoside_torsion_new;
	utility::vector1< utility::vector1< Real > > suite_torsion, suite_torsion_new;

	Real const rep_score_cutoff = 100;
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

	Size kT_id = 1;
	Size new_kT_id = 0;
	Real kT_sys = kT_sys_list[kT_id];
	Real score_new, score_saved, new_kT, beta_old, beta_new,log_prob, rep_score;
	Size bin (0), n_accpet (0), n_exchange(0), n_accp_exch (0);
	bool is_kT_inf = (kT_sys > 999);
	Real sample_range = kT_sys * 24.0 / double(n_residue);


	ozstream output_rmsd;
	output_rmsd.open("rmsd_vs_energy.out");

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


//		if ( option[ o2prime_trials ]() ) pack::rotamer_trials( pose, *o2prime_pack_scorefxn, o2prime_pack_task );

		rep_score = (*rep_scorefxn)( pose_full );
		if (is_check_clash && rep_score > rep_score_cutoff) {
			score_new = rep_score;
		} else {
			score_new = (*scorefxn)( pose );
		}

		if (is_kT_inf || score_new <= score || numeric::random::rg().uniform() < exp( (score - score_new) / kT_sys )) {
			score = score_new;
			nucleoside_torsion = nucleoside_torsion_new;
			suite_torsion = suite_torsion_new;
			++n_accpet;

			if (score < lowest_score) {
				lowest_score = score;
				lowest_pose = pose;
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

		//Try changing kT
		if (numeric::random::rg().uniform() < exchange_rate_) {
			++n_exchange;
			(numeric::random::rg().uniform() < 0.5) ? new_kT_id = kT_id + 1 : new_kT_id = kT_id - 1;
			if (new_kT_id < 1 || new_kT_id > n_temp) continue;
			new_kT = kT_sys_list[new_kT_id];
			beta_old = (kT_sys > 999) ? 0 : (1.0 / kT_sys);
			beta_new = (new_kT > 999) ? 0 : (1.0 / new_kT);
			log_prob = - (beta_new - beta_old) * score + (weight_list[new_kT_id] - weight_list[kT_id]);
			if (log_prob > 0 || numeric::random::rg().uniform() < exp (log_prob) ) { //check if we want to exchange kT
				++n_accp_exch;
				kT_id = new_kT_id;
				kT_sys = new_kT;
				is_kT_inf = (kT_sys > 999);
				sample_range = kT_sys * 24.0 / double(n_residue);
			}
		}

	}

	std::cout << "Total number of rotamers applied: " << num_cycle << std::endl;
	std::cout << "sampling stdev: " << sample_range << std::endl;
	std::cout << "accept rate:" << (double(n_accpet) / num_cycle) << std::endl;
	std::cout << "exchange rate:" << (double(n_exchange) / num_cycle) << std::endl;
	std::cout << "exchange accept rate:" << (double(n_accp_exch) / n_exchange) << std::endl;
	std::cout << "Time in sampler: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

	lowest_pose.dump_pdb("lowest.pdb");

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
}
//////////////////////////////////
//////////////////////////////////


//////////////////////////////////
/*
void
one_chain_SWA_cluster(){
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace pack;
	using namespace pack::task;
	using namespace utility::io;
	using namespace chemical::rna;

	//////////////////////
	//SWA cluster centers
	Real cluster_center [9] [9] =
	{
	{ 84.15, 64.58, 212.55, 304.78, 297.25, 181.88, 59.32, 152.01, 122.36 },
	{ 84.37, 68.52, 216.09, 299.35, 290.05, 170.60, 53.62, 85.48, 79.70 },
	{ 84.19, 73.48, 204.91, 288.81, 154.12, 210.17, 165.99, 85.21, 67.51 },
	{ 84.38, 68.80, 196.64, 298.59, 159.83, 197.71, 179.87, 151.84, 91.93 },
	{ 84.34, 63.77, 202.80, 300.08, 295.52, 199.65, 53.62, 146.75, 308.71 },
	{ 84.43, 69.10, 203.71, 293.24, 155.14, 189.24, 173.08, 84.72, 246.89 },
	{ 85.06, 67.14, 201.27, 297.07, 159.67, 189.54, 186.98, 147.72, 272.36 },
	{ 84.51, 81.53, 209.94, 38.16, 151.99, 169.39, 52.11, 151.41, 117.12 },
	{ 83.83, 76.83, 217.93, 288.53, 297.70, 178.78, 60.59, 153.16, 111.76 }
	};
	Real lowest_score_list [9] = {999, 999, 999, 999, 999, 999, 999, 999, 999};
	Size counts_list [9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
	/////////////////////
	Pose pose, pose_full, lowest_pose;
	Real lowest_score = 9999;
	setup_one_chain_pose( pose );
	setup_one_chain_pose( pose_full, false );

	Pose const ideal_pose = pose;
	utility::vector1< Pose > lowest_pose_list;
	for (Size i = 1; i <= 9; ++i) {
		Pose pose_1 = pose;
		lowest_pose_list.push_back(pose_1);
	}

	pose.dump_pdb("ideal.pdb");
	pose_full.dump_pdb("ideal_full.pdb");
	Size const n_residue = pose.n_residue();
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	//Setup scoring function
	std::string const force_field_name = option[force_field];
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( force_field_name );
	ScoreFunctionOP rep_scorefxn = new ScoreFunction;
	rep_scorefxn->set_weight( fa_rep, scorefxn -> get_weight( fa_rep ) );
	scorefxn -> show(pose);

	std::string const outfile = option[ out::file::o ] ();
	Size const num_cycle = option[ n_sample ]();
	bool const is_check_clash = option[ check_clash ] ();
	bool const is_add_virt_res = option[ add_virt_res ] ();
	bool const output_rmsd_vs_energy = option[ rmsd_vs_energy ] ();
	Real const kT_sys = option[ kT ] ();
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

	// initialize for o2prime rotamer trials.
	PackerTaskOP o2prime_pack_task =  pack::task::TaskFactory::create_packer_task( pose );
	ScoreFunctionOP o2prime_pack_scorefxn = new ScoreFunction;
	if ( option[ o2prime_trials ]() ) {
		initialize_o2prime_pack(pose, scorefxn, o2prime_pack_scorefxn, o2prime_pack_task );
	}

	// Get ready for main loop.
	clock_t const time_start( clock() );

	utility::vector1< Size > hist;
	utility::vector1< Real > nucleoside_torsion, nucleoside_torsion_new;
	utility::vector1< utility::vector1< Real > > suite_torsion, suite_torsion_new;

	Real const rep_score_cutoff = 100;
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


//		if ( option[ o2prime_trials ]() ) pack::rotamer_trials( pose, *o2prime_pack_scorefxn, o2prime_pack_task );

		Real const rep_score = (*rep_scorefxn)( pose_full );
		Real score_new;
		if (is_check_clash && rep_score > rep_score_cutoff) {
			score_new = rep_score;
		} else {
			score_new = (*scorefxn)( pose );
		}

		if (is_kT_inf || score_new <= score || numeric::random::rg().uniform() < exp( (score - score_new) / kT_sys )) {
			score = score_new;
			nucleoside_torsion = nucleoside_torsion_new;
			suite_torsion = suite_torsion_new;
			++n_accpet;


			if (score < lowest_score) {
				lowest_score = score;
				lowest_pose = pose;
			}

			//////Find if near SWA clustering center///////
			Size n = 0;
			for (;n < 9; ++n) {
				Size m = 0;
				for (;m < 9; ++m) {
					Real angle;
					if (m < 2) {
						angle = nucleoside_torsion[m + 1];
					} else {
						angle = suite_torsion[1][m - 1];
					}
					//std::cout << angle << ' ' << cluster_center[n][m] << std::endl;
					if ( std::abs( angle - cluster_center[n][m] ) > 20 ) {
						//std::cout << std::endl;
						break;
					}
				}
				if (m == 9)	break;
			}

			if (n != 9) {
				++counts_list[n];
				if (score < lowest_score_list[n]) {
					lowest_score_list[n] = score;
					lowest_pose_list[n+1] = pose;
				}
			}
		  //////////////////////////////////////////////////
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

	lowest_pose.dump_pdb("lowest.pdb");

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

	ozstream out_SWA;
	out_SWA.open( "SWA_clustering.profile" );
	out_SWA << "counts lowest_E" << std::endl;
	for (Size i = 0; i < 9; ++i) {
		out_SWA << counts_list[i] << ' ' << lowest_score_list[i] << std::endl;
		std::ostringstream oss;
		oss << "SWA_cluster_lowest_"  << (i+1) << ".pdb";
		if (lowest_score_list[i] < 100) {
			lowest_pose_list[i+1].dump_pdb(oss.str());
		}
	}
	out_SWA.close();

}
*/
/////////////////////////////////////
//Convert torsion id into torsions
utility::vector1< Real > id2torsion(Size torsion_id) {
	utility::vector1< Real > torsion_list;
	//Ordering: epsilon, zeta, alpha, beta, gamma, chi1, chi2, delta1, delta2
	static utility::vector1< Real > const  A_form_torsion_set = get_suite_ideal_A_form_torsions();
	static const Real delta_north = rna_fitted_torsion_info.delta_north();
	static const Real delta_south = rna_fitted_torsion_info.delta_south();
	Size remainder;
	Real angle;
	for (Size i = 1; i <= 7; ++i) {
		remainder = torsion_id % 12;
		torsion_id = torsion_id / 12;
		if (i <= 5) {
			angle = A_form_torsion_set[i] + remainder * 30;
		} else {
			angle = A_form_torsion_set[7] + remainder * 30;
		}

		if (angle > 180) {
			angle -= 360;
		} else if (angle < -180) {
			angle += 360;
		}

		torsion_list.push_back(angle);
	}

	for (Size i = 1; i <= 2; ++i) {
		remainder = torsion_id % 2;
		torsion_id = torsion_id / 2;
		(remainder == 0) ? angle = delta_north : angle = delta_south;
		torsion_list.push_back(angle);
	}

	return torsion_list;
}
//Convert torsions in torsion_id (start with 0)
Size
torsion2id ( utility::vector1< Real > const & nucleoside_torsion,
             utility::vector1< Real > const & suite_torsion) {

	static utility::vector1< Real > const A_form_torsion = get_suite_ideal_A_form_torsions();
	static const Real delta_north = rna_fitted_torsion_info.delta_north();
	static const Real delta_south = rna_fitted_torsion_info.delta_south();
	Size id = 0;
	Size multiplier = 1;
	Real angle_diff;
	for (Size i = 1; i <= 5; ++i) {
		angle_diff = suite_torsion[i] - A_form_torsion[i] + 15;
		if (angle_diff > 360) {
			angle_diff -= 360;
		} else if (angle_diff < 0) {
			angle_diff += 360;
		}
		id += std::floor(angle_diff / 30.0) * multiplier;
		multiplier *= 12;
	}

	//chi1
	angle_diff = nucleoside_torsion[2] - A_form_torsion[7] + 15;
	if (angle_diff > 360) {
		angle_diff -= 360;
	} else if (angle_diff < 0) {
		angle_diff += 360;
	}
	id += std::floor(angle_diff / 30.0) * multiplier;
	multiplier *= 12;

	//chi2
	angle_diff = suite_torsion[7] - A_form_torsion[7] + 15;
	if (angle_diff > 360) {
		angle_diff -= 360;
	} else if (angle_diff < 0) {
		angle_diff += 360;
	}
	id += std::floor(angle_diff / 30.0) * multiplier;
	multiplier *= 12;

	//delta1
	if ( abs (nucleoside_torsion[1] - delta_south) < 1) id += multiplier;
	multiplier *= 2;

	//delta2
	if (abs (suite_torsion[6] - delta_south) < 1) id += multiplier;
	multiplier *= 2;

	return id;
}

//Convert torsions to 1d torsion_list
utility::vector1< Real >
torsion_convert ( utility::vector1< Real > const & nucleoside_torsion,
                  utility::vector1< Real > const & suite_torsion) {

	utility::vector1< Real >  torsion_list;
	for (Size i = 1; i <= 5; ++i) {
		torsion_list.push_back(suite_torsion[i]);
	}

	//chi1
	torsion_list.push_back(nucleoside_torsion[2]);

	//chi2
	torsion_list.push_back(suite_torsion[7]);

	//delta1
	torsion_list.push_back(nucleoside_torsion[1]);

	//delta2
	torsion_list.push_back(suite_torsion[6]);

	return torsion_list;
}


bool
sort_mine( std::pair <Size, Size> const i, std::pair <Size, Size> const j) {
	return (i.first > j.first);
}

/*
void test() {
	static utility::vector1< Real > const A_form_torsion = get_suite_ideal_A_form_torsions();
	static const Real delta_north = rna_fitted_torsion_info.delta_north();
	static const Real delta_south = rna_fitted_torsion_info.delta_south();
	utility::vector1< Real > nucleoside_torsion, suite_torsion;
	nucleoside_torsion.push_back(delta_north);
	nucleoside_torsion.push_back(155.22);
	suite_torsion.push_back(A_form_torsion[1]);
	suite_torsion.push_back(100.3);
	suite_torsion.push_back(300.3);
	suite_torsion.push_back(27.5);
	suite_torsion.push_back(98.33);
	suite_torsion.push_back(delta_south);
	suite_torsion.push_back(330.3);

	for (Size i = 1; i <= suite_torsion.size(); ++i) {
		std::cout << suite_torsion[i] << ' ';
	}
	for (Size i = 1; i <= nucleoside_torsion.size(); ++i) {
		std::cout << nucleoside_torsion[i] << ' ';
	}
	std::cout << std::endl;

	Size id = torsion2id(nucleoside_torsion, suite_torsion);
	std::cout << id << std::endl;
	utility::vector1< Real > torsion_list = id2torsion(id);
	for (Size i = 1; i <= torsion_list.size(); ++i) {
		std::cout << torsion_list[i] << ' ';
	}
	std::cout << std::endl;
}

void
one_chain_torsion_cluster(){
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace pack;
	using namespace pack::task;
	using namespace utility::io;
	using namespace chemical::rna;


	//Initialize one_chain_torsion_clustering
	utility::vector1< Real > lowest_score_torsion;
	utility::vector1< std::pair <Size, Size> > counts_torsion;
	utility::vector1< Size > used_id;
	utility::vector1< Size >::iterator id_it;
	utility::vector1< utility::vector1< Real > > lowest_torsion;
	Size n_id (0), curr_id(99999), pre_id(99999), curr_index(0);
	//////////////////////////////
	//Initialize suite test
	protocols::farna::RNA_suite_list const suite_list;
	utility::vector1 <protocols::farna::suite_info> const & all_suites = suite_list.full_list();
	utility::vector1 <std::string> suite_names;
	utility::vector1 <Size> suite_counts;
	for (Size i = 1; i <= all_suites.size(); ++i) {
		suite_names.push_back( all_suites[i].name );
		suite_counts.push_back( 0 );
	}
	suite_names.push_back( "!!" );
	suite_counts.push_back( 0 );
	Size suite_index = 1;
	//////////////////////////////////////////

	Pose pose, pose_full, lowest_pose;
	Real lowest_score = 9999;
	setup_one_chain_pose( pose );
	setup_one_chain_pose( pose_full, false );

	Pose const ideal_pose = pose;

	pose.dump_pdb("ideal.pdb");
	pose_full.dump_pdb("ideal_full.pdb");
	Size const n_residue = pose.n_residue();
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	//Setup scoring function
	std::string const force_field_name = option[force_field];
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( force_field_name );
	ScoreFunctionOP rep_scorefxn = new ScoreFunction;
	rep_scorefxn->set_weight( fa_rep, scorefxn -> get_weight( fa_rep ) );
	scorefxn -> show(pose);

	std::string const outfile = option[ out::file::o ] ();
	Size const num_cycle = option[ n_sample ]();
	bool const is_check_clash = option[ check_clash ] ();
	bool const is_add_virt_res = option[ add_virt_res ] ();
	bool const output_rmsd_vs_energy = option[ rmsd_vs_energy ] ();
	Real const kT_sys = option[ kT ] ();
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

	// initialize for o2prime rotamer trials.
	PackerTaskOP o2prime_pack_task =  pack::task::TaskFactory::create_packer_task( pose );
	ScoreFunctionOP o2prime_pack_scorefxn = new ScoreFunction;
	if ( option[ o2prime_trials ]() ) {
		initialize_o2prime_pack(pose, scorefxn, o2prime_pack_scorefxn, o2prime_pack_task );
	}

	// Get ready for main loop.
	clock_t const time_start( clock() );

	utility::vector1< Size > hist;
	utility::vector1< Real > nucleoside_torsion, nucleoside_torsion_new;
	utility::vector1< utility::vector1< Real > > suite_torsion, suite_torsion_new;

	Real const rep_score_cutoff = 100;
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


//		if ( option[ o2prime_trials ]() ) pack::rotamer_trials( pose, *o2prime_pack_scorefxn, o2prime_pack_task );

		Real const rep_score = (*rep_scorefxn)( pose_full );
		Real score_new;
		if (is_check_clash && rep_score > rep_score_cutoff) {
			score_new = rep_score;
		} else {
			score_new = (*scorefxn)( pose );
		}

		if (is_kT_inf || score_new <= score || numeric::random::rg().uniform() < exp( (score - score_new) / kT_sys )) {
			score = score_new;
			nucleoside_torsion = nucleoside_torsion_new;
			suite_torsion = suite_torsion_new;
			++n_accpet;

			//Calculate suite
			std::pair< std::string, std::pair <Size, Real> > suite = protocols::farna::suite_assign(pose, 2);
			suite_index = suite.second.first;
			if (suite_index == 0) suite_index = suite_names.size();
			++suite_counts[suite_index];

			//Cluster the torsions//
			curr_id = torsion2id(nucleoside_torsion, suite_torsion[1]);
			if (curr_id == pre_id) {
				counts_torsion[curr_index].first++;
				if (score < lowest_score_torsion[curr_index]) {
					lowest_score_torsion[curr_index] = score;
					lowest_torsion[curr_index] = torsion_convert(nucleoside_torsion, suite_torsion[1]);
				}
			} else {
				id_it = find (used_id.begin(), used_id.end(), curr_id);
				if ( id_it != used_id.end() ) {
					curr_index = std::distance(used_id.begin(), id_it) + 1;
					pre_id = curr_id;
					counts_torsion[curr_index].first++;
					if (score < lowest_score_torsion[curr_index]) {
						lowest_score_torsion[curr_index] = score;
						lowest_torsion[curr_index] = torsion_convert(nucleoside_torsion, suite_torsion[1]);
					}
				} else {
					++n_id;
					used_id.push_back(curr_id);
					counts_torsion.push_back( std::pair<Size, Size> (1, n_id) );
					lowest_score_torsion.push_back( score );
					lowest_torsion.push_back( torsion_convert(nucleoside_torsion, suite_torsion[1]) );
				}
			}

			////////////////////////
		} else {
			nucleoside_torsion_new = nucleoside_torsion;
			suite_torsion_new = suite_torsion;
			++suite_counts[suite_index];
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

	//output_torsion_clusters
	ozstream out_tor;
	out_tor.open( "torsion_cluster.out" );
	std::sort(counts_torsion.begin(), counts_torsion.end(), sort_mine );
	for (Size i = 1; i <= counts_torsion.size(); ++i) {
		out_tor << counts_torsion[i].first << " : ";
		Size const index = counts_torsion[i].second;
		out_tor << lowest_score_torsion[index] << " : ";
		utility::vector1< Real > cluster_center = id2torsion( used_id[index] );
		for (Size j = 1; j <= cluster_center.size(); ++j) {
			out_tor << cluster_center[j] << ' ';
		}
		out_tor << ": ";

		for (Size j = 1; j <= lowest_torsion[index].size(); ++j) {
			out_tor << lowest_torsion[index][j] << ' ';
		}
		out_tor << std::endl;
	}

	//output_suite_counts
	ozstream out_suite;
	out_suite.open( "suite_counts.out" );
	for (Size i = 1; i <= suite_names.size(); ++i) {
		out_suite << suite_names[i] << " : " << suite_counts[i] << std::endl;
	}
	out_suite.close();
}
*/
//////////////////////////////////
void
torsion2decoy () {
	utility::vector1< Real > const torsion = option[ torsion_list ] ();
	utility::vector1< Real > suite_torsion;
	utility::vector1< Real > nucleoside_torsion;

	core::pose::Pose pose;
	setup_one_chain_pose( pose );

	for (Size i = 1; i <= 5; ++i) {
		suite_torsion.push_back(torsion[i]);
	}
	suite_torsion.push_back(torsion[9]);
	suite_torsion.push_back(torsion[7]);
	nucleoside_torsion.push_back(torsion[8]);
	nucleoside_torsion.push_back(torsion[6]);
	apply_nucleoside_torsion(nucleoside_torsion, pose, 1);
	apply_suite_torsion( suite_torsion, pose, 1 );
	std::pair< std::string, std::pair <Size, Real> > suite = protocols::farna::suite_assign(pose, 2);
	std::cout << suite.first << ' ' << suite.second.second << std::endl;
	pose.dump_pdb("decoy.pdb");
}
/*
//////////////////////////////////
void
unpack_binary () {
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace pack;
	using namespace pack::task;
	using namespace utility::io;
	using namespace chemical::rna;

	std::string const outfile = option[ out::file::o ] ();
	std::string const infile  = option[  in::file::s ] () [1];
	utility::vector1 <RNA_scores> scores_list;

	read_scores( scores_list, infile);

	ozstream out;
	out.open( outfile );

	for (Size i = 1; i <= scores_list.size(); ++i) {
		RNA_scores const & scores = scores_list[i];
		out << scores.first << ' ' << scores.second[0] << ' ' << scores.second[1] << ' ' << scores.second[2] << ' ' <<
		scores.second[3] << std::endl;
	}
}
*/
//////////////////////////////////
void*
my_main( void* ) {
	std::string const algorithm_name = option[algorithm];


	if ( algorithm_name == "one_chain_MC" ) {
		one_chain_MC_sampling();
	} else if ( algorithm_name == "one_chain_ST" ) {
		one_chain_ST_MC();
/*
  } else if ( algorithm_name == "one_chain_cluster_SWA_test" ) {
		one_chain_SWA_cluster();

  } else if ( algorithm_name == "one_chain_torsion_cluster" ) {
		one_chain_torsion_cluster();

	} else if ( algorithm_name == "torsion2decoy" ) {
		torsion2decoy();
	} else if ( algorithm_name == "unpack_binary" ) {
		unpack_binary();
*/
	}

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}
//////////////////////////////////
int
main( int argc, char * argv [] ) {
    try {
	utility::vector1< Size > blank_size_vector;
	utility::vector1< Real > blank_size_vector_real;

	NEW_OPT(out_scores_prefix, "output the scores", "");
	NEW_OPT(force_field, "score_file", "farna/rna_hires_07232011_with_intra_base_phosphate");
	NEW_OPT( seq, "sequence to model", "" );
	NEW_OPT( algorithm, "Specify algorithm to execute", "");
	NEW_OPT( score_cutoff, "Cutoff in energy score", 99999.99 );
	NEW_OPT( n_sample, "Sample number for Random sampling", 0 );
	NEW_OPT( sampling_stddev, "Sampling standard deviation in degree", 0.0 );
	NEW_OPT( kT, "kT of simulation in RU", 9999.99 );
	NEW_OPT( check_clash, "check clahsed conformers", false );
	//NEW_OPT( output_binary, "", false );
	//NEW_OPT( output_binary_prefix, "", "" );
	NEW_OPT( rmsd_vs_energy, "", true );
	NEW_OPT( add_virt_res, "Append a virtual residue", false );
	NEW_OPT( o2prime_trials, "do rotamer trials", false );
	NEW_OPT( kT_list, "list of kT for ST", blank_size_vector_real );
	NEW_OPT( ST_weight_list, "list of ST weights", blank_size_vector_real );
	NEW_OPT( torsion_list, "list of torsions", blank_size_vector_real );
	NEW_OPT( exchange_rate, "", 0.2 );
	//////////////////////////////
	// setup
	//////////////////////////////
	core::init::init(argc, argv);

	//////////////////////////////
	// end of setup
	//////////////////////////////

	protocols::viewer::viewer_main( my_main );
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;

    }
