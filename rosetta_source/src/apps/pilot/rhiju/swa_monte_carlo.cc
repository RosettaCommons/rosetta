// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file swa_monte_Carlo.cc
/// @author Rhiju Das (rhiju@stanford.edu)

// libRosetta headers
#include <core/types.hh>

// do we need all these?
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/rms_util.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/types.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/init.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
///////////////////////////////////////////////////
#include <protocols/idealize/idealize.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <core/pose/MiniPose.hh>
#include <core/pose/util.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/import_pose.hh>

//////////////////////////////////////////////////////////
#include <protocols/viewer/viewers.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <protocols/rna/RNA_LoopCloser.hh>
#include <protocols/swa/rna/StepWiseRNA_OutputData.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_PoseSetup.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_PoseSetup.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters_Setup.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>  //Test
#include <cctype>
#include <iomanip>
#include <map>
#include <cstdlib>
#include <ctime>
#include <unistd.h>


#include <list>
#include <stdio.h>
#include <math.h>

using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

// SWA Monte Carlo -- July 2, 2012 -- Rhiju Das
//
// TO DO
//
//  Clean up pose setup, inherited from SWA code -- all these options and setup functions
//    should go into their own .cc file.
//
//  Generalize addition/deletion code to handle ends with 3' termini -- code
//    below only handles building back from 3' termini (for historical reasons)
//
//  Add screen for building new base [base atr/rep] -- like SWA
//
//  Set up GAAA tetraloop run [should be faster test case], with buildup
//   from either end.
//
//  Set up chainbreak variants when all residues are formed.
//
//  Set up loop closer move (perhaps just analytical loop close move?)
//
//  Set up constraints  when gap is 1,2, etc.
//
//  encapsulate -- move into a namespace? lay out plans for others?
//


// A lot of these options should be placed into an 'official' namespace
// and the SWA RNA pose setup should go to its own function
OPT_KEY( Boolean, do_not_sample_multiple_virtual_sugar)
OPT_KEY( Boolean, sample_ONLY_multiple_virtual_sugar)
OPT_KEY( Boolean, skip_sampling)
OPT_KEY( Boolean, skip_clustering)
OPT_KEY( Boolean, minimize_and_score_native_pose)
OPT_KEY( Integer, num_pose_minimize)
OPT_KEY( Boolean, combine_long_loop_mode)
OPT_KEY( Integer, job_queue_ID)
OPT_KEY( String, filter_output_filename)
OPT_KEY( Boolean, filter_for_previous_contact)
OPT_KEY( Boolean, filter_for_previous_clash)
OPT_KEY( Boolean, filterer_undercount_ribose_rotamers)
OPT_KEY( Boolean, combine_helical_silent_file)
OPT_KEY( Boolean, exclude_alpha_beta_gamma_sampling)
OPT_KEY( Boolean, debug_eplison_south_sugar_mode)
OPT_KEY( Boolean, rebuild_bulge_mode)
OPT_KEY( Boolean, sampler_include_torsion_value_in_tag)
OPT_KEY( Boolean, sampler_extra_anti_chi_rotamer)
OPT_KEY( Boolean, sampler_extra_syn_chi_rotamer)
OPT_KEY( Boolean, sampler_extra_beta_rotamer)
OPT_KEY( Boolean, sampler_extra_epsilon_rotamer)
OPT_KEY( Boolean, sample_both_sugar_base_rotamer)
OPT_KEY( Boolean, reinitialize_CCD_torsions)
OPT_KEY( Boolean, PBP_clustering_at_chain_closure)
OPT_KEY( Boolean, finer_sampling_at_chain_closure)
OPT_KEY( StringVector, 	VDW_rep_screen_info)
OPT_KEY( Real, 	VDW_rep_alignment_RMSD_CUTOFF)
OPT_KEY( Boolean, graphic )
OPT_KEY( Real, Real_parameter_one )
OPT_KEY( Boolean, add_lead_zero_to_tag )
OPT_KEY( Boolean, distinguish_pucker )
OPT_KEY( Boolean, include_syn_chi  )
OPT_KEY( Boolean, sampler_allow_syn_pyrimidine )
OPT_KEY( IntegerVector, native_virtual_res  )
OPT_KEY( Real, whole_struct_cluster_radius )
OPT_KEY( Real, suite_cluster_radius )
OPT_KEY( Real, loop_cluster_radius )
OPT_KEY( StringVector, alignment_res )
OPT_KEY( Boolean, filter_user_alignment_res )
OPT_KEY( IntegerVector, native_alignment_res )
OPT_KEY( StringVector, jump_point_pairs )
OPT_KEY( Boolean, floating_base )
OPT_KEY( Boolean, parin_favorite_output )
OPT_KEY( IntegerVector, sample_res )
OPT_KEY( IntegerVector, input_res )
OPT_KEY( IntegerVector, input_res2 )
OPT_KEY( IntegerVector, missing_res )
OPT_KEY( IntegerVector, missing_res2 )
OPT_KEY( IntegerVector, cutpoint_open )
OPT_KEY( Integer, cutpoint_closed )
OPT_KEY( IntegerVector, fixed_res )
OPT_KEY( IntegerVector, minimize_res )
OPT_KEY( IntegerVector, virtual_res )
OPT_KEY( IntegerVector, bulge_res )
OPT_KEY( IntegerVector, terminal_res )
OPT_KEY( IntegerVector, rmsd_res )
OPT_KEY( Boolean, centroid_screen )
OPT_KEY( Boolean, allow_base_pair_only_centroid_screen )
OPT_KEY( Boolean, VDW_atr_rep_screen )
OPT_KEY( Boolean, sampler_perform_o2star_pack )
OPT_KEY( Boolean, fast )
OPT_KEY( Boolean, medium_fast )
OPT_KEY( Boolean, allow_bulge_at_chainbreak )
OPT_KEY( Boolean, VERBOSE )
OPT_KEY( Boolean, sampler_native_rmsd_screen )
OPT_KEY( Real, sampler_native_screen_rmsd_cutoff )
OPT_KEY( Real, score_diff_cut )
OPT_KEY( Boolean, clusterer_perform_score_diff_cut )
OPT_KEY( String, 	algorithm)
OPT_KEY( Integer, sampler_num_pose_kept)
OPT_KEY( Integer, clusterer_num_pose_kept)
OPT_KEY( Boolean, recreate_silent_struct )
OPT_KEY( Boolean, clusterer_use_best_neighboring_shift_RMSD )
OPT_KEY( Boolean, allow_chain_boundary_jump_partner_right_at_fixed_BP )
OPT_KEY( Boolean, allow_fixed_res_at_moving_res )
OPT_KEY( Boolean, clusterer_rename_tags )
OPT_KEY( Boolean, simple_append_map )
OPT_KEY( IntegerVector, global_sample_res_list )
OPT_KEY( IntegerVector, force_syn_chi_res_list )
OPT_KEY( IntegerVector, force_north_ribose_list )
OPT_KEY( IntegerVector, force_south_ribose_list )
OPT_KEY( IntegerVector, protonated_H1_adenosine_list )
OPT_KEY( Boolean,  output_pdb )
OPT_KEY( String, 	start_silent)
OPT_KEY( String, 	start_tag)
OPT_KEY( Boolean,  simple_full_length_job_params )
OPT_KEY( Real, sampler_cluster_rmsd )
OPT_KEY( Boolean, 	output_extra_RMSDs)
OPT_KEY( Boolean, 	integration_test)
OPT_KEY( Boolean, 	add_virt_root ) //For Fang's electron density code.

OPT_KEY( Real, stddev_small )
OPT_KEY( Real, stddev_large )
OPT_KEY( Real, output_score_cutoff )
OPT_KEY( Real, kT )
OPT_KEY( Integer, n_sample )
OPT_KEY( Boolean, skip_randomize );
OPT_KEY( Boolean, sample_all_o2star );
OPT_KEY( Boolean, do_add_delete );
OPT_KEY( Boolean, presample_added_residue );
OPT_KEY( Integer, presample_internal_cycles );
OPT_KEY( Boolean, forward_build ); //temporary -- delete this soon!


//////////////////////////////////////////////////
// copied from fang's turner_test_one_chain.cc
static const scoring::rna::RNA_FittedTorsionInfo rna_fitted_torsion_info;
static numeric::random::RandomGenerator RG(245099);  // <- Magic number, do not change it!

void
sample_near_suite_torsion(utility::vector1< Real > & torsion_list, Real const stddev) {
	static const Real delta_north = rna_fitted_torsion_info.ideal_delta_north();
	static const Real delta_south = rna_fitted_torsion_info.ideal_delta_south();

	torsion_list[1] += RG.gaussian() * stddev;
	torsion_list[2] += RG.gaussian() * stddev;
	torsion_list[3] += RG.gaussian() * stddev;
	torsion_list[4] += RG.gaussian() * stddev;
	torsion_list[5] += RG.gaussian() * stddev;

}


//////////////////////////////////////////////////
// copied from fang's turner_test_one_chain.cc
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
										 Size const moving_suite ){

	using namespace id;
	using namespace scoring::rna;

	pose.set_torsion( TorsionID( moving_suite, id::BB, 5 ), torsion_set[1] );   //epsilon
	pose.set_torsion( TorsionID( moving_suite, id::BB, 6 ), torsion_set[2] );   //zeta
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, 1 ), torsion_set[3] ); //alpha
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, 2 ), torsion_set[4] ); //beta
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, 3 ), torsion_set[5] ); //gamma

}


//////////////////////////////////
void
apply_nucleoside_torsion_Aform(
													pose::Pose & pose,
													Size const moving_res ){

	utility::vector1< Real > ideal_A_form_torsions;
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info.ideal_delta_north() );
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[1].center );
	apply_nucleoside_torsion( ideal_A_form_torsions, pose, moving_res );
}
//////////////////////////////////
void
apply_suite_torsion_Aform(
													pose::Pose & pose,
													Size const moving_suite ){

	utility::vector1< Real > ideal_A_form_torsions;
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_epsilon_north()[1].center );
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_minus()[1].center );
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_alpha()[1].center );
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_beta()[1].center );
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_gamma()[1].center );
	apply_suite_torsion( ideal_A_form_torsions, pose, moving_suite );
}


///////////////////////////////////////////////////
utility::vector1< Real>
get_suite_torsion( pose::Pose const & pose, Size const moving_suite ){

	using namespace id;

	utility::vector1< Real > torsion_set;
	torsion_set.push_back(	pose.torsion( TorsionID( moving_suite, id::BB, 5 ) ) );   //epsilon
	torsion_set.push_back(	pose.torsion( TorsionID( moving_suite, id::BB, 6 ) ) );   //zeta
	torsion_set.push_back(	pose.torsion( TorsionID( moving_suite+1, id::BB, 1 ) ) ); //alpha
	torsion_set.push_back(	pose.torsion( TorsionID( moving_suite+1, id::BB, 2 ) ) ); //beta
	torsion_set.push_back(	pose.torsion( TorsionID( moving_suite+1, id::BB, 3 ) ) ); //gamma

	return torsion_set;
}

///////////////////////////////////////////////////
utility::vector1< Real>
get_nucleoside_torsion( pose::Pose const & pose, Size const moving_nucleoside ){

	using namespace id;

	utility::vector1< Real > torsion_set;
	torsion_set.push_back(	pose.torsion( TorsionID( moving_nucleoside, id::BB, 4 ) ) );  //delta
	torsion_set.push_back(	pose.torsion( TorsionID( moving_nucleoside, id::CHI, 1 ) ) ); //chi

	return torsion_set;
}
///////////////////////////////////////////////////
void
sample_near_suite_torsion( pose::Pose & pose, Size const moving_suite, Real const sample_range){
	utility::vector1< Real> torsion_set = get_suite_torsion( pose, moving_suite );
	sample_near_suite_torsion( torsion_set, sample_range );
	apply_suite_torsion( torsion_set, pose, moving_suite );
}
///////////////////////////////////////////////////
void
crankshaft_alpha_gamma( pose::Pose & pose, Size const moving_suite, Real const sample_range){

	using namespace id;

	TorsionID alpha_torsion_id( moving_suite+1, id::BB, 1 );
	TorsionID gamma_torsion_id( moving_suite+1, id::BB, 3 );

	Real alpha = pose.torsion( alpha_torsion_id );
	Real gamma = pose.torsion( gamma_torsion_id );
	Real const perturb = RG.gaussian() * sample_range;

	alpha += perturb;
	gamma -= perturb;

	pose.set_torsion( alpha_torsion_id, alpha );
	pose.set_torsion( gamma_torsion_id, gamma );

}
///////////////////////////////////////////////////
void
sample_near_nucleoside_torsion( pose::Pose & pose, Size const moving_res, Real const sample_range){
	utility::vector1< Real> torsion_set = get_nucleoside_torsion( pose, moving_res );
	sample_near_nucleoside_torsion( torsion_set, sample_range);
	apply_nucleoside_torsion( torsion_set, pose, moving_res );
}

///////////////////////////////////////////////////
void
sample_near_o2star_torsion( pose::Pose & pose, Size const moving_res, Real const sample_range){
	id::TorsionID o2star_torsion_id( moving_res, id::CHI, 4 );
	Real o2star_torsion = pose.torsion( o2star_torsion_id ); //get
	o2star_torsion += RG.gaussian() * sample_range; //sample_near
	pose.set_torsion( o2star_torsion_id, o2star_torsion ); // apply
}



//////////////////////////////////////////////////////////////////////////////////////
Size const
get_random_o2star_residue( pose::Pose & pose ){
	// pick at random from whole pose -- a quick initial stab.
	Size const o2star_num = int( pose.total_residue() * RG.uniform() ) + 1;
	return o2star_num;
}

//////////////////////////////////////////////////////////////////////////////////////
// This could be made smarter -- could go over nucleoside *and* suite.
Size const
get_random_o2star_residue_near_moving_residue( pose::Pose & pose, utility::vector1< Size > const moving_res_list ){

	// should be better -- actually look at o2star's that might be engaged in interactions with moving nucleoside
	utility::vector1< bool > residues_allowed_to_be_packed( pose.total_residue(), false );
	scoring::EnergyGraph const & energy_graph( pose.energies().energy_graph() );

	Distance DIST_CUTOFF( 4.0 );

	for (Size k = 1; k <= moving_res_list.size(); k++ ){
		Size const i = moving_res_list[ k ];

		for( graph::Graph::EdgeListConstIter
					 iter = energy_graph.get_node( i )->const_edge_list_begin();
				 iter != energy_graph.get_node( i )->const_edge_list_end();
				 ++iter ){

			Size j( (*iter)->get_other_ind( i ) );

			// check for potential interaction of o2* of this new residue and any atom in moving residue.
			Vector const & o2star_other = pose.residue( j ).xyz( " O2*" );
			for ( Size n = 1; n <= pose.residue( i ).natoms(); n++ ){
				if ( ( pose.residue( i ).xyz( n ) - o2star_other ).length() < DIST_CUTOFF ) {
					residues_allowed_to_be_packed[ j ] = true;
					break;
				}
			}

			// check for potential interaction of o2* of moving residue and any atom in this new residue
			if (residues_allowed_to_be_packed[ i ]) continue;

			Vector const & o2star_i = pose.residue( i ).xyz( " O2*" );
			for ( Size n = 1; n <= pose.residue( j ).natoms(); n++ ){
				if ( ( pose.residue( j ).xyz( n ) - o2star_i ).length() < DIST_CUTOFF ) {
					residues_allowed_to_be_packed[ i ] = true;
					break;
				}
			}

		}
	}

	utility::vector1< Size > res_list;
	for ( Size n = 1; n <= pose.total_residue(); n++ ) {
		if ( residues_allowed_to_be_packed[ n ] ) {
			res_list.push_back( n );
		}
	}
	if (res_list.size()==0) return 0; //nothing to move!

	Size const o2star_idx_within_res_list = int(  res_list.size() * RG.uniform() ) + 1;
	Size const o2star_num = res_list[ o2star_idx_within_res_list ];
	return o2star_num;
}


//////////////////////////////////////////////////////////////////////////////////////
// Copied from swa_rna_main
utility::vector1< core::Size >
get_fixed_res(core::Size const nres){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::swa::rna;

	utility::vector1< Size > actual_fixed_res_list;
	actual_fixed_res_list.clear();

	utility::vector1< core::Size > const fixed_res_list = option[ fixed_res  ]();
	utility::vector1< core::Size > const minimize_res_list= option[ minimize_res ]();

	if(fixed_res_list.size()!=0 && minimize_res_list.size()!=0 ){
		utility_exit_with_message( "User Cannot specify both  fixed_res and minimize_res!" );
	}


	if( fixed_res_list.size()!=0  ){
		actual_fixed_res_list=fixed_res_list;

	}else if( minimize_res_list.size()!=0){

		for(Size seq_num=1; seq_num<=nres; seq_num++){
			if( Contain_seq_num( seq_num, minimize_res_list) ) continue;
			actual_fixed_res_list.push_back(seq_num);
		}

	}else{ //here I am being a little stringent and require user specify one of these option. Could just return empty list...
		utility_exit_with_message( "User did not specify both fixed res and minimize_res!" );
	}

	return actual_fixed_res_list;
}


//////////////////////////////////////////////////////////////////////////////////////
// Copied from swa_rna_main
utility::vector1< core::Size >
get_input_res(core::Size const nres , std::string const pose_num){


	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::swa::rna;

	utility::vector1< core::Size > input_res_list;
	utility::vector1< core::Size > missing_res_list;

	if(pose_num=="1"){
		input_res_list= option[ input_res ]();
		missing_res_list= option[ missing_res ]();
	}else if(pose_num=="2"){
		input_res_list= option[ input_res2 ]();
		missing_res_list= option[ missing_res2 ]();
	}else{
		utility_exit_with_message( "Invalid pose_num " + pose_num + ", must by either 1 or 2 !" );
	}


	if( input_res_list.size()!=0 && missing_res_list.size()!=0 ){
		utility_exit_with_message( "User Cannot specify both input_res" + pose_num + " and missing_res" + pose_num + "!" );
	}

	utility::vector1< core::Size > actual_input_res_list;
	actual_input_res_list.clear();

	if( input_res_list.size()!=0){
		actual_input_res_list=input_res_list;

	}else if( missing_res_list.size()!=0){

		for(Size seq_num=1; seq_num<=nres; seq_num++){
			if( Contain_seq_num( seq_num, missing_res_list) ) continue;
			actual_input_res_list.push_back(seq_num);
		}

	}else{ //did not specify both input_res and missing_res, return empty list
		std::cout << "user did not specify both input_res" << pose_num << " and missing_res" << pose_num << std::endl;
	}

	return actual_input_res_list;

}

///////////////////////////////////////////////////////////////////////////
// Copied from swa_rna_main
core::scoring::ScoreFunctionOP
create_scorefxn(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;


	std::string score_weight_file;

	Size num_score_weight_file=0;

	if ( option[ basic::options::OptionKeys::score::weights ].user() ) {
		score_weight_file= option[ basic::options::OptionKeys::score::weights ]();
		std::cout << "User passed in score:weight option: " << score_weight_file << std::endl;
		num_score_weight_file++;
	}


	if(num_score_weight_file==0){
		//rna_loop_hires_04092010.wts is same as 5X_linear_quarter_fa_stack_and_adjust_bulge_ss_benchmark.wts
		//change default from single_strand_benchmark to 5X_linear_chainbreak_single_strand_benchmark on May 24, 2010
		//change default to 5X_linear_quarter_fa_stack_and_adjust_bulge_ss_benchmark.wts" on April 9th, 2011
		//score_weight_file="rna_loop_hires_04092010.wts";
		utility_exit_with_message("User to need to pass in score:weights"); //Remove the default weight on Sept 28, 2011 Parin S.
	}

	if(num_score_weight_file>1){
		std::cout << "num_score_weight_file (inputted by user)=" << num_score_weight_file << std::endl;
		utility_exit_with_message("num_score_weight_file>1");
	}

	core::scoring::ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( score_weight_file );


	std::cout << "---------score function weights----------" << std::endl;
	scorefxn->show(std::cout);
	std::cout << "-----------------------------------------" << std::endl;


	return scorefxn;
}


void
setup_copy_DOF_input(protocols::swa::rna::StepWiseRNA_PoseSetupOP & stepwise_rna_pose_setup){

	/////////////////////////////////////////////////////////////////////////////////////////
	// StepWisePoseSetup should create the starting pose.
	// This class might eventually be united with the protein StepWisePoseSetup.
	utility::vector1< std::string > input_tags;
	utility::vector1< std::string > silent_files_in;

	if ( option[ in::file::s ].user() ) {
		// Then any pdbs that need to be read in from disk.
		utility::vector1< std::string > const	pdb_tags_from_disk( option[ in::file::s ]() );
		for ( Size n = 1; n <= pdb_tags_from_disk.size(); n++ ) {
			input_tags.push_back( pdb_tags_from_disk[ n ] );
		}
	}

	if(input_tags.size() > 2 ){
		utility_exit_with_message( "input_tags.size() > 2!!" );
	}

	std::cout << "Input structures for COPY DOF" << std::endl;
	for(Size n=1; n<=input_tags.size(); n++){
		if(n<=silent_files_in.size()){
			std::cout << "silent_file tag= " << input_tags[n] << " silent_file= " << silent_files_in[n] << std::endl;
		}else{
			std::cout << "input_tag= " << input_tags[n] << std::endl;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	stepwise_rna_pose_setup->set_input_tags( input_tags);
	stepwise_rna_pose_setup->set_silent_files_in( silent_files_in);


}


//////////////////////////////////////////////////////////////////////////////////////
// Copied from swa_rna_main -- removed stuff to check silent files.
// probably don't need to set all these options... anyway.
//
protocols::swa::rna::StepWiseRNA_JobParametersOP
setup_rna_job_parameters(){

	using namespace protocols::swa::rna;
	using namespace ObjexxFCL;
	///////////////////////////////
	// Read in sequence.
	if ( !option[ in::file::fasta ].user() ) utility_exit_with_message( "Must supply in::file::fasta!" );
	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
	std::string const full_sequence = fasta_sequence->sequence();
	core::Size const nres=full_sequence.length();

	if ( !option[ sample_res ].user() ) utility_exit_with_message( "Must supply sample_res!" );


	/////////////////////////////////////////////////////

	StepWiseRNA_JobParameters_Setup stepwise_rna_job_parameters_setup( option[ sample_res ](), /*the first element of moving_res_list is the sampling_res*/
																								 										 full_sequence,
																								 										 get_input_res(nres, "1" ),
																								 										 get_input_res(nres, "2" ),
																								 										 option[ cutpoint_open ](),
																								 										 option[ cutpoint_closed ]() );
	stepwise_rna_job_parameters_setup.set_simple_append_map( option[ simple_append_map]() );
	stepwise_rna_job_parameters_setup.set_allow_fixed_res_at_moving_res( option[ allow_fixed_res_at_moving_res ]() ); //Hacky just to get Hermann Duplex working. Need to called before set_fixed_res
	stepwise_rna_job_parameters_setup.set_fixed_res( get_fixed_res(nres) );
	stepwise_rna_job_parameters_setup.set_terminal_res( option[ terminal_res ]() );
	stepwise_rna_job_parameters_setup.set_rmsd_res_list( option[ rmsd_res ]() );
	stepwise_rna_job_parameters_setup.set_jump_point_pair_list( option[ jump_point_pairs ]() ); //Important!: Need to be called after set_fixed_res
	stepwise_rna_job_parameters_setup.set_alignment_res( option[ alignment_res ]() );
	stepwise_rna_job_parameters_setup.set_filter_user_alignment_res( option[ filter_user_alignment_res ]() );
	stepwise_rna_job_parameters_setup.set_native_alignment_res( option[ native_alignment_res ]() );

	stepwise_rna_job_parameters_setup.set_global_sample_res_list( option[ global_sample_res_list ]() ); //March 20, 2011

	stepwise_rna_job_parameters_setup.set_force_syn_chi_res_list( option[ force_syn_chi_res_list]() ); //April 29, 2011
	stepwise_rna_job_parameters_setup.set_force_north_ribose_list( option[ force_north_ribose_list ]() ); //April 29, 2011
	stepwise_rna_job_parameters_setup.set_force_south_ribose_list( option[ force_south_ribose_list ]() ); //April 29, 2011
	stepwise_rna_job_parameters_setup.set_protonated_H1_adenosine_list( option[ protonated_H1_adenosine_list ]() ); //May 02, 2011

	stepwise_rna_job_parameters_setup.set_allow_chain_boundary_jump_partner_right_at_fixed_BP( option[ allow_chain_boundary_jump_partner_right_at_fixed_BP ]() ); //Hacky just to get Square RNA working.

	stepwise_rna_job_parameters_setup.set_output_extra_RMSDs( option[ output_extra_RMSDs ]() );
	stepwise_rna_job_parameters_setup.set_add_virt_res_as_root( option[ add_virt_root]() );


	stepwise_rna_job_parameters_setup.set_skip_complicated_stuff( true ); // new by Rhiju.

	stepwise_rna_job_parameters_setup.apply();

	return stepwise_rna_job_parameters_setup.job_parameters();

}



//////////////////////////////////////////////////////////////////////////////////////////
// Copied from swa_rna_main
protocols::swa::rna::StepWiseRNA_PoseSetupOP
setup_pose_setup_class(protocols::swa::rna::StepWiseRNA_JobParametersOP & job_parameters, bool const copy_DOF=true){

  using namespace core::pose;
  using namespace core::chemical;
  using namespace core::kinematics;
  using namespace core::scoring;
	using namespace protocols::swa::rna;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	// Read in native_pose.
	PoseOP native_pose;
	if (option[ in::file::native ].user() ) {
		native_pose = PoseOP( new Pose );
		import_pose::pose_from_pdb( *native_pose, *rsd_set, option[ in::file::native ]() );
		std::cout << "native_pose->fold_tree(): " << native_pose->fold_tree();
		std::cout << "native_pose->annotated_sequence(true): " << native_pose->annotated_sequence( true ) << std::endl;
		protocols::rna::make_phosphate_nomenclature_matches_mini( *native_pose);
	}

	StepWiseRNA_PoseSetupOP stepwise_rna_pose_setup = new StepWiseRNA_PoseSetup( job_parameters);
	stepwise_rna_pose_setup->set_copy_DOF(copy_DOF);

	if(copy_DOF==true){
		setup_copy_DOF_input(stepwise_rna_pose_setup);
	}


	stepwise_rna_pose_setup->set_virtual_res( option[ virtual_res ]() );
	stepwise_rna_pose_setup->set_bulge_res( option[ bulge_res ]() );
	stepwise_rna_pose_setup->set_native_pose( native_pose );
	stepwise_rna_pose_setup->set_native_virtual_res( option[ native_virtual_res]() );
	stepwise_rna_pose_setup->set_rebuild_bulge_mode( option[rebuild_bulge_mode]() );
	stepwise_rna_pose_setup->set_output_pdb( option[ output_pdb ]() );
	stepwise_rna_pose_setup->set_apply_virtual_res_variant_at_dinucleotide( false );
	stepwise_rna_pose_setup->set_align_to_native( true );

	return stepwise_rna_pose_setup;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Following assumed that pose is already properly aligned to native pose!
void
output_silent_file( pose::Pose const & pose, Size const n_accept, Size const count,
										utility::vector1< Size > working_res_list,
										std::map< Size, Size > sub_to_full,
										pose::Pose const & native_pose,
										core::io::silent::SilentFileDataOP silent_file_data,
										std::string const & silent_file ){

  using namespace core::io::silent;
  using namespace core::conformation;

	// useful to keep track of what's a working residue and what's not.
	utility::vector1< bool > is_working_res;
	for ( Size i = 1; i <= pose.total_residue(); i++ ) is_working_res.push_back( false );
	for ( Size n = 1; n <= working_res_list.size(); n++ ) is_working_res[ working_res_list[n] ] = true;

	std::string const tag = "S_"+lead_zero_string_of(n_accept,6);
	BinaryRNASilentStruct s( pose, tag );

	// could initialize this once somewhere else.
	utility::vector1< std::string > next_suite_atoms;
	next_suite_atoms.push_back(" P  ");
	next_suite_atoms.push_back(" O1P");
	next_suite_atoms.push_back(" O2P");
	next_suite_atoms.push_back(" O5*");

	Real dev( 0.0 );
	Real rmsd( 0.0 );
	Size natoms( 0 );

	for ( Size n = 1; n <= working_res_list.size(); n++ ){
		Size const i      = working_res_list[ n ];
		Size const i_full = sub_to_full[ i ];

		Residue const & rsd        = pose.residue( i );
		Residue const & rsd_native = native_pose.residue( i_full );
		runtime_assert( rsd.aa() == rsd_native.aa() );

		for ( Size j = 1; j <= rsd.nheavyatoms(); j++ ){
			if ( rsd.is_virtual( j ) ) continue;  // note virtual phosphates on 5'-ends of loops.
			std::string atom_name = rsd.atom_name( j );
			//			std::cout << "RMSD: " << i << atom_name << std::endl;
			if ( !rsd_native.has( atom_name ) ) continue;
			Size const j_full = rsd_native.atom_index( atom_name );
			dev += ( rsd_native.xyz( j_full ) - rsd.xyz( j ) ).length_squared();
			natoms++;
		}

		// also add in atoms in next suite, if relevant (and won't be covered later in rmsd calc.)
		if ( i < pose.total_residue() && !is_working_res[ i+1 ] &&  !pose.fold_tree().is_cutpoint(i) ){
			Size const i_next      = i+1;
			Size const i_next_full = sub_to_full[ i+1 ];
			runtime_assert( i_next_full == i_full + 1 ); //better be a connection in both the pose & native pose!

			Residue const & rsd_next        = pose.residue( i_next );
			Residue const & rsd_next_native = native_pose.residue( i_next_full );
			runtime_assert( rsd_next.aa() == rsd_next_native.aa() );

			for (Size k = 1; k <= next_suite_atoms.size(); k++ ){
				std::string atom_name = next_suite_atoms[ k ];
				//std::cout << "RMSD: " << i+1 << atom_name << std::endl;
				runtime_assert( rsd_next.has( atom_name ) );
				runtime_assert( rsd_next_native.has( atom_name ) );
				Size const j = rsd_next.atom_index( atom_name );
				Size const j_full = rsd_next_native.atom_index( atom_name );
				dev += ( rsd_next_native.xyz( j_full ) - rsd_next.xyz( j ) ).length_squared();
				natoms++;
			}

		}

	}
	if ( natoms > 0 ) rmsd = std::sqrt( dev / static_cast<Real>( natoms ) );

	s.add_energy( "rms",rmsd );
	s.add_string_value( "count", ObjexxFCL::fmt::I(9,count) );

	std::string built_res = "";
	for ( Size n = 1; n <= working_res_list.size(); n++ ) built_res += string_of(working_res_list[n]);
	s.add_string_value( "built_res", built_res);

	silent_file_data->write_silent_struct( s, silent_file, true /*score_only*/ );
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// can go into a namespace later...
enum MovingResidueCase { NONE=0, CHAIN_TERMINUS_5PRIME, CHAIN_TERMINUS_3PRIME, INTERNAL, FLOATING_BASE };

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
random_torsion_move( pose::Pose & pose,
										 utility::vector1< Size > const & moving_res_list,
										 std::string & move_type,
										 Real const & sample_range ){

	using namespace ObjexxFCL::fmt;

	Size const random_idx = int( RG.uniform() * moving_res_list.size() ) + 1;
	Size const i = moving_res_list[ random_idx ];

	Size const & nres( pose.total_residue() );
	kinematics::FoldTree const & fold_tree( pose.fold_tree() );

	MovingResidueCase moving_residue_case;
	if ( i == nres || fold_tree.is_cutpoint( i ) ){ // could be a 5' chain terminus
		if ( i == 1 || fold_tree.is_cutpoint( i-1 ) ){
			moving_residue_case = FLOATING_BASE; // don't know how to handle this yet.
		} else {
			moving_residue_case = CHAIN_TERMINUS_3PRIME;
		}
	} else if ( i == 1 || fold_tree.is_cutpoint( i-1 ) ){
		moving_residue_case = CHAIN_TERMINUS_5PRIME;
	} else {
		moving_residue_case = INTERNAL;
	}


	if ( moving_residue_case == CHAIN_TERMINUS_3PRIME  || moving_residue_case == CHAIN_TERMINUS_5PRIME ){

		// an edge residue -- change both its nucleoside & suite -- can go crazy.
		Size const nucleoside_num = i;
		sample_near_nucleoside_torsion( pose, nucleoside_num, sample_range);

		Size suite_num( 0 );
		if ( moving_residue_case == CHAIN_TERMINUS_3PRIME ) suite_num = i-1;
		else suite_num = i;

		if ( RG.uniform() < 0.5) {
			sample_near_suite_torsion( pose, suite_num, sample_range);
			move_type += "-nuc-suite";
		} else {
			crankshaft_alpha_gamma( pose, suite_num, sample_range);
			move_type += "-nuc-crank";
		}

	} else {
		runtime_assert( moving_residue_case == INTERNAL ); // cannot handle floating base yet.

		// don't do anything super-crazy -- either do previous suite, current nucleoside, or next suite.
		Real const random_number = RG.uniform();

		if ( random_number < 0.6 ){

			Size suite_num( 0 );
			if ( RG.uniform() < 0.5) {
				suite_num= i-1;
			} else {
				suite_num= i;
			}

			if ( RG.uniform() < 0.5) {
				sample_near_suite_torsion( pose, suite_num, sample_range);
				//move_type += "-suite" + string_of(suite_num);
				move_type += "-suite";
			} else {
				crankshaft_alpha_gamma( pose, suite_num, sample_range);
				//move_type += "-suite" + string_of(suite_num);
				move_type += "-crank";
			}
		} else {
			Size const nucleoside_num = i;
			sample_near_nucleoside_torsion( pose, nucleoside_num, sample_range);
			//move_type += "-nuc" + string_of(nucleoside_num);
			move_type += "-nuc";
		}

	}

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1<Size>
reorder_after_delete( utility::vector1<Size> & moving_res_list,
											Size const & res_to_delete ){

	utility::vector1< Size > moving_res_list_new;

	for (Size i = 1; i <= moving_res_list.size(); i++ ){
		Size const n = moving_res_list[ i ];
		if ( n < res_to_delete ) moving_res_list_new.push_back( n );
		else if (n > res_to_delete ) moving_res_list_new.push_back( n-1 );
	}

	return moving_res_list_new;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::map<Size,Size>
reorder_after_delete( std::map<Size,Size> & sub_to_full,
											Size const & res_to_delete ){

	std::map< Size, Size > sub_to_full_new;

	for ( std::map< Size, Size >::const_iterator it = sub_to_full.begin(); it != sub_to_full.end(); ++it ) {
		Size const n = it->first;
		Size const m = it->second;
		if ( n < res_to_delete ) sub_to_full_new[ n ] = m;
		else if ( n > res_to_delete ) sub_to_full_new[ n-1 ] = m;
	}

	return sub_to_full_new;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1<Size>
reorder_after_insert( utility::vector1<Size> & moving_res_list,
											 Size const & res_to_add ){

	utility::vector1< Size > moving_res_list_new;

	for (Size i = 1; i <= moving_res_list.size(); i++ ){
		Size const n = moving_res_list[ i ];
		if ( n < res_to_add ) moving_res_list_new.push_back( n );
	}
	moving_res_list_new.push_back( res_to_add );
	for (Size i = 1; i <= moving_res_list.size(); i++ ){
		Size const n = moving_res_list[ i ];
		if ( n >= res_to_add ) moving_res_list_new.push_back( n+1 );
	}

	return moving_res_list_new;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
std::map<Size,Size>
reorder_after_insert( std::map<Size,Size> & sub_to_full,
											Size const & res_to_add ){

	std::map< Size, Size > sub_to_full_new;

	for ( std::map< Size, Size >::const_iterator it = sub_to_full.begin(); it != sub_to_full.end(); ++it ) {
		Size const n = it->first;
		Size const m = it->second;
		if ( n < res_to_add )  sub_to_full_new[ n ] = m;
		if ( n >= res_to_add ) sub_to_full_new[ n+1 ] = m;
	}
	sub_to_full_new[ res_to_add ] = sub_to_full[ res_to_add ]-1;

	return  sub_to_full_new;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
swa_rna_sample()
{

  using namespace core::pose;
  using namespace core::chemical;
  using namespace core::kinematics;
  using namespace core::scoring;
  using namespace core::io::silent;
	using namespace protocols::swa::rna;
	using namespace protocols::moves;

	Output_title_text("Enter swa_rna_sample()");

	clock_t const time_start( clock() );

	// pose_setup stuff copied from swa_rna_main
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	std::cout << "Total time to setup ResidueTypeSet: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;

	core::scoring::ScoreFunctionOP scorefxn=create_scorefxn();

	///////////////////////////////
	StepWiseRNA_JobParametersOP	job_parameters = setup_rna_job_parameters();
	StepWiseRNA_JobParametersCOP job_parameters_COP( job_parameters );

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	StepWiseRNA_PoseSetupOP stepwise_rna_pose_setup = setup_pose_setup_class(job_parameters);

  Pose pose;
	stepwise_rna_pose_setup->apply( pose );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 800, 800 );

	stepwise_rna_pose_setup->setup_native_pose( pose ); //NEED pose to align native_pose to pose.

	pose.dump_pdb( "START.pdb" );
	job_parameters->working_native_pose()->dump_pdb( "working_native.pdb" );
	job_parameters->working_native_pose()->dump_pdb( "working_native.pdb" );

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Try Fang's Monte Carlo machinery
	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	Size const num_cycle = option[ n_sample ]();
	Real const sample_range_small = option[ stddev_small ]();
	Real const sample_range_large = option[ stddev_large ]();

	utility::vector1< Size > moving_res_list = job_parameters->working_moving_res_list();

	// is this in the right order or what? Total hack for now.
	if ( moving_res_list.size() > 1 && moving_res_list[2] < moving_res_list[1]){
		utility::vector1< Size > moving_res_list_new;
		for ( Size n = moving_res_list.size(); n >= 1; n-- ) moving_res_list_new.push_back( moving_res_list[ n ] );
		moving_res_list = moving_res_list_new;
	}
	std::cout << "MOVING_RES ";
	for (Size i = 1; i <= moving_res_list.size(); i++ ) std::cout << ' ' <<  moving_res_list[i];
	std::cout << std::endl;

	std::map< Size, Size > sub_to_full = job_parameters->sub_to_full();
	std::string full_sequence = job_parameters->full_sequence();
	pose::Pose const & native_pose = *( stepwise_rna_pose_setup->get_native_pose() );

	std::string const silent_file = option[ out::file::silent ]();
	SilentFileDataOP silent_file_data = new SilentFileData;

	( *scorefxn )( pose ); //score it for first silent output
	output_silent_file( pose, 0, 0, moving_res_list, sub_to_full, native_pose, silent_file_data, silent_file );
	Real const output_score_cutoff_ = option[ output_score_cutoff ]();
	bool const do_add_delete_ = option[ do_add_delete ]();
	bool const presample_added_residue_ = option[ presample_added_residue ]();
	Size const internal_cycles_ = option[ presample_internal_cycles ](); // totally arbitrary


	std::string move_type( "" );
	if ( ! option[ skip_randomize ]() ){
		std::cout << "randomizing... ";
		for (Size count = 1; count <= 1000; count++) {
			move_type = "lrg";
			random_torsion_move( pose, moving_res_list, move_type, sample_range_large );
		}
		std::cout << " done. " << std::endl;
	}

	MonteCarloOP monte_carlo_ = new MonteCarlo( pose, *scorefxn, option[ kT ]() );

	bool accepted( true );
	Size n_accept( 0 );
	Size o2star_res( 0 );
	std::map< Size, Size > sub_to_full_new;
	utility::vector1< Size > moving_res_list_new;

	for (Size count = 1; count <= num_cycle; count++) {

		if ( count % 1000 == 0 ) {
			std::cout << "On " << count << " of " << num_cycle << " trials." << std::endl;
		}

		Real const random_number = RG.uniform();

		move_type = "";
		moving_res_list_new = moving_res_list;
		sub_to_full_new = sub_to_full;

		if ( random_number < 0.01 && do_add_delete_ ) {

			Real const random_number2 = RG.uniform();
			if ( random_number2 < 0.5 ) {
				///////////////////////////////////
				// Deletion
				///////////////////////////////////
				// try to delete a sampled residue
				if ( moving_res_list.size() > 1 ){
					move_type = "delete";

					std::cout << "Before delete: " << (*scorefxn)( pose ) << std::endl;

					bool const forward_build_ = option[ forward_build ](); // later, both forward and backward will be possible -- and determined automatically
					if ( forward_build_ ){

						Size const res_to_delete = moving_res_list[ moving_res_list.size() ];
						pose.delete_polymer_residue( res_to_delete ); // only works for fragment built backward off 5' fragment
						moving_res_list_new = reorder_after_delete( moving_res_list, res_to_delete );
						sub_to_full_new = reorder_after_delete( sub_to_full, res_to_delete );


					} else {
						Size const res_to_delete = moving_res_list[1];
						pose.delete_polymer_residue( res_to_delete ); // only works for fragment built backward off 5' fragment
						moving_res_list_new = reorder_after_delete( moving_res_list, res_to_delete );
						sub_to_full_new = reorder_after_delete( sub_to_full, res_to_delete );
						pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", res_to_delete );

					}

					std::cout << "After delete: " << (*scorefxn)( pose ) << std::endl << std::endl;

				}

			} else {
				///////////////////////////////////
				// Addition
				///////////////////////////////////
				// try to add a residue that is supposed to be sampled.

				move_type = "add";
				std::cout << "Before add: " << (*scorefxn)( pose ) << std::endl;

				Size suite_num( 0 ), nucleoside_num( 0 ); // will record which new dofs added.

				bool did_addition( false );
				//pose.dump_pdb( "before_add.pdb" );

				bool const forward_build_ = option[ forward_build ](); // later, both forward and backward will be possible -- and determined automatically
				if ( forward_build_ ){

					Size const res_to_build_off = moving_res_list[ moving_res_list.size() ]; // for now, just build off 3' fragment.

					if ( res_to_build_off < pose.total_residue() &&
							 sub_to_full[ res_to_build_off ] < sub_to_full[ res_to_build_off+1 ] -1 ){  // need to fix for general case [e.g. cutpoint]
						Size const res_to_add = res_to_build_off + 1;

						char newrestype = full_sequence[ (sub_to_full[ res_to_build_off ] + 1) - 1 ];
						//std::cout << "I want to add: " << newrestype << std::endl;

						chemical::AA my_aa = chemical::aa_from_oneletter_code( newrestype );
						chemical::ResidueTypeCAPs const & rsd_type_list( rsd_set->aa_map( my_aa ) );
						// iterate over rsd_types, pick one.
						chemical::ResidueType const & rsd_type = *rsd_type_list[1];
						core::conformation::ResidueOP new_rsd = conformation::ResidueFactory::create_residue( rsd_type );

						pose.append_polymer_residue_after_seqpos( *new_rsd, res_to_build_off, true /*build ideal geometry*/ );

						suite_num = res_to_add-1;
						nucleoside_num = res_to_add;

						did_addition = true;
					}
				} else {
					Size const res_to_add = moving_res_list[1]; // for now, just build off 5' fragment.
					if ( sub_to_full[ res_to_add ] > 1 ){  // again, need to fix for general case.

						char newrestype = full_sequence[ (sub_to_full[ res_to_add ] - 1) - 1 ];
						//std::cout << "I want to add: " << newrestype << std::endl;

						chemical::AA my_aa = chemical::aa_from_oneletter_code( newrestype );
						chemical::ResidueTypeCAPs const & rsd_type_list( rsd_set->aa_map( my_aa ) );
						// iterate over rsd_types, pick one.
						chemical::ResidueType const & rsd_type = *rsd_type_list[1];
						core::conformation::ResidueOP new_rsd = conformation::ResidueFactory::create_residue( rsd_type );

						pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", res_to_add ); // got to be safe.

						pose.prepend_polymer_residue_before_seqpos( *new_rsd, res_to_add, true /*build ideal geometry*/ );
						pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", res_to_add );

						// initialize with a random torsion... ( how about an A-form + perturbation ... or go to a 'reasonable' rotamer)
						suite_num = res_to_add;
						nucleoside_num = res_to_add;

						did_addition = true;
					}
				}


				if ( did_addition ){
					moving_res_list_new = reorder_after_insert( moving_res_list, nucleoside_num );
					sub_to_full_new = reorder_after_insert( sub_to_full, nucleoside_num );

					apply_suite_torsion_Aform( pose, suite_num );
					apply_nucleoside_torsion_Aform( pose, nucleoside_num );
					sample_near_suite_torsion( pose, suite_num, sample_range_large);
					sample_near_nucleoside_torsion( pose, nucleoside_num, sample_range_large);

					///////////////////////////////////
					// Presampling added residue
					///////////////////////////////////
					if ( presample_added_residue_ ){
						std::cout << "presampling added residue! " << nucleoside_num << std::endl;
						MonteCarloOP monte_carlo_internal = new MonteCarlo( pose, *scorefxn, option[ kT ]() );

						for ( Size count_internal = 1; count_internal <= internal_cycles_; count_internal++ ){

							sample_near_suite_torsion( pose, suite_num, sample_range_large);
							sample_near_nucleoside_torsion( pose, nucleoside_num, sample_range_large);
							monte_carlo_internal->boltzmann( pose, move_type );
							sample_near_suite_torsion( pose, suite_num, sample_range_small);
							sample_near_nucleoside_torsion( pose, nucleoside_num, sample_range_small);
							monte_carlo_internal->boltzmann( pose, move_type );
							//std::cout << "During presampling: " << (*scorefxn)( pose );
						}
					}
					//std::cout << pose.annotated_sequence() << std::endl;
					std::cout << "After add: " << (*scorefxn)( pose ) << std::endl << std::endl;
					//					pose.dump_pdb( "after_add.pdb" );
				} else {
					move_type = ""; // no move!
				}
			}

		} else if ( random_number  < 0.8 ){

			Real const random_number2 = RG.uniform();
			if ( random_number2  < 0.5 ){
				move_type = "sml";
				random_torsion_move( pose, moving_res_list, move_type, sample_range_small );
			} else {
				move_type = "lrg";
				random_torsion_move( pose, moving_res_list, move_type, sample_range_large );
			}
		} else{
			// perhaps should also move 2'-OH torsions?
			if ( option[ sample_all_o2star ]() ){
				o2star_res = get_random_o2star_residue( pose );
			} else {
				// warning -- following might lead to weird 'hysteresis' effects since it picks
				// o2star to sample based on what's near moving residue.
				o2star_res = get_random_o2star_residue_near_moving_residue( pose, moving_res_list );
			}
			if ( o2star_res > 0 ) {
				Real const random_number2 = RG.uniform();
				if ( random_number2  < 0.5 ){
					//move_type = "sml-o2star"+string_of( o2star_res);
					move_type = "sml-o2star";
					sample_near_o2star_torsion( pose, o2star_res, sample_range_small);
				} else {
					//move_type = "lrg-o2star"+string_of( o2star_res );
					move_type = "lrg-o2star";
					sample_near_o2star_torsion( pose, o2star_res, sample_range_large);
				}
			}
		}

		if ( move_type.size() == 0 ){
			count--; continue; // slight hack -- try monte carlo cycle again. dangerous -- might cause infinite loop
		}

		accepted = monte_carlo_->boltzmann( pose, move_type );

		if ( accepted ) { //slight pain in the ass. Maybe we should keep these cached in the pose somewhere
			moving_res_list = moving_res_list_new;
			sub_to_full     = sub_to_full_new;
		}
		//std::cout << "Score: " << (*scorefxn)( pose ) << std::endl;
		Real const current_score = (*scorefxn)( pose );

		// Need to fix following to not be dependent on job_parameters
		if (accepted && current_score < output_score_cutoff_ ) output_silent_file( pose, n_accept, count, moving_res_list, sub_to_full, native_pose, silent_file_data, silent_file );

	}

	output_silent_file( pose, n_accept, num_cycle, moving_res_list, sub_to_full, native_pose, silent_file_data, silent_file );

	pose.dump_pdb( "FINAL.pdb" );
	monte_carlo_->show_counters();
	monte_carlo_->recover_low( pose );


	//output_silent_file( pose, n_accept, num_cycle+1, moving_res_list, sub_to_full, native_pose, silent_file_data, silent_file );
	if ( option[ out::file::o].user() ) { // makes life easier on cluter.
		pose.dump_pdb( option[ out::file::o ]() );
	} else {
		pose.dump_pdb( "LOW.pdb" );
	}
	std::cout << "Total time for monte carlo: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;


}



///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	clock_t const my_main_time_start( clock() );

  using namespace basic::options;

	std::string algorithm_input = option[algorithm];

	swa_rna_sample();

	protocols::viewer::clear_conformation_viewers();

	std::cout << "Total time took to run algorithm (" << algorithm_input << "): " << static_cast<Real>( clock() - my_main_time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;

	std::cout << "JOB_SUCCESSFULLY_COMPLETED" << std::endl;

  exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
  using namespace basic::options;

	utility::vector1< Size > blank_size_vector;
	utility::vector1< std::string > blank_string_vector;

	//////////////General/////////////////////////////
	NEW_OPT( graphic, "Turn graphic on/off", true);
	NEW_OPT( Real_parameter_one, "free_variable for testing purposes ", 0.0);
	NEW_OPT( distinguish_pucker, "distinguish pucker when cluster:both in sampler and clusterer", true);
	NEW_OPT( output_pdb, "output_pdb: If true, then will dump the pose into a PDB file at different stages of the stepwise assembly process.", false); //Sept 24, 2011
	NEW_OPT( algorithm, "Specify algorithm to execute", "");

	//////////////Job_Parameters///////////
	NEW_OPT( sample_res, "residues to build, the first element is the actual sample res while the other are the bulge residues", blank_size_vector );
	NEW_OPT( input_res, "Residues already present in starting pose_1", blank_size_vector );
	NEW_OPT( input_res2, "Residues already present in starting  pose_2", blank_size_vector );
	NEW_OPT( missing_res, "Residues missing in starting pose_1, alternative to input_res", blank_size_vector );
	NEW_OPT( missing_res2, "Residues missing in starting pose_2, alternative to input_res2", blank_size_vector );
	NEW_OPT( rmsd_res, "residues that will be use to calculate rmsd (for clustering as well as RMSD to native_pdb if specified)", blank_size_vector );
	NEW_OPT( alignment_res , "align_res_list", blank_string_vector );
	NEW_OPT( global_sample_res_list, "A list of all the nucleotide to be build/sample over the entire dag.", blank_size_vector); //March 20, 2011

	NEW_OPT( cutpoint_open, "optional: chainbreak in full sequence", blank_size_vector );
	NEW_OPT( cutpoint_closed, "optional: cutpoint at which to apply chain closure", 0 );
	NEW_OPT( jump_point_pairs , "optional: extra jump_points specified by the user for setting up the fold_tree ", blank_string_vector );

	NEW_OPT( native_virtual_res , " optional: native_virtual_res ", blank_size_vector );
	NEW_OPT( native_alignment_res , "optional: native_alignment_res ", blank_size_vector );
	NEW_OPT( fixed_res, "optional: residues to be held fixed in minimizer", blank_size_vector );
	NEW_OPT( minimize_res, "optional: residues to be minimize in minimizer, alternative to fixed_res", blank_size_vector );
	NEW_OPT( virtual_res, "optional: residues to be made virtual", blank_size_vector );
	NEW_OPT( terminal_res, "optional: residues that are not allowed to stack during sampling", blank_size_vector );
	NEW_OPT( bulge_res, "optional: residues to be turned into a bulge variant", blank_size_vector );
	NEW_OPT( force_syn_chi_res_list, "optional: sample only syn chi for the res in sampler.", blank_size_vector); //April 29, 2011
	NEW_OPT( force_north_ribose_list, "optional: sample only north ribose for the res in sampler.", blank_size_vector); //April 29, 2011
	NEW_OPT( force_south_ribose_list, "optional: sample only south ribose for the res in sampler.", blank_size_vector); //April 29, 2011
	NEW_OPT( protonated_H1_adenosine_list, "optional: protonate_H1_adenosine_list", blank_size_vector); //May 02, 2011

	//////////////Pose setup///////
	NEW_OPT( job_queue_ID, " swa_rna_sample()/combine_long_loop mode: Specify the tag pair in filter_output_filename to be read in and imported (start from 0!)" , 0);

	///////////////Sampler////////////
	NEW_OPT( rebuild_bulge_mode, "rebuild_bulge_mode", false);
	NEW_OPT( floating_base , " floating_base ", false ); //DO NOT CHANGE TO TRUE, since single-nucleotide sampling need this to be false! April 9th, 2011

	//////////////CombineLongLoopFilterer/////////////
	NEW_OPT( filter_output_filename, "CombineLongLoopFilterer: filter_output_filename", "filter_struct.txt"); //Sept 12, 2010
	NEW_OPT( filter_for_previous_contact, "CombineLongLoopFilterer: filter_for_previous_contact", false); //Sept 12, 2010
	NEW_OPT( filter_for_previous_clash, "CombineLongLoopFilterer: filter_for_previous_clash", false); //Sept 12, 2010
	NEW_OPT( combine_helical_silent_file, "CombineLongLoopFilterer: combine_helical_silent_file", false); //Nov 27, 2010

	//////////////post_rebuild_bulge_assembly//////
	NEW_OPT( start_silent, "start_silent", ""); //Oct 22, 2011
	NEW_OPT( start_tag, "start_tag", ""); //Oct 22, 2011

	///////The options below are for testing purposes. Please do not make any changes without first consulting/////////////
	///////Parin Sripakdeevong (sripakpa@stanford.edu) or Rhiju Das (rhiju@stanford.edu) //////////////////////////////////
	//////////////General/////////////////////////////
	NEW_OPT( VERBOSE, "VERBOSE", false );
	NEW_OPT( parin_favorite_output , " parin_favorite_output ", true ); //Change to true on Oct 10, 2010
	NEW_OPT( integration_test , " integration_test ", false ); //March 16, 2012


	//////////////Job_Parameters///////////
	NEW_OPT( filter_user_alignment_res, " filter_user_alignment_res ", true ); //General want this to be true except for special cases! June 13, 2011
	NEW_OPT( simple_append_map , "simple_append_map", false);
	NEW_OPT( add_virt_root, "add_virt_root", false); //For Fang's electron density code.
	NEW_OPT( allow_chain_boundary_jump_partner_right_at_fixed_BP, "mainly just to get Hermann nano-square RNA modeling to work", false);
	NEW_OPT( allow_fixed_res_at_moving_res, "mainly just to get Hermann Duplex modeling to work", false); //Nov 15, 2010
	NEW_OPT( simple_full_length_job_params, "simple_full_length_job_params", false); //Oct 31, 2011
	NEW_OPT( output_extra_RMSDs, "output_extra_RMSDs", false); //March 16, 2012

	///////////////Sampler////////////
	NEW_OPT( sampler_cluster_rmsd, " Clustering rmsd of conformations in the sampler", 0.5); //DO NOT CHANGE THIS!
	NEW_OPT( skip_sampling, "no sampling step in rna_swa residue sampling", false );
	NEW_OPT( do_not_sample_multiple_virtual_sugar, " Samplerer: do_not_sample_multiple_virtual_sugar " , false);
	NEW_OPT( sample_ONLY_multiple_virtual_sugar, " Samplerer: sample_ONLY_multiple_virtual_sugar " , false);
	NEW_OPT( filterer_undercount_ribose_rotamers, "Undercount all ribose_rotamers as 1 count", false); //July 29, 2011
	NEW_OPT( exclude_alpha_beta_gamma_sampling, "Speed up the debug eplison south sugar mode", false);
	NEW_OPT( debug_eplison_south_sugar_mode, "Check why when eplison is roughly -160 and pucker is south, energy is not favorable", false);
	NEW_OPT( sampler_extra_anti_chi_rotamer, "Samplerer: extra_anti_chi_rotamer", false);
	NEW_OPT( sampler_extra_syn_chi_rotamer, "Samplerer: extra_syn_chi_rotamer", false);
	NEW_OPT( sampler_extra_beta_rotamer, "Samplerer: extra_beta_rotamer", false);
	NEW_OPT( sampler_extra_epsilon_rotamer, "Samplerer: extra_epsilon_rotamer", true); //Change this to true on April 9, 2011
	NEW_OPT( sample_both_sugar_base_rotamer, "Samplerer: Super hacky for SQAURE_RNA", false);
	NEW_OPT( reinitialize_CCD_torsions, "Samplerer: reinitialize_CCD_torsions: Reinitialize_CCD_torsion to zero before every CCD chain closure", false);
	NEW_OPT( PBP_clustering_at_chain_closure, "Samplerer: PBP_clustering_at_chain_closure", false);
	NEW_OPT( finer_sampling_at_chain_closure, "Samplerer: finer_sampling_at_chain_closure", false); //Jun 9, 2010
	NEW_OPT( sampler_include_torsion_value_in_tag, "Samplerer:include_torsion_value_in_tag", true);
	NEW_OPT( include_syn_chi, "include_syn_chi", true); //Change to true on Oct 10, 2010
	NEW_OPT( sampler_allow_syn_pyrimidine, "sampler_allow_syn_pyrimidine", false); //Nov 15, 2010
	NEW_OPT( fast, "quick runthrough for debugging", false );
	NEW_OPT( medium_fast, "quick runthrough for debugging (keep more poses and not as fast as fast option)", false );
	NEW_OPT( centroid_screen, "centroid_screen", true);
	NEW_OPT( allow_base_pair_only_centroid_screen, "allow_base_pair_only_centroid_screen", false); //This only effect floating base sampling + dinucleotide.. deprecate option
	NEW_OPT( sampler_perform_o2star_pack, "perform O2* hydrogen packing inside StepWiseRNA_ResidueSampler", true );
	NEW_OPT( allow_bulge_at_chainbreak, "Allow sampler to replace chainbreak res with virtual_rna_variant if it looks have bad fa_atr score.", true );

	NEW_OPT( add_lead_zero_to_tag, "Add lead zero to clusterer output tag ", false);
	NEW_OPT( n_sample, "Sample number for Random sampling", 0 );
	NEW_OPT( stddev_small, "Sampling standard deviation in degree", 5.0 );
	NEW_OPT( stddev_large, "Sampling standard deviation in degree", 40.0 );
	NEW_OPT( output_score_cutoff, "Score cutoff for output to disk", 0.0 );
	NEW_OPT( kT, "kT of simulation in RU", 2.0 );
	NEW_OPT( skip_randomize, "do not randomize...", false );
	NEW_OPT( sample_all_o2star, "do not focus o2star sampling at residue of interest...", false );
	NEW_OPT( do_add_delete, "try add & delete moves...", false );
	NEW_OPT( presample_added_residue, "when adding a residue, do a little monte carlo to try to get it in place", false );
	NEW_OPT( presample_internal_cycles, "when adding a residue, number of monte carlo cycles", 100 );
	NEW_OPT( forward_build, "TEMPORARY - REMOVE THIS SOON!", false );



	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  ////////////////////////////////////////////////////////////////////////////
  // setup
  ////////////////////////////////////////////////////////////////////////////
  core::init(argc, argv);


  ////////////////////////////////////////////////////////////////////////////
  // end of setup
  ////////////////////////////////////////////////////////////////////////////

  protocols::viewer::viewer_main( my_main );

}



