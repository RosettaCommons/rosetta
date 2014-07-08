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
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/stepwise/sampling/util.hh>
#include <protocols/stepwise/StepWiseLegacyClusterer.hh>
#include <protocols/stepwise/sample_generators/StepWisePoseSampleGenerator.hh>
#include <protocols/stepwise/sample_generators/StepWisePoseCombineSampleGenerator.hh>
#include <protocols/stepwise/RigidBodySampler.hh>
#include <protocols/stepwise/InputStreamWithResidueInfo.hh>
//#include <protocols/stepwise/sampling/rna/rigid_body_settings.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_BaseSugarRotamer.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/sampling/rna/util.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_RotamerGeneratorWrapper.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_RotamerGeneratorWrapper.fwd.hh>
#include <protocols/stepwise/sampling/rna/RNA_AnalyticLoopCloser.hh>
#include <protocols/stepwise/sampling/rna/RNA_LoopCloseSampler.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/file/gzip_util.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

//#include <armadillo>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

#include <core/util/Tracer.hh>
using core::util::T;
#include <time.h>

// option key includes

#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/cluster.OptionKeys.gen.hh>

using namespace core;
using namespace protocols;
using namespace options::OptionKeys;
using ObjexxFCL::format::A;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using io::pdb::dump_pdb;
using utility::vector1;
using utility::tools::make_vector1;
typedef  numeric::xyzMatrix< Real > Matrix;

OPT_KEY( Real, xyz_sample )
OPT_KEY( Real, box_radius )
OPT_KEY( Integer, bin_sample )
OPT_KEY( IntegerVector, minimize_sidechain_res )
OPT_KEY( Boolean, minimize_jump )
OPT_KEY( Boolean, output_all )
OPT_KEY( Boolean, split_silent_files )
OPT_KEY( Real, temperature )
OPT_KEY( Integer, n_sample )
OPT_KEY( Integer, n_sample_beta )
OPT_KEY( Integer, num_torsion_list1 )
OPT_KEY( Real, RBangle_range )
OPT_KEY( Real, RBangle_increment )
OPT_KEY( Real, torsion_range )
OPT_KEY( Real, torsion_increment )
OPT_KEY( Real, score_cutoff )
OPT_KEY( Real, contact_cutoff )
OPT_KEY( Real, steric_dist_cutoff )
OPT_KEY( Real, score_diff_cut )
OPT_KEY( Integer, min_contacts )
OPT_KEY( Integer, fixed_pair_state_number )
OPT_KEY( Integer, new_pair_state_number )
OPT_KEY( Integer, num_new_pair_states )
OPT_KEY( Integer, min_hbonds )
OPT_KEY( Real, fa_rep_cutoff )
OPT_KEY( Real, chi1 )
OPT_KEY( Real, chi2 )
OPT_KEY( Boolean, virtualize_phosphate )
OPT_KEY( Boolean, superimpose_over_all_res )
OPT_KEY( Boolean, south1 )
OPT_KEY( Boolean, south2 )
OPT_KEY( Boolean, syn_chi1 )
OPT_KEY( Boolean, syn_chi2 )
OPT_KEY( Boolean, o2prime_trials )
OPT_KEY( Boolean, all_new_pair_states )
OPT_KEY( Boolean, only_positive_Z )
OPT_KEY( Boolean, cycle_axes )
OPT_KEY( Boolean, do_not_rotate_base2 )
OPT_KEY( Boolean, cluster_poses )
OPT_KEY( Boolean, filter_base_stack_direction )
OPT_KEY( Boolean, expand_chi )
OPT_KEY( Boolean, base_doublet_rmsd )
OPT_KEY( Boolean, cluster_rigid_body_settings )
OPT_KEY( Boolean, assign_to_clusters )
OPT_KEY( Boolean, two_base_pairs )
OPT_KEY( Boolean, two_base_pairs_via_loop_closure )
OPT_KEY( Boolean, test_ideal )
OPT_KEY( Boolean, fine_torsions )
OPT_KEY( Boolean, super_fine_torsions )
OPT_KEY( Boolean, gzip_out )
OPT_KEY( Boolean, close_loop_test )
OPT_KEY( Boolean, check_determinant )
OPT_KEY( Boolean, reverse_rbs )
OPT_KEY( Boolean, switch_chainbreak )
OPT_KEY( Boolean, dinucleotide )
OPT_KEY( Boolean, delta_chi_correction )
OPT_KEY( Boolean, finely_sample_base_pair )
OPT_KEY( Boolean, base_pair_to_base_pair )
OPT_KEY( Boolean, just_output_score )
OPT_KEY( Boolean, force_antiparallel_bases )
OPT_KEY( Boolean, force_parallel_bases )
OPT_KEY( Boolean, center_around_native )
OPT_KEY( Boolean, ignore_o2prime_hbonds_in_filter )
OPT_KEY( Boolean, assign_WC_edges )
OPT_KEY( Boolean, reverse_doublet )
OPT_KEY( Real, x_min )
OPT_KEY( Real, x_max )
OPT_KEY( Real, x_increment )
OPT_KEY( Real, y_min )
OPT_KEY( Real, y_max )
OPT_KEY( Real, y_increment )
OPT_KEY( Real, z_min )
OPT_KEY( Real, z_max )
OPT_KEY( Real, z_increment )
OPT_KEY( Real, alpha_min )
OPT_KEY( Real, alpha_max )
OPT_KEY( Real, alpha_increment )
OPT_KEY( Real, cosbeta_min )
OPT_KEY( Real, cosbeta_max )
OPT_KEY( Real, cosbeta_increment )
OPT_KEY( Real, gamma_min )
OPT_KEY( Real, gamma_max )
OPT_KEY( Real, gamma_increment )
OPT_KEY( Real, rmsd_cutoff )
OPT_KEY( Real, rep_cutoff )
OPT_KEY( String,  seq )
OPT_KEY( String,  rigid_body_samples )
OPT_KEY( String,  reference_rigid_body_samples )
OPT_KEY( String,  reference_rigid_body_samples_fixed_pair )
OPT_KEY( String,  reference_rigid_body_samples_new_pair )
OPT_KEY( String,  input_base1_torsion_M_v_lists )
OPT_KEY( String,  output_base1_torsion_M_v_lists )
OPT_KEY( String,  input_base2_torsion_M_v_lists )
OPT_KEY( String,  output_base2_torsion_M_v_lists )


///////////////////////////////////////////////////////////////
utility::vector1< Real >
get_suite_ideal_A_form_torsions(){

	using namespace chemical::rna;
	static RNA_FittedTorsionInfo rna_fitted_torsion_info;

	utility::vector1< Real >  ideal_A_form_torsions;

	ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_epsilon_north()[1].center );
	ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_minus()[1].center );
  ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_alpha()[1].center );
  ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_beta()[1].center );
  ideal_A_form_torsions.push_back( rna_fitted_torsion_info.gaussian_parameter_set_gamma()[1].center );

	return ideal_A_form_torsions;
}

///////////////////////////////////////////////////////////////
void
apply_ideal_A_form_torsions( pose::Pose & pose ){

	using namespace chemical::rna;
	using namespace id;
	using namespace protocols::swa;
	using namespace protocols::stepwise::sampling::rna;

	// For now assume delta, chi are at ideal values.
	static RNA_FittedTorsionInfo rna_fitted_torsion_info;
	StepWiseRNA_BaseSugarRotamerOP base_sugar_rotamer = new StepWiseRNA_BaseSugarRotamer( ANTI, NORTH, rna_fitted_torsion_info, 20.0, 3 );

	// This is kind of ugly. Need to know that the second BaseSugarRotamer  is the "ideal" one.
	// Anyway, stick with this for now, then consult with Parin on better design.
	base_sugar_rotamer->get_next_rotamer();
	base_sugar_rotamer->get_next_rotamer();

	for ( Size i = 1; i <= pose.total_residue(); i++ ){
		pose.set_torsion( TorsionID( i, id::BB, DELTA ), base_sugar_rotamer->delta() );
		pose.set_torsion( TorsionID( i, id::CHI, NU2 - NUM_RNA_MAINCHAIN_TORSIONS ), base_sugar_rotamer->nu2() );
		pose.set_torsion( TorsionID( i, id::CHI, NU1 - NUM_RNA_MAINCHAIN_TORSIONS ), base_sugar_rotamer->nu1() );
	}

	utility::vector1< Real > torsion_set = get_suite_ideal_A_form_torsions();
	pose.set_torsion( TorsionID( 1,   id::BB, EPSILON ),  torsion_set[1] );
	pose.set_torsion( TorsionID( 1,   id::BB, ZETA ),     torsion_set[2] );
	pose.set_torsion( TorsionID( 1,   id::BB, ALPHA ),    torsion_set[3] );
	pose.set_torsion( TorsionID( 1,   id::BB, BETA ),     torsion_set[4] );
	pose.set_torsion( TorsionID( 1,   id::BB, GAMMA ),    torsion_set[5] );
	//		pose.set_torsion( TorsionID( moving_suite+1, id::BB, ALPHA ),    torsion_set[3] );
	//		pose.set_torsion( TorsionID( moving_suite+1, id::BB, BETA ),     torsion_set[4] );
	//		pose.set_torsion( TorsionID( moving_suite+1, id::BB, GAMMA ),    torsion_set[5] );
}


///////////////////////////////////////////////////////////////////////////////
void
apply_south_syn_to_dinucleotide_pose( pose::Pose & pose ){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical::rna;
	using namespace id;

	RNA_FittedTorsionInfo rna_fitted_torsion_info;

	// fix the pucker. Probably should create a function "apply_ideal_c3endo_sugar_coords".
	if ( option[ south1 ]() ) {
		apply_ideal_c2endo_sugar_coords( pose, 1 );
	}
	if ( option[ south2 ]() ) {
		apply_ideal_c2endo_sugar_coords( pose, 2 );
	}


	// fix the chi angles
	if ( option[ syn_chi1 ]() ){
		if( option[ south1 ]() ){
			pose.set_torsion( TorsionID( 1,   id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ),
												rna_fitted_torsion_info.gaussian_parameter_set_chi_south()[2].center );
		} else {
			pose.set_torsion( TorsionID( 1,   id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ),
												rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[2].center );
		}
	} else {
		if( option[ south1 ]() ){
			pose.set_torsion( TorsionID( 1,   id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ),
												rna_fitted_torsion_info.gaussian_parameter_set_chi_south()[1].center );
		} else {
			pose.set_torsion( TorsionID( 1,   id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ),
												rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[1].center );
		}
	}

	if ( option[ syn_chi2 ]() ){
		if( option[ south2 ]() ){
			pose.set_torsion( TorsionID( 2,   id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ),
												rna_fitted_torsion_info.gaussian_parameter_set_chi_south()[2].center );
		} else {
			pose.set_torsion( TorsionID( 2,   id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ),
												rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[2].center );
		}
	} else {
		if( option[ south2 ]() ){
			pose.set_torsion( TorsionID( 2,   id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ),
												rna_fitted_torsion_info.gaussian_parameter_set_chi_south()[1].center );
		} else {
			pose.set_torsion( TorsionID( 2,   id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ),
												rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[1].center );
		}
	}

	if ( option[ chi1 ].user() ){
			pose.set_torsion( TorsionID( 1,   id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ),
												option[ chi1 ]() );
	}
	if ( option[ chi2 ].user() ){
			pose.set_torsion( TorsionID( 2,   id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ),
												option[ chi2 ]() );
	}

	std::cout << "### CHI ANGLES: " << pose.torsion( TorsionID( 1, id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ) ) <<
		" " << pose.torsion( TorsionID( 2, id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ) ) << std::endl;


}

///////////////////////////////////////////////////////////////
void
initialize_base_pair( pose::Pose & pose,
											utility::vector1< Size >  & moving_res1,
											utility::vector1< Size >  & moving_res2
											){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace pose;
	using namespace id;
	using namespace protocols::swa;
	using namespace protocols::stepwise::sampling::rna;
	using namespace chemical::rna;

	ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	// Create C-G base pair.
	make_pose_from_sequence( pose, option[ seq ](),	*rsd_set );

	apply_ideal_A_form_torsions( pose );

	apply_south_syn_to_dinucleotide_pose( pose );

	FoldTree f( 2 );
	f.new_jump( 1, 2, 1);
	f.set_jump_atoms( 1,
										chemical::rna::chi1_torsion_atom( pose.residue( 1) ),
										chemical::rna::chi1_torsion_atom( pose.residue( 2) )   );
	pose.fold_tree( f );
	//	pose.dump_pdb( "start.pdb" );

	// Virtualize backbone
	if ( option[ virtualize_phosphate ]() ) {
		add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );
		add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 2 );
	} else {
		add_variant_type_to_pose_residue( pose, "VIRTUAL_BACKBONE_EXCEPT_C1PRIME", 1 );
		add_variant_type_to_pose_residue( pose, "VIRTUAL_BACKBONE_EXCEPT_C1PRIME", 2 );
	}
	//	pose.dump_pdb( "virtualize.pdb" );

	moving_res1 = make_vector1( 1 );
	moving_res2 = make_vector1( 2 );

	translate_and_rotate_residue_to_origin( pose, 1, moving_res1, false /*do not rotate*/ );
	translate_and_rotate_residue_to_origin( pose, 2, moving_res2, option[ do_not_rotate_base2 ]() );

}



//////////////////////////////////////////////////////////////////////////////////////
Matrix
cycle( Matrix const & M0 ){
	Matrix M;
	for ( Size k = 1; k <= 3; k++ ){
		M(k,1) = M0(k,3);
		M(k,2) = M0(k,2);
		M(k,3) = M0(k,1);
	}
	return M;
}


//////////////////////////////////////////////////////////////////////////////////////
protocols::stepwise::RigidBodySamplerOP
initialize_rigid_body_sampler( utility::vector1< Size > const & moving_res1,
															 utility::vector1< Size > const & moving_res2 ){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace io::silent;
	using namespace protocols::swa;

	//	RNA_CentroidInfo rna_centroid_info;
	//	Vector centroid1 = rna_centroid_info.get_base_centroid( pose.residue(1) );
	//	Stub s1 = rna_centroid_info.get_base_coordinate_system( pose.residue(1), centroid1 );
	//	Matrix axes = s1.M;
	//	if ( option[ cycle_axes ]() ) axes = cycle( axes );

	// Hey is it necessary to input centroid?
	// Also, are axes necessary? The code basically assumes these are the origin and the identity matrix.
	rigid::RigidBodySamplerOP rigid_body_sampler = new rigid::RigidSampler( moving_res1, moving_res2 );

	//////////////////////////////////////////////////////////////
	// Initialize from command line.
	//////////////////////////////////////////////////////////////
	rigid_body_sampler->set_n_sample_alpha_full_range( option[ n_sample ]() );
	rigid_body_sampler->set_n_sample_cosbeta_full_range( option[ n_sample_beta ]() );
	rigid_body_sampler->set_n_sample_gamma_full_range( option[ n_sample ]() );

	if ( option[ alpha_min ].user() ){
		 rigid_body_sampler->set_alpha_min( option[ alpha_min ]() );
		 rigid_body_sampler->set_alpha_max( option[ alpha_max ]() );
		 rigid_body_sampler->set_alpha_increment( option[ alpha_increment ]() );

		 rigid_body_sampler->set_cosbeta_min( option[ cosbeta_min ]() );
		 rigid_body_sampler->set_cosbeta_max( option[ cosbeta_max ]() );
		 rigid_body_sampler->set_cosbeta_increment( option[ cosbeta_increment ]() );

		 rigid_body_sampler->set_gamma_min( option[ gamma_min ]() );
		 rigid_body_sampler->set_gamma_max( option[ gamma_max ]() );
		 rigid_body_sampler->set_gamma_increment( option[ gamma_increment ]() );
	}

	rigid_body_sampler->set_translation_sample( option[ box_radius ](), option[ xyz_sample ]() );
	if ( option[ x_min ].user() ){
		rigid_body_sampler->set_x_min( option[ x_min ]() );
		rigid_body_sampler->set_x_max( option[ x_max ]() );
		rigid_body_sampler->set_x_increment( option[ x_increment ]() );

		rigid_body_sampler->set_y_min( option[ y_min ]() );
		rigid_body_sampler->set_y_max( option[ y_max ]() );
		rigid_body_sampler->set_y_increment( option[ y_increment ]() );

		rigid_body_sampler->set_z_min( option[ z_min ]() );
		rigid_body_sampler->set_z_max( option[ z_max ]() );
		rigid_body_sampler->set_z_increment( option[ z_increment ]() );
	}


	rigid_body_sampler->set_score_cutoff( option[ score_cutoff ]() );

	rigid_body_sampler->set_min_num_contacts( option[ min_contacts ]() );
	rigid_body_sampler->set_contact_cutoff( option[ contact_cutoff ]() );
	rigid_body_sampler->set_steric_dist_cutoff( option[ steric_dist_cutoff ]() );

	if ( option[ force_antiparallel_bases ]() ) {
		rigid_body_sampler->force_antiparallel();
		rigid_body_sampler->force_coplanar();
	}
	if ( option[ force_parallel_bases ]() ) {
		rigid_body_sampler->force_parallel();
		rigid_body_sampler->force_coplanar();
	}

	rigid_body_sampler->set_min_hbonds( option[ min_hbonds ]() );
	rigid_body_sampler->set_ignore_o2prime_hbonds_in_filter( option[ ignore_o2prime_hbonds_in_filter ]() );
	rigid_body_sampler->set_assign_WC_edges( option[ assign_WC_edges ]() );
	rigid_body_sampler->set_fa_rep_cutoff( option[ fa_rep_cutoff ]() );

	ScoreFunctionOP scorefxn = get_score_function();
	rigid_body_sampler->set_score_function( scorefxn );
	rigid_body_sampler->set_o2prime_trials( option[ o2prime_trials ]() );

	SilentFileDataOP sfd;
	if ( option[ out::file::silent ].user() ){
		sfd = new SilentFileData;
		rigid_body_sampler->set_silent_file_data( sfd );
	}

	return rigid_body_sampler;

}

///////////////////////////////////////////////////////////////
void
output_stuff( protocols::stepwise::RigidBodySamplerOP rigid_body_sampler ){

	using namespace options;
	using namespace options::OptionKeys;

	std::string const outfile = option[ out::file::o ]();
	utility::io::ozstream out( outfile);
	std::cout << "outputting results to " << outfile << std::endl;
	rigid_body_sampler->output_results( out );
	out.close();

	std::string const histogram_outfile = outfile+".histogram";
	utility::io::ozstream hout( histogram_outfile);
	std::cout << "outputting histogram to " << histogram_outfile << std::endl;
	rigid_body_sampler->output_histogram( hout );
	hout.close();
}

////////////////////////////////////////////////////////////////
void
cluster_silent_file_data( io::silent::SilentFileDataOP & sfd ){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace io::silent;
	using namespace protocols::swa;

	std::string const silent_file = option[ out::file::silent ]();
	std::cout << "About to cluster..." << std::endl;
	// Cluster lowest energy states and output.
	protocols::stepwise::StepWiseLegacyClusterer stepwise_clusterer( sfd );
	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer.set_max_decoys( max_decoys );
	stepwise_clusterer.set_cluster_by_all_atom_rmsd( true );
	stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );

	if ( !option[ superimpose_over_all_res ]() ){
		utility::vector1< Size > moving_res2 = make_vector1( 2 );
		stepwise_clusterer.set_calc_rms_res( moving_res2 );
	}

	Real cluster_radius( 0.25 );
	if ( option[ OptionKeys::cluster::radius ].user() ) cluster_radius = option[ OptionKeys::cluster::radius ]();
	stepwise_clusterer.set_cluster_radius( cluster_radius	);
	stepwise_clusterer.set_score_diff_cut( option[ score_diff_cut ]() );

	std::cout << "Clustering..." << std::endl;
	stepwise_clusterer.cluster();

	sfd = stepwise_clusterer.silent_file_data();

}

///////////////////////////////////////////////////////////////
void
adjust_pose_chi( core::pose::Pose & pose,
								 Size const res,
								 Real const delstd ){

	using namespace core::chemical::rna;
	using namespace protocols::stepwise::sampling::rna;
	using namespace core::id;

	static RNA_FittedTorsionInfo rna_fitted_torsion_info;

	PuckerState pucker = Get_residue_pucker_state( pose, res );

	Real const current_chi = pose.chi( res );
	bool is_syn = ( numeric::principal_angle_degrees( current_chi ) < 0.0 );
	Size const which_chi = ( is_syn ?  2 : 1 );

	//	std::cout << "CHECK_THIS. RES: " << res << " PUCKER: " << pucker << " CHI: " << which_chi << std::endl;

	if ( pucker == NORTH ){
		pose.set_torsion( TorsionID( res,   id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ),
											rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[which_chi].center +
											delstd * rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[which_chi].width );
	} else {
		pose.set_torsion( TorsionID( res,   id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS ),
											rna_fitted_torsion_info.gaussian_parameter_set_chi_south()[which_chi].center +
											delstd * rna_fitted_torsion_info.gaussian_parameter_set_chi_south()[which_chi].width );
	}

}

///////////////////////////////////////////////////////////////
void
expand_chi_for_silent_structs( core::io::silent::SilentFileDataOP & sfd, core::scoring::ScoreFunctionOP scorefxn ){

	using namespace core::io::silent;

	SilentFileDataOP sfd_new = new SilentFileData;

	Size count( 0 );
	for ( core::io::silent::SilentFileData::iterator iter = sfd->begin(),
					end = sfd->end(); iter != end; ++iter ) {

		pose::Pose pose;
		iter->fill_pose( pose );

		for ( int i = -1; i <= 1; i++ ){
			adjust_pose_chi( pose, 1, static_cast<Real>(i) );

			for ( int j = -1; j <= 1; j++ ){
				adjust_pose_chi( pose, 2, static_cast<Real>(j) );

				(*scorefxn)( pose );
				count++;
				std::string tag = "S_" + ObjexxFCL::lead_zero_string_of( count, 5 );
				BinarySilentStruct s( pose, tag ); // will this copy in scores?
				sfd_new->add_structure( s );

			}
		}

	}

	sfd = sfd_new;

}


///////////////////////////////////////////////////////////////
void
split_silent_file_data( core::io::silent::SilentFileDataOP & sfd,
												utility::vector1< core::io::silent::SilentFileDataOP > & split_sfds )

{
	using namespace core::io::silent;
	split_sfds.clear();
	for ( core::io::silent::SilentFileData::iterator iter = sfd->begin(),
					end = sfd->end(); iter != end; ++iter ) {
		SilentFileDataOP sfd_new= new SilentFileData;
		sfd_new->add_structure(  *iter );
		split_sfds.push_back( sfd_new );
	}
}

///////////////////////////////////////////////////////////////
void
define_states_test(){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace io::silent;
	using namespace protocols::swa;


	///////////////////////////////////////////////
	// initialization
	///////////////////////////////////////////////
	Pose pose;
	utility::vector1< Size > moving_res1, moving_res2;
	initialize_base_pair( pose, moving_res1, moving_res2 );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );
	rigid::RigidBodySamplerOP rigid_body_sampler = initialize_rigid_body_sampler( moving_res1, moving_res2 );

	///////////////////////////////////////////////
	// Rigid body sample. Keep track of total number of states so that we can extract a Kd
	// Use input parameters to define fineness of sampling -- will look for convergence.
	// Save lowest energy states.
	///////////////////////////////////////////////
	if ( option[ rigid_body_samples ].user() ){
		rigid_body_sampler->apply_input_samples( pose, option[ rigid_body_samples ]() );
	} else {
		rigid_body_sampler->do_the_sampling( pose );
	}

	///////////////////////////////////////////////
	// output!
	///////////////////////////////////////////////
	if ( option[ out::file::o ].user() )  output_stuff( rigid_body_sampler );
	std::string silent_tag = option[ out::file::silent ]();

	///////////////////////////////////////////////
	// Clustering!
	///////////////////////////////////////////////
	SilentFileDataOP sfd = rigid_body_sampler->silent_file_data();
	if ( sfd ){

		if ( option[ cluster_poses ]() )	cluster_silent_file_data( sfd );

		if ( option[ split_silent_files ]() ) {
			utility::vector1< SilentFileDataOP > split_sfds;
			split_silent_file_data( sfd, split_sfds );

			for ( Size n = 1; n <= split_sfds.size(); n++ ) {
				sfd = split_sfds[ n ];
				if ( option[ expand_chi ]() ) expand_chi_for_silent_structs( sfd, rigid_body_sampler->score_function() );

				Size pos( silent_tag.find( ".out" ) );
				std::string filename = silent_tag.substr( 0, pos ) + "_" + ObjexxFCL::lead_zero_string_of(n,5)+".out";

				sfd->write_all( filename );
			}
		} else {
			if ( option[ expand_chi ]() ) expand_chi_for_silent_structs( sfd, rigid_body_sampler->score_function() );
			sfd->write_all( silent_tag );
		}


	}

	// Print out Kd.
	//Real const xyz_increment = option[ xyz_sample ]();
	//Real const concentration = (1.0 / 6.022e-4);
	//std::cout << "sampling concentration " << concentration << std::endl;

	//Real const additional_probability = Z * xyz_increment * xyz_increment * xyz_increment * alpha_increment / ( N_SAMPLE * N_SAMPLE_COSBETA * N_SAMPLE );
	//	std::cout << "additional_probability: " <<  additional_probability << std::endl;

	//Real const Kd = concentration * (1.0 / additional_probability);
	//	std::cout << "Kd " << Kd << std::endl;

}

///////////////////////////////////////////////////////////////
Real
msd_base_doublet(	Matrix const & MA,
									Matrix const & MB,
									Vector const & vA,
									Vector const & vB,
									Matrix const & moments ){

	Real dev2( 0.0 );

	dev2 += ( vA - vB ).length_squared(); // center of mass.
	dev2 += 2 * ( moments.trace() -  ( moments * MA.transposed() * MB ).trace() );


	return dev2;

}

///////////////////////////////////////////////////////////////
Real
rmsd_base_doublet( Matrix const & MA,
									 Matrix const & MB,
									 Vector const & vA,
									 Vector const & vB,
									 Matrix const & moments ){
	return std::sqrt( msd_base_doublet( MA, MB, vA, vB, moments ) );
}


///////////////////////////////////////////////////////////////
Real
msd_base_doublet_symmetric( Matrix const & MA,
														Matrix const & MB,
														Vector const & vA,
														Vector const & vB,
														Matrix const & moments1,
														Matrix const & moments2 ){
	Real dev2( 0.0 );

	// moments2 because we assume M1, M2, v1, v2 define rigid body transformation
	// of second base -- first base is assumed to be at origin of coordinate system.
	dev2 += msd_base_doublet( MA, MB, vA, vB, moments2 );

	Matrix MAT = MA.transposed();
	Matrix MBT = MB.transposed();
	Vector vAT = -1.0 * MAT * vA;
	Vector vBT = -1.0 * MBT * vB;
	dev2 += msd_base_doublet( MAT, MBT, vAT, vBT, moments1 );

	return dev2;
}

///////////////////////////////////////////////////////////////
Real
rmsd_base_doublet_symmetric( Matrix const & MA,
														 Matrix const & MB,
														 Vector const & vA,
														 Vector const & vB,
														 Matrix const & moments1,
														 Matrix const & moments2 ){

	Real const dev2 = msd_base_doublet_symmetric( MA, MB, vA, vB, moments1, moments2 );

	if ( dev2 <= 0.0 ) return 0.0;
	return std::sqrt( dev2 );

}

///////////////////////////////////////////////////////////////
Real
rmsd_over_base( conformation::Residue const & rsd1, conformation::Residue const & rsd2 ){

	Real dev2( 0.0 );
	Size natoms( 0 );

	for ( Size i = rsd1.first_sidechain_atom()+1; i <= rsd1.nheavyatoms(); i++ ) {
		if ( rsd1.is_virtual( i ) ) continue;
		dev2 += ( rsd1.xyz( i ) - rsd2.xyz( i ) ).length_squared();
		natoms++;
	}

	return std::sqrt( dev2 / natoms );

}

///////////////////////////////////////////////////////////////
Matrix
calculate_moments( conformation::Residue const & rsd ){

	using namespace conformation;
	using namespace protocols::swa;

	// This assumes that base is centered at 0.0,
	// and that it is aligned to desired x, y, z axes.

	// So, let's do the alignment.
	pose::Pose pose_scratch;
	pose_scratch.append_residue_by_bond( rsd );
	translate_and_rotate_residue_to_origin( pose_scratch, 1 );
	Residue const & rsd_scratch = pose_scratch.residue( 1 );

	Matrix M;
	Size natoms( 0 );
	for ( Size n = rsd_scratch.first_sidechain_atom()+1; n <= rsd_scratch.nheavyatoms(); n++ ) {
		if ( rsd_scratch.is_virtual( n ) ) continue;
		for ( Size i = 1; i <= 3; i++ ){
			for ( Size j = 1; j <= 3; j++ ){
				M( i, j ) += rsd_scratch.xyz( n )( i ) * rsd_scratch.xyz( n )( j );
			}
		}
		//		std::cout << "IS THIS A BASE ATOM? " << rsd_scratch.atom_name( n ) << std::endl;
		natoms++;
	}

	//	std::cout << "HEY! NATOMS: " << natoms << std::endl;

	return M/ natoms;

}

///////////////////////////////////////////////////////////////
void
base_doublet_rmsd_test(){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace io::silent;
	using namespace protocols::swa;


	///////////////////////////////////////////////
	// initialization
	///////////////////////////////////////////////
	Pose pose;
	utility::vector1< Size > moving_res1, moving_res2;
	initialize_base_pair( pose, moving_res1, moving_res2 );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );
	rigid::RigidBodySamplerOP rigid_body_sampler = initialize_rigid_body_sampler( moving_res1, moving_res2 );

	Pose pose1 = pose;
	Pose pose2 = pose;

	Vector axis1( 1.0, 0.0, 0.0 );
	Vector axis2( 0.0, 1.0, 0.0 );
	Vector axis3( 0.0, 0.0, 1.0 );
	Matrix M1, M2;

	// Define alpha, beta, gamma, x, y, z of pose 1
	Real alpha1( 50.0 ), beta1( 50.0 ), gamma1( 190.0 ), x1( -0.5 ), y1( 3.0 ), z1( 0.8 );
	rigid_body_sampler->apply_rigid_body_settings( pose1, pose, alpha1, beta1, gamma1, x1, y1, z1 );
  create_euler_rotation( M1, alpha1, beta1, gamma1, axis1, axis2, axis3 );

	// Define alpha, beta, gamma, x, y, z of pose 2
	Real alpha2( -5.0 ), beta2( 95.0 ), gamma2( 10.0 ), x2( -2.0 ), y2( -2.0 ), z2( -2.0 );
	rigid_body_sampler->apply_rigid_body_settings( pose2, pose, alpha2, beta2, gamma2, x2, y2, z2 );
	create_euler_rotation( M2, alpha2, beta2, gamma2, axis1, axis2, axis3 );

	// Calculate standard rmsd.
	std::cout << "RMSD over base atoms on residue 2: " << rmsd_over_base( pose1.residue(2), pose2.residue(2) ) << std::endl;

	// Calculate moments.
	Matrix moments = calculate_moments( pose.residue( 2 ) );

	// Calculate rmsd based on moments.
	std::cout << "RMSD based on Euler angles and c.o.m.: " << rmsd_base_doublet( M1, M2, Vector( x1, y1, z1), Vector( x2, y2, z2 ), moments ) << std::endl;
}


////////////////////////////////////////////////////////////////
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

////////////////////////////////////////////////////////////////
void
cluster_rigid_body_settings_test(){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace io::silent;
	using namespace protocols::swa;
	using namespace utility::io;

	std::string const outfile = option[ out::file::o ]();
	std::string const silent_file = option[ out::file::silent ]();

	// set up starting pose.
	Pose pose;
	utility::vector1< Size > moving_res1, moving_res2;
	initialize_base_pair( pose, moving_res1, moving_res2 );

	// determine moments.
	Matrix moments1 = calculate_moments( pose.residue( 1 ) );
	Matrix moments2 = calculate_moments( pose.residue( 2 ) );

	// read in starting list of rigid body settings.
	utility::vector1< utility::vector1< Real > > input_rigid_body_settings;
	std::string const infile = option[ rigid_body_samples ]();
	std::cout << "Reading rigid body settings from " << infile << std::endl;
	read_rigid_body_settings( infile, input_rigid_body_settings );

	// sort by energy.
	std::cout << "Setting up vectors and matrices for " <<  input_rigid_body_settings.size() << " rigid body settings. " << std::endl;
	std::list< std::pair< Real, Size > > all_energies;
	for ( Size n = 1; n <= input_rigid_body_settings.size(); n++ ){
		all_energies.push_back( std::make_pair( input_rigid_body_settings[n][7], n ) );
	}
	all_energies.sort();

	// save center-of-mass vectors, and rotation matrices, in order of energy.
	utility::vector1< Size > sort_index;
	utility::vector1< Vector > sort_v;
	utility::vector1< Matrix > sort_M;
	Vector axis1( 1.0, 0.0, 0.0 );
	Vector axis2( 0.0, 1.0, 0.0 );
	Vector axis3( 0.0, 0.0, 1.0 );
	Matrix M;
	for ( 	std::list< std::pair< Real, Size > >::iterator iter = all_energies.begin();
					iter != all_energies.end(); ++iter ) {

		Size const & i = iter->second;
		sort_index.push_back( i );

		utility::vector1< Real > const & rbs = input_rigid_body_settings[ i ];
		sort_v.push_back( Vector( rbs[4], rbs[5], rbs[6] ) );

		create_euler_rotation( M, rbs[1], rbs[2], rbs[3], axis1, axis2, axis3 );
		sort_M.push_back( M );
	}

	Real const rmsd_cutoff_ = option[ rmsd_cutoff ]();

	// new list of CLUSTERED rigid body settings.
	utility::vector1< Size > cluster_index;
	std::cout << "Clustering " << sort_index.size() << " models " << std::endl;
	Real rmsd( 0.0 );
	for ( Size n = 1; n <= sort_index.size(); n++ ){

		bool found_neighbor( false );
		Vector const & v1( sort_v[ n ] );
		Matrix const & M1( sort_M[ n ] );

		for ( Size j = 1; j <= cluster_index.size(); j++ ){
			Vector const & v2( sort_v[ cluster_index[j] ] );
			Matrix const & M2( sort_M[ cluster_index[j] ] );

			Real const rmsd = rmsd_base_doublet_symmetric( M1, M2, v1, v2, moments1, moments2 );
			//rmsd = rmsd_base_doublet( M1, M2, v1, v2, moments2 );

			if ( rmsd < rmsd_cutoff_ ) {
				//std::cout << "[ " << rmsd << " inside cutoff " << rmsd_cutoff << " ] ";
				found_neighbor = true;
				break;
			}
		}

		Real const & energy = input_rigid_body_settings[ sort_index[ n ] ][ 7 ];

		//		std::cout << "CHECKING " << sort_index[ n ] << " with score " << energy << " against list of size " << cluster_index.size();

		if ( !found_neighbor ){
			//std::cout << "... added " << std::endl;
			cluster_index.push_back( n );
			std::cout << "Adding cluster " << cluster_index.size() << " with energy " <<
				energy  << ". At model " << n << " out of " << sort_index.size() << std::endl;
		} else {
			//std::cout << "... not added. " << std::endl;
		}

	}

	// output CLUSTERED list, and silent files.
	// output silent file that has representative poses.
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );
	rigid::RigidBodySamplerOP rigid_body_sampler = initialize_rigid_body_sampler( moving_res1, moving_res2 );
	SilentFileData sfd;
	ozstream out( outfile );
	Pose pose_start = pose;
	ScoreFunctionOP scorefxn = get_score_function();

	std::cout << "Outputting  " << cluster_index.size() << " clusters to " << outfile << std::endl;
	for ( Size i = 1; i <= cluster_index.size(); i++ ){

		utility::vector1< Real > const & rbs = input_rigid_body_settings[ sort_index[ cluster_index[ i ]  ] ];
		for ( Size n = 1; n <= 8; n++ )	out << ' ' << rbs[ n ];
		out << std::endl;

		rigid_body_sampler->apply_rigid_body_settings( pose, pose_start, rbs[1],rbs[2],rbs[3],rbs[4],rbs[5],rbs[6] );
		(*scorefxn)(pose);
		std::string const tag = "S_" + ObjexxFCL::lead_zero_string_of( sort_index[ cluster_index[ i ] ], 6 );
		BinarySilentStruct s( pose, tag ); // this is RNA-centric -- could make it OK for proteins.
		sfd.write_silent_struct( s, silent_file, false );

	}

	out.close();

}


////////////////////////////////////////////////////////////////////////////
// This permits calculation of Kd of forming
// particular base pair arrangement, and is
// important in defining deltaG_fold for base pairs or base pair stacks.
////////////////////////////////////////////////////////////////////////////
void
finely_sample_base_pair_test(){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace io::silent;
	using namespace protocols::swa;
	using namespace utility::io;

	// set up pose
	Pose pose;
	utility::vector1< Size > moving_res1, moving_res2;
	initialize_base_pair( pose, moving_res1, moving_res2 );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );


	// read in input rigid body setting
	utility::vector1< utility::vector1< Real > > input_rigid_body_settings;
	std::string const infile = option[ rigid_body_samples ]();
	std::cout << "Reading rigid body settings from " << infile << std::endl;
	read_rigid_body_settings( infile, input_rigid_body_settings );
	utility::vector1< Real > const & rbs = input_rigid_body_settings[ option[ fixed_pair_state_number ]() ];
	Real const alpha_center =  rbs[1];
	Real const beta_center  =  rbs[2];
	Real const gamma_center =  rbs[3];
	Real const x_center     =  rbs[4];
	Real const y_center     =  rbs[5];
	Real const z_center     =  rbs[6];

	Pose pose_start = pose;

	// set up rigid body sampler
	rigid::RigidBodySamplerOP rigid_body_sampler = initialize_rigid_body_sampler( moving_res1, moving_res2 );

	Real const box_size      = option[ box_radius ]();
	Real const xyz_increment = option[ xyz_sample ]();
	Real const degree_range     = option[ RBangle_range ]();
	Real const degree_increment = option[ RBangle_increment ]();

	rigid_body_sampler->set_x_min( x_center - box_size );
	rigid_body_sampler->set_x_max( x_center + box_size );
	rigid_body_sampler->set_x_increment( xyz_increment );

	rigid_body_sampler->set_y_min( y_center - box_size );
	rigid_body_sampler->set_y_max( y_center + box_size );
	rigid_body_sampler->set_y_increment( xyz_increment );

	rigid_body_sampler->set_z_min( z_center - box_size );
	rigid_body_sampler->set_z_max( z_center + box_size );
	rigid_body_sampler->set_z_increment( xyz_increment );

	rigid_body_sampler->set_alpha_min( alpha_center - degree_range );
	rigid_body_sampler->set_alpha_max( alpha_center + degree_range );
	rigid_body_sampler->set_alpha_increment( degree_increment );

	rigid_body_sampler->set_gamma_min( gamma_center - degree_range );
	rigid_body_sampler->set_gamma_max( gamma_center + degree_range );
	rigid_body_sampler->set_gamma_increment( degree_increment );

	rigid_body_sampler->set_alpha_min( alpha_center - degree_range );
	rigid_body_sampler->set_alpha_max( alpha_center + degree_range );
	rigid_body_sampler->set_alpha_increment( degree_increment );

	rigid_body_sampler->set_rmsd_cutoff( option[ rmsd_cutoff ]() );

	Real const cosbeta_range = -cos( radians( beta_center ) ) + cos( radians( beta_center - degree_range ) );
	Real cosbeta_min = cos( radians( beta_center ) ) - cosbeta_range;
	//	if ( cosbeta_min < -1.0 ) cosbeta_min = -1.0;
	Real cosbeta_max = cos( radians( beta_center ) ) + cosbeta_range;
	//	if ( cosbeta_max >  1.0 ) cosbeta_max = 1.0;
	Real const cosbeta_increment = cosbeta_range * (degree_increment/degree_range);

	rigid_body_sampler->set_cosbeta_min( cosbeta_min );
	rigid_body_sampler->set_cosbeta_max( cosbeta_max );
	rigid_body_sampler->set_cosbeta_increment( cosbeta_increment );

	std::cout << "COSBETA sampled: " << cosbeta_min << " to " << cosbeta_max << "  in increments of: " << cosbeta_increment;

	// save the "starting" (central) configuration.
	rigid_body_sampler->apply_rigid_body_settings( pose, pose_start, rbs[1],rbs[2],rbs[3],rbs[4],rbs[5],rbs[6] );
	PoseOP native_pose = new Pose;
	*native_pose = pose;
	rigid_body_sampler->set_native_pose( native_pose );
	rigid_body_sampler->save_silent_struct( pose, "START" );

	// do the sampling
	pose = pose_start;
	rigid_body_sampler->do_the_sampling( pose );

	// output energies & phase space volumes of discovered poses.
	std::string const silent_file = option[ out::file::silent ]();
	rigid_body_sampler->output_silent_file( silent_file, option[ just_output_score ]() );


}



////////////////////////////////////////////////////////////////
void
fill_v_and_M( utility::vector1< utility::vector1< Real > > & rigid_body_settings,
							utility::vector1< Vector > & v_list,
							utility::vector1< Matrix > & M_list
							){

	using namespace protocols::swa;

	Vector axis1( 1.0, 0.0, 0.0 );
	Vector axis2( 0.0, 1.0, 0.0 );
	Vector axis3( 0.0, 0.0, 1.0 );
	Matrix M;
	for ( Size i = 1; i <= rigid_body_settings.size(); i++ ){

		utility::vector1< Real > const & rbs = rigid_body_settings[ i ];
		v_list.push_back( Vector( rbs[4], rbs[5], rbs[6] ) );

		create_euler_rotation( M, rbs[1], rbs[2], rbs[3], axis1, axis2, axis3 );
		M_list.push_back( M );
	}

}

////////////////////////////////////////////////////////////////
void
assign_rigid_body_settings_to_clusters_test(){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace io::silent;
	using namespace protocols::swa;
	using namespace utility::io;

	std::string const outfile = option[ out::file::o ]();

	// set up starting pose.
	Pose pose;
	utility::vector1< Size > moving_res1, moving_res2;
	initialize_base_pair( pose, moving_res1, moving_res2 );

	// determine moments.
	Matrix moments1 = calculate_moments( pose.residue( 1 ) );
	Matrix moments2 = calculate_moments( pose.residue( 2 ) );

	// read in starting list of rigid body settings.
	utility::vector1< utility::vector1< Real > > input_rigid_body_settings, reference_rigid_body_settings;

	std::string const infile = option[ rigid_body_samples ]();
	std::cout << "Reading rigid body settings from " << infile << std::endl;
	read_rigid_body_settings( infile, input_rigid_body_settings );

	std::string const infile_reference = option[ reference_rigid_body_samples ]();
	std::cout << "Reading rigid body settings from " << infile_reference << std::endl;
	read_rigid_body_settings( infile_reference, reference_rigid_body_settings );

	// save center-of-mass vectors, and rotation matrices, in order of energy.
	utility::vector1< Vector > v, v_reference;
	utility::vector1< Matrix > M, M_reference;
	fill_v_and_M( input_rigid_body_settings, v, M );
	fill_v_and_M( reference_rigid_body_settings, v_reference, M_reference );

	Real const rmsd_cutoff_ = option[ rmsd_cutoff ]();

	////////////////////////////////////////////////////////////////////////
	// Assign clusters.
	////////////////////////////////////////////////////////////////////////
	utility::vector1< utility::vector1< Size > > cluster_indices;
	utility::vector1< Size > blank_vector;
	for ( Size n = 1; n <= reference_rigid_body_settings.size(); n++ ) cluster_indices.push_back( blank_vector );

	std::cout << "Clustering " << input_rigid_body_settings.size() << " models into " << cluster_indices.size() << " clusters." << std::endl;

	for ( Size n = 1; n <= input_rigid_body_settings.size(); n++ ){

		Size closest_neighbor( 0 );
		Real closest_rmsd( 99999.9 );

		Vector const & v1( v[ n ] );
		Matrix const & M1( M[ n ] );

		for ( Size j = 1; j <= reference_rigid_body_settings.size(); j++ ){
			Vector const & v2( v_reference[ j ] );
			Matrix const & M2( M_reference[ j ] );

			Real const rmsd = rmsd_base_doublet_symmetric( M1, M2, v1, v2, moments1, moments2 );
			//rmsd = rmsd_base_doublet( M1, M2, v1, v2, moments2 );

			//if ( n == j ) std::cout << n << ' ' << j << "  rmsd: " << rmsd << std::endl;

			if ( rmsd < closest_rmsd || j == 1){
				closest_rmsd = rmsd;
				closest_neighbor = j;
			}

		}

		if ( closest_rmsd < rmsd_cutoff_ ) {
			cluster_indices[ closest_neighbor ].push_back( n );
		}

	}


	////////////////////////////////////////////////////////////////////////
	// Output
	////////////////////////////////////////////////////////////////////////

	std::cout << "Outputting  " << cluster_indices.size() << " clusters to " << outfile+".cluster*" << std::endl;

	for ( Size i = 1; i <= cluster_indices.size(); i++ ){

		utility::vector1< Size > const & cluster_index = cluster_indices[ i ];

		if ( cluster_index.size() == 0 ) continue;

		std::string const outfile_cluster = outfile + ".cluster" + lead_zero_string_of( i, 5 );

		ozstream out( outfile_cluster );

		for ( Size j = 1; j <= cluster_index.size(); j++ ){
			utility::vector1< Real > const & rbs = input_rigid_body_settings[ cluster_index[j] ];
			for ( Size n = 1; n <= 8; n++ )	out << ' ' << rbs[ n ];
			out << std::endl;
		}

		out.close();
	}

}

///////////////////////////////////////////////////////////////
void
setup_two_base_pair_pose( pose::Pose & pose ){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace kinematics;
	using namespace pose;
	using namespace io::silent;
	using namespace id;
	using namespace protocols::swa;
	using namespace protocols::stepwise::sampling::rna;
	using namespace chemical::rna;

	ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	make_pose_from_sequence( pose, "ccgg", *rsd_set );

	FoldTree f( 4 );
	f.new_jump( 1, 4, 2 );
	f.set_jump_atoms( 1,
										chemical::rna::chi1_torsion_atom( pose.residue( 1 ) ),
										chemical::rna::chi1_torsion_atom( pose.residue( 4 ) )   );
	pose.fold_tree( f );

	add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );
	add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 3 );

	utility::vector1< Size > strand1_res = make_vector1( 1, 2 );
	utility::vector1< Size > strand2_res = make_vector1( 3, 4 );
	translate_and_rotate_residue_to_origin( pose, 1, strand1_res );
	translate_and_rotate_residue_to_origin( pose, 4, strand2_res );


	//Need to setup starting base pair with user-input rigid body setting.
	utility::vector1< utility::vector1< Real > > reference_rigid_body_settings;
	std::string const infile_reference = option[ reference_rigid_body_samples_fixed_pair ]();
	std::cout << "Reading rigid body settings from " << infile_reference << std::endl;
	read_rigid_body_settings( infile_reference, reference_rigid_body_settings );

	rigid::RigidBodySamplerOP rigid_body_sampler = initialize_rigid_body_sampler( strand1_res, strand2_res );
	utility::vector1< Real > const & rbs = reference_rigid_body_settings[ option[ fixed_pair_state_number ]() ];
	rigid_body_sampler->apply_rigid_body_settings( pose, pose, rbs[1],rbs[2],rbs[3],rbs[4],rbs[5],rbs[6] );

	// For now assume delta, chi are at ideal values.
	RNA_FittedTorsionInfo rna_fitted_torsion_info;
	StepWiseRNA_BaseSugarRotamerOP base_sugar_rotamer = new StepWiseRNA_BaseSugarRotamer( ANTI, NORTH, rna_fitted_torsion_info, 20.0, 3 );

	// This is kind of ugly. Need to know that the second BaseSugarRotamer  is the "ideal" one.
	// Anyway, stick with this for now, then consult with Parin on better design.
	std::cout << "NEXT ROTAMER: " << base_sugar_rotamer->get_next_rotamer() << std::endl;
	std::cout << "NEXT ROTAMER: " << base_sugar_rotamer->get_next_rotamer() << std::endl;

	for ( Size i = 1; i <= pose.total_residue(); i++ ){
		pose.set_torsion( TorsionID( i, id::BB, DELTA ), base_sugar_rotamer->delta() );
		pose.set_torsion( TorsionID( i, id::CHI, NU2 - NUM_RNA_MAINCHAIN_TORSIONS ), base_sugar_rotamer->nu2() );
		pose.set_torsion( TorsionID( i, id::CHI, NU1 - NUM_RNA_MAINCHAIN_TORSIONS ), base_sugar_rotamer->nu1() );
	}

}

///////////////////////////////////////////////////////////////
Real
initialize_fa_rep( pose::Pose const & pose,
									 utility::vector1< Size > const & moving_suites,
									 scoring::ScoreFunctionOP rep_scorefxn ) {

	using namespace pose;
	using namespace kinematics;
	using namespace scoring;
	using namespace protocols::swa;

	Pose pose_expand = pose;

	for ( Size n = 1; n <= moving_suites.size(); n++ ){
		Size const jump_at_moving_suite = make_cut_at_moving_suite( pose_expand, moving_suites[n] );
		Jump j = pose_expand.jump( jump_at_moving_suite );
		j.set_translation( Vector( 1.0e4 * n, 0.0, 0.0 ) );
		pose_expand.set_jump( jump_at_moving_suite, j );
	}

	(*rep_scorefxn)( pose_expand );
	EnergyMap const & energy_map=pose_expand.energies().total_energies();
	return energy_map[ fa_rep ] * rep_scorefxn->get_weight( fa_rep );

}


///////////////////////////////////////////////////////////////
bool
check_clash( pose::Pose & pose,
						 Real const & fa_rep_score_baseline,
						 Real const & rep_cutoff_,
						 scoring::ScoreFunctionOP rep_scorefxn ){

	using namespace scoring;
	using namespace options;
	using namespace options::OptionKeys;

	(*rep_scorefxn)( pose );
	EnergyMap const & energy_map=pose.energies().total_energies();
	Real const fa_rep_score = energy_map[ fa_rep ] * rep_scorefxn->get_weight( fa_rep );

	//	std::cout << fa_rep_score << " " << fa_rep_score_baseline << std::endl;

	if ( (fa_rep_score - fa_rep_score_baseline) > rep_cutoff_ ) return false;

	static Real const tolerance( 1.0e-3 );
	if ( (fa_rep_score - fa_rep_score_baseline) < -1.0 * tolerance ) {
		std::cout << fa_rep_score << " " << fa_rep_score_baseline << std::endl;
		//		utility_exit_with_message( "Weird fa_rep?" );
	}

	return true;
}

///////////////////////////////////////////////////////////////
void
save_torsions( pose::Pose const & pose,
							 Size const & moving_suite,
							 utility::vector1< utility::vector1< Real > > & torsion_list ){

	using namespace id;
	using namespace chemical::rna;

	utility::vector1< Real > torsion_set;

	torsion_set.push_back( pose.torsion( TorsionID( moving_suite,   id::BB, EPSILON ) ) );
	torsion_set.push_back( pose.torsion( TorsionID( moving_suite,   id::BB, ZETA ) ) );
	torsion_set.push_back( pose.torsion( TorsionID( moving_suite+1, id::BB, ALPHA ) ) );
	torsion_set.push_back( pose.torsion( TorsionID( moving_suite+1, id::BB, BETA ) ) );
	torsion_set.push_back( pose.torsion( TorsionID( moving_suite+1, id::BB, GAMMA ) ) );
	torsion_list.push_back( torsion_set );
}



///////////////////////////////////////////////////////////////
void
save_torsion_M_v( pose::Pose const & pose,
									Size const & moving_suite,
									Size const & moving_base,
									utility::vector1< utility::vector1< Real > > & torsion_list,
									utility::vector1< Matrix > & M_list,
									utility::vector1< Vector > & v_list ) {

	using namespace kinematics;
	using namespace protocols::swa;

	save_torsions( pose, moving_suite, torsion_list );

	static Vector centroid;
	static Matrix M;
	get_base_centroid_and_rotation_matrix( pose, moving_base, centroid, M);
	M_list.push_back( M );
	v_list.push_back( centroid );

}

////////////////////////////////////////////////////////////////////////
Size
get_bin_size(){

	using namespace options;
	using namespace options::OptionKeys;

	Size bin_size = option[ fine_torsions ]() ?  10 : 20; /*instead of default 20.0 -- so 32x more sampling*/;
	if ( option[ super_fine_torsions ]() ) bin_size = 5;
	if ( option[ bin_sample ].user() ) bin_size = option[ bin_sample ]();

	return bin_size;

}

///////////////////////////////////////////////////////////////
void
sample_new_base_in_two_base_pair_pose( pose::Pose & pose,
																			 Size const which_strand,
																			 utility::vector1< utility::vector1< Real > > & torsion_list,
																			 utility::vector1< Matrix > & M_list,
																			 utility::vector1< Vector > & v_list,
																			 scoring::ScoreFunctionOP rep_scorefxn	) {

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace pose;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace protocols::stepwise::sampling::rna;

	///////////////////////////////////////////////////////////////////////////////////
	// Sample conformation [epsilon,zeta,alpha,beta,gamma] of suite in second strand
	///////////////////////////////////////////////////////////////////////////////////
	Pose pose_start = pose;

	Size moving_suite, moving_base;
	if ( which_strand == 2 ){
		moving_suite = 3;
		moving_base = 3;
		add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESIDUE", 2 ); // this residue is not in the game.
	} else {
		assert( which_strand == 1);
		moving_suite = 1;
		moving_base = 2;
		add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESIDUE", 3 ); // this residue is not in the game.
	}

	utility::vector1< Size > suite_res_list = make_vector1( moving_suite );
	//	suite_res_list.push_back( 1 );

	Real const bin_size = get_bin_size();
	StepWiseRNA_RotamerGeneratorWrapperOP rotamer_generator = new StepWiseRNA_RotamerGeneratorWrapper( pose,
																																																			 suite_res_list,
																																																			 false /*sample_sugar_and_base1*/,
																																																			 false /*sample_sugar_and_base2*/,
																																																			 bin_size );
	Size count( 0 );
	Size count_saved( 0 );

	Real const fa_rep_score_baseline = initialize_fa_rep( pose, make_vector1( moving_suite ), rep_scorefxn );
	Real const rep_cutoff_ = option[ rep_cutoff ]();

	while( rotamer_generator->has_another_rotamer() ){

		count++;

		utility::vector1< Torsion_Info > current_rotamer = rotamer_generator->get_next_rotamer();
		apply_rotamer( pose, current_rotamer);

		// Disallow steric clashes.
		if ( !check_clash( pose, fa_rep_score_baseline, rep_cutoff_, rep_scorefxn  ) ) continue;
		save_torsion_M_v( pose, moving_suite, moving_base, torsion_list, M_list, v_list );

		// Save torsion angles, centroid, and base coordinate system.
		//save_suite_and_rb_info( pose, moving_suite, 3 /* rigid_body_base */, torsion_list, M_list, v_list );
		count_saved++;

	}
	std::cout << "Total number of rotamers applied: " << count << std::endl;
	std::cout << "Total number that passed cuts:    " << count_saved << std::endl;

	pose = pose_start;
}

///////////////////////////////////////////////////////////////
void
input_torsion_M_v_lists( 	utility::vector1< utility::vector1< Real > > & torsion_list,
													utility::vector1< Matrix > & M_list,
													utility::vector1< Vector > & v_list,
													std::string const infile ){

	using namespace utility::io;

	std::cout << "Readin from file: " << infile << std::endl;

	torsion_list.clear();
	M_list.clear();
	v_list.clear();

	utility::vector1< Real > torsion_ids = make_vector1( 0.0, 0.0, 0.0, 0.0, 0.0 ); //epsilon,zeta,alpha,beta,gamma
	Vector v( 0.0, 0.0, 0.0 );
	Matrix M;

	utility::vector1< Real > torsion_set;
	for ( Size n = 1; n <= 5; n++ ) torsion_set.push_back( 0.0 );

	izstream input_stream( infile );
	while( input_stream >> torsion_set[1] ){

		input_stream >> torsion_set[2] >> torsion_set[3] >> torsion_set[4] >> torsion_set[5];

		for ( Size i = 1; i <= 3; i++ ){
			for ( Size k = 1; k <= 3; k++ ){
				input_stream >> M(i,k);
			}
		}

		input_stream >> v(1);
		input_stream >> v(2);
		input_stream >> v(3) >> skip;

		torsion_list.push_back( torsion_set );
		M_list.push_back( M );
		v_list.push_back( v );

	}

	std::cout << "Read in: " << torsion_list.size() << " samples from " << infile << std::endl;

}

///////////////////////////////////////////////////////////////
void
output_torsion_M_v_lists( 	utility::vector1< utility::vector1< Real > > const & torsion_list,
														utility::vector1< Matrix > const & M_list,
														utility::vector1< Vector > const & v_list,
														std::string const outfile ){

	using namespace utility::io;

	ozstream out( outfile );
	if ( not out ) utility_exit_with_message( "Could not make outfile" );

	for ( Size n = 1; n <= torsion_list.size(); n++ ){

		for ( Size i = 1; i <= torsion_list[n].size(); i++ ) out << ' ' << torsion_list[n][i];

		for ( Size i = 1; i <= 3; i++ ){
			for ( Size k = 1; k <= 3; k++ ){
				out << ' ' << M_list[n](i,k);
			}
		}

		for ( Size i = 1; i <= 3; i++ ) out << ' ' << v_list[n](i);
		out << std::endl;
	}

}


///////////////////////////////////////////////////////////////
void
apply_suite_torsions( utility::vector1< Real > const & torsion_set, pose::Pose & pose,  Size const moving_suite  ){

	using namespace id;
	using namespace chemical::rna;

	pose.set_torsion( TorsionID( moving_suite,   id::BB, EPSILON ),  torsion_set[1] );
	pose.set_torsion( TorsionID( moving_suite,   id::BB, ZETA ),     torsion_set[2] );
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, ALPHA ),    torsion_set[3] );
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, BETA ),     torsion_set[4] );
	pose.set_torsion( TorsionID( moving_suite+1, id::BB, GAMMA ),    torsion_set[5] );

}




///////////////////////////////////////////////////////////////
void
do_the_match(
						 utility::vector1< utility::vector1< utility::vector1< Real > > > & strand1_strand2_info_for_each_cluster,
						 utility::vector1< Size > const & strand2_index,
						 utility::vector1< Size > const & reference_index,
						 pose::Pose & pose,
						 Matrix const & moments1, Matrix const & moments2,
						 utility::vector1< Real > strand1_torsion_set,
						 Matrix const & M1, Vector const & v1,
						 utility::vector1< utility::vector1< Real > > const & torsion_list,
						 utility::vector1< Matrix > const & M_list,
						 utility::vector1< Vector > const & v_list,
						 utility::vector1< utility::vector1< Real > > const & reference_rigid_body_settings,
						 pose::Pose const & ideal_pose,
						 scoring::ScoreFunctionOP scorefxn,
						 std::string const & silent_file,
						 Real const rmsd_cutoff_,
						 Real const rep_cutoff_,
						 Real const fa_rep_score_baseline,
						 scoring::ScoreFunctionOP rep_scorefxn ){

	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace io::silent;
	using namespace id;
	using namespace protocols::swa;
	using namespace utility::io;
	using namespace chemical::rna;
	using namespace protocols::stepwise::sampling::rna;
	using namespace options;
	using namespace options::OptionKeys;


	//////////////////////////////////////////////////////////////////////////////////////////////////
	// Find torsion angle combinations that might fit into a cluster. Build those poses
	// and save energies.
	//////////////////////////////////////////////////////////////////////////////////////////////////
	for ( Size m = 1; m <= strand2_index.size(); m++ ){

		Size const i = strand2_index[ m ];

		if ( m % 1000 == 0 ) std::cout << "On torsion_set " << m << " out of " << strand2_index.size() << std::endl;

		Matrix const & M2 = M_list[ i ];
		Vector const & v2 = v_list[ i ];
		//		std::cout << "V2" << ' ' <<  v2(1) << ' ' <<  v2(2) << ' ' <<  v2(3) << std::endl;

		//Oooh, linear algebra.
		Matrix const M2_in_M1_frame = M1.transposed() * M2;
		Vector const v2_in_v1_frame = M1.transposed() * ( v2 - v1 );
		//		std::cout << "V2_IN_V1_frame" << ' ' <<  v2_in_v1_frame(1) << ' ' <<  v2_in_v1_frame(2) << ' ' <<  v2_in_v1_frame(3) << std::endl;

		// How close is this rigid body arrangement to any of the reference states?
		Real closest_rmsd( 0.0 );
		Size closest_neighbor( 0 );

		for ( Size n = 1; n <= reference_index.size(); n++ ){

			Size const j = reference_index[ n ];

			// probably should have these all precomputed.
			utility::vector1< Real > const & rbs = reference_rigid_body_settings[ j ];
			Matrix M_ref;
			create_euler_rotation( M_ref, rbs[1],rbs[2],rbs[3] );
			Vector const v_ref( rbs[4], rbs[5], rbs[6] );

			//			std::cout << "V_REF" << ' ' <<  v_ref(1) << ' ' <<  v_ref(2) << ' ' <<  v_ref(3) << std::endl;

			Real const rmsd = rmsd_base_doublet_symmetric( M2_in_M1_frame, M_ref, v2_in_v1_frame, v_ref,
																										 moments1, moments2 );

			if ( rmsd < closest_rmsd || n == 1){
				closest_rmsd = rmsd;
				closest_neighbor = j;
			}
		}

		if ( closest_rmsd < rmsd_cutoff_ ) {

			utility::vector1< Real > strand1_strand2_info =  strand1_torsion_set;
			for ( Size n = 1; n <= torsion_list[ i ].size(); n++ ) strand1_strand2_info.push_back( torsion_list[ i ][ n ] );
			apply_suite_torsions( torsion_list[ i ], pose, 3 );

			strand1_strand2_info.push_back( closest_rmsd );

			// Clash_check!
			if ( !check_clash( pose, fa_rep_score_baseline, rep_cutoff_, rep_scorefxn ) ) continue;

			Real const score = ( *scorefxn )( pose );
			strand1_strand2_info.push_back( score );

			//			if( score > 100.0 )	scorefxn->show( std::cout, pose );

			if ( silent_file.size() > 0  ){

				std::string const tag = "S_" + ObjexxFCL::lead_zero_string_of( i, 6 ) + "_" +
					ObjexxFCL::lead_zero_string_of( closest_neighbor, 6 );
				BinarySilentStruct s( pose, tag ); // this is RNA-centric -- could make it OK for proteins.
				Real const rmsd_to_ideal = all_atom_rmsd( pose, ideal_pose );
				s.add_energy( "all_rms", rmsd_to_ideal );

				s.add_energy( "log_vol", 10.0 * log( radians( get_bin_size() ) ) );

				static SilentFileData sfd;
				sfd.write_silent_struct( s, silent_file, option[ just_output_score ]()  );
			}

			// Save pose into silent file?
			strand1_strand2_info_for_each_cluster[ closest_neighbor ].push_back( strand1_strand2_info );
		}

	}
	}


///////////////////////////////////////////////////////////////
void
brute_force_matcher( 	utility::vector1< utility::vector1< utility::vector1< Real > > > & strand1_strand2_info_for_each_cluster,
											pose::Pose & pose,
											utility::vector1< Real > strand1_torsion_set,
											utility::vector1< utility::vector1< Real > > const & torsion_list,
											utility::vector1< Matrix > const & M_list,
											utility::vector1< Vector > const & v_list,
											utility::vector1< utility::vector1< Real > > const & reference_rigid_body_settings,
											pose::Pose const & ideal_pose,
											scoring::ScoreFunctionOP scorefxn,
											scoring::ScoreFunctionOP rep_scorefxn,
											std::string const & silent_file ){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace protocols::swa;

	Vector v1; 	Matrix M1;
	get_base_centroid_and_rotation_matrix( pose, 2, v1, M1);

	Matrix moments1 = calculate_moments( pose.residue( 2 ) );
	Matrix moments2 = calculate_moments( pose.residue( 3 ) );

	Real const rmsd_cutoff_ = option[ rmsd_cutoff ]();
	utility::vector1< Size > reference_index, strand2_index;
	for ( Size i = 1; i <= torsion_list.size(); i++ ) strand2_index.push_back( i );
	for ( Size i = 1; i <= reference_rigid_body_settings.size(); i++ ) reference_index.push_back( i );

	Real const fa_rep_score_baseline = initialize_fa_rep( pose, make_vector1( 1, 3 ), rep_scorefxn );
	Real const rep_cutoff_ = 3 * option[ rep_cutoff ](); //to allow weak base/pose and base/base clashes.

	do_the_match( strand1_strand2_info_for_each_cluster,
								strand2_index, reference_index,
								pose,
								moments1, moments2,
								strand1_torsion_set, M1, v1,
								torsion_list, M_list, v_list,
								reference_rigid_body_settings,
								ideal_pose,
								scorefxn,
								silent_file,
								rmsd_cutoff_,
								rep_cutoff_,
								fa_rep_score_baseline,
								rep_scorefxn );

}



///////////////////////////////////////////////////////////////
void
initialize_for_grid_matcher(	Real & box_size, Real & box_spacing,
															ObjexxFCL::FArray3D< utility::vector1< Size > > & grid_strand2_index,
															utility::vector1< Vector > const & v_list ){

	box_size = 20.0;
	box_spacing = 0.5;

	int const numbins = 2 * static_cast< int >( box_size / box_spacing );
	grid_strand2_index.dimension( numbins, numbins, numbins );

	//centroids of strand2 base, based on sampling torsions in strand2
	std::cout << "Filling grid_strand2_index... " << std::endl;
	for ( Size n = 1; n <= v_list.size(); n++ ){

		Real const & x = v_list[n](1);
		int const x_index = floor( ( x + box_size ) / box_spacing );
		if ( x_index < 1 || x_index > numbins ) continue;

		Real const & y = v_list[n](2);
		int const y_index = floor( ( y + box_size ) / box_spacing );
		if ( y_index < 1 || y_index > numbins ) continue;

		Real const & z = v_list[n](3);
		int const z_index = floor( ( z + box_size ) / box_spacing );
		if ( z_index < 1 || z_index > numbins ) continue;

		grid_strand2_index( x_index, y_index, z_index ).push_back( n );

	}


}



///////////////////////////////////////////////////////////////
void
grid_matcher( 	utility::vector1< utility::vector1< utility::vector1< Real > > > & strand1_strand2_info_for_each_cluster,
								pose::Pose & pose,
								utility::vector1< Real > strand1_torsion_set,
								utility::vector1< utility::vector1< Real > > const & torsion_list,
								utility::vector1< Matrix > const & M_list,
								utility::vector1< Vector > const & v_list,
								utility::vector1< utility::vector1< Real > > const & reference_rigid_body_settings,
								pose::Pose const & ideal_pose,
								scoring::ScoreFunctionOP scorefxn,
								scoring::ScoreFunctionOP rep_scorefxn,
								std::string const & silent_file,
								Real const box_size, Real const box_spacing,
								ObjexxFCL::FArray3D< utility::vector1< Size > > const &  grid_strand2_index
								){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace io::silent;
	using namespace id;
	using namespace protocols::swa;
	using namespace utility::io;
	using namespace chemical::rna;
	using namespace protocols::stepwise::sampling::rna;

	Vector v1; 	Matrix M1;
	get_base_centroid_and_rotation_matrix( pose, 2, v1, M1);

	Real const rmsd_cutoff_ = option[ rmsd_cutoff ]();
	Matrix moments1 = calculate_moments( pose.residue( 2 ) );
	Matrix moments2 = calculate_moments( pose.residue( 3 ) );

	Real const fa_rep_score_baseline = initialize_fa_rep( pose, make_vector1( 1, 3 ), rep_scorefxn );
	Real const rep_cutoff_ = 3 * option[ rep_cutoff ](); //to allow weak base/pose and base/base clashes.

	// Make a grid of potential base2 centroid locations.
	int const numbins = 2 * static_cast< int >( box_size / box_spacing );

	//////////////////////////////////////////////////////////////////////////
	// possible base2 centroid locations, based on location of base1 and
	// known reference rigid body arrangements of base2/base1.
	//////////////////////////////////////////////////////////////////////////
	FArray3D< utility::vector1< Size > > grid_reference_index( numbins, numbins, numbins );
	std::cout << "Filling reference_index... " << std::endl;
	Vector const & axis1 = M1.col_x();
	Vector const & axis2 = M1.col_y();
	Vector const & axis3 = M1.col_z();
	//
	// Can divide rmsd^2 by 2.0 because rmsd involves the sum of squares of
	//  two rmsds (superimpose on base1, deviation of base2;  plus superimpose on base2, deviation of base 1).
	// This should still be a rigorous upper bound! WAIT, NOT WORKING!!
	//
	Real const rmsd_cutoff2 = (rmsd_cutoff_ * rmsd_cutoff_); //  / 2.0;

	for ( Size n = 1; n <= reference_rigid_body_settings.size(); n++ ){

		utility::vector1< Real > const & rbs = reference_rigid_body_settings[n];

		Vector const v2_centroid = v1 + ( rbs[4] * axis1 + rbs[5]*axis2 + rbs[6] * axis3 );
		Real const & x = v2_centroid( 1 );
		Real const & y = v2_centroid( 2 );
		Real const & z = v2_centroid( 3 );

		if ( std::abs( x ) > box_size ) continue;
		if ( std::abs( y ) > box_size ) continue;
		if ( std::abs( z ) > box_size ) continue;

		// Need to fill in a *neighborhood* around this location.
		// Be conservative -- allow for worst case scenario.
		int const i_min = floor( ( x - rmsd_cutoff_ + box_size )/ box_spacing );
		int const i_max = floor( ( x + rmsd_cutoff_ + box_size )/ box_spacing ) + 1;

		int const j_min = floor( ( y - rmsd_cutoff_ + box_size )/ box_spacing );
		int const j_max = floor( ( y + rmsd_cutoff_ + box_size )/ box_spacing ) + 1;

		int const k_min = floor( ( z - rmsd_cutoff_ + box_size )/ box_spacing );
		int const k_max = floor( ( z + rmsd_cutoff_ + box_size )/ box_spacing ) + 1;

		FArray3D< bool > mask_box( (i_max-i_min+1), (j_max-j_min+1), (k_max-k_min+1), false );

		int i_center = floor( ( x + box_size ) / box_spacing ) - i_min + 1;
		int j_center = floor( ( y + box_size ) / box_spacing ) - j_min + 1;
		int k_center = floor( ( z + box_size ) / box_spacing ) - k_min + 1;
		mask_box( i_center, j_center, k_center ) = true;

		for ( int i = i_min; i <= i_max; i++ ){
			if ( i < 1 || i > numbins ) continue;
			Real const x_mask = -box_size + (i-1)*box_spacing;

			for ( int j = j_min; j <= j_max; j++ ){
				if ( j < 1 || j > numbins ) continue;
				Real const y_mask = -box_size + (j-1)*box_spacing;

				for ( int k = k_min; k <= k_max; k++ ){
					if ( k < 1 || k > numbins ) continue;
					Real const z_mask = -box_size + (k-1)*box_spacing;

					Real const dist2 = ( Vector( x_mask, y_mask, z_mask ) - v2_centroid ).length_squared();
					if ( dist2 < rmsd_cutoff2 ){
						for ( int o1 = 0; o1 <= 1; o1++ ){
							for ( int o2 = 0; o2 <= 1; o2++ ){
								for ( int o3 = 0; o3 <= 1; o3++ ){
									// this doesn't quite look right -- are the indices ever zero?
									mask_box( i - i_min + o1,
														j - j_min + o2,
														k - k_min + o3 ) = true;
								}
							}
						}
					}
				}
			}
		}

		for ( int i = i_min; i <= i_max; i++ ){
			if ( i < 1 || i > numbins ) continue;

			for ( int j = j_min; j <= j_max; j++ ){
				if ( j < 1 || j > numbins ) continue;

				for ( int k = k_min; k <= k_max; k++ ){
					if ( k < 1 || k > numbins ) continue;

					if ( !mask_box( i - i_min + 1,
													j - j_min + 1,
													k - k_min + 1 ) ) continue;

					grid_reference_index( i, j, k ).push_back( n );

				}
			}
		}

	}


	std::cout << "Checking matches... " << std::endl;
	for ( int i = 1; i <= numbins; i++ ){
		for ( int j = 1; j <= numbins; j++ ){
			for ( int k = 1; k <= numbins; k++ ){

					utility::vector1< Size > const & reference_index = grid_reference_index( i, j, k );
					if ( reference_index.size() == 0 ) continue;

					utility::vector1< Size > const & strand2_index = grid_strand2_index( i, j, k );
					if ( strand2_index.size() == 0 ) continue;

					do_the_match( strand1_strand2_info_for_each_cluster,
												strand2_index, reference_index,
												pose,
												moments1, moments2,
												strand1_torsion_set, M1, v1,
												torsion_list, M_list, v_list,
												reference_rigid_body_settings,
												ideal_pose,
												scorefxn,
												silent_file,
												rmsd_cutoff_,
												rep_cutoff_,
												fa_rep_score_baseline,
												rep_scorefxn );
			}
		}
	}

}



///////////////////////////////////////////////////////////////
void
output_strand1_strand2_info( std::string const & outfile_prefix,
														 utility::vector1< utility::vector1< utility::vector1< Real > > > const & strand1_strand2_info_for_each_cluster ){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace utility::io;

	//	std::string const outfile = outfile_prefix + ".cluster" + ObjexxFCL::lead_zero_string_of( n, 5 );
	std::string const outfile = outfile_prefix;
	ozstream out( outfile );

	for ( Size n = 1; n <= strand1_strand2_info_for_each_cluster.size(); n++ ){

		utility::vector1< utility::vector1< Real > > const & info = strand1_strand2_info_for_each_cluster[ n ];

		if (info.size() == 0 ) continue;

		std::cout << "Cluster " << n << ": Outputting " << info.size() << " torsion sets to " << outfile << std::endl;

		for ( Size i = 1; i <= info.size(); i++ ){
			out << n;
			for ( Size j = 1; j <= info[ i ].size(); j++ ) out << ' ' << info[i][j];
			out << std::endl;
		}

	}

	out.close();

	if ( option[ gzip_out ]() ) utility::file::gzip( outfile, true /*overwrite*/ );

}


///////////////////////////////////////////////////////////////
void
two_base_pairs_test(){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace io::silent;
	using namespace id;
	using namespace protocols::swa;
	using namespace utility::io;
	using namespace chemical::rna;
	using namespace protocols::stepwise::sampling::rna;

	/////////////////////////////////////////////
	/////////////////////////////////////////////
	// CONVENTION:
	//  1 == starting base pair, strand 1
	//  2 == new base pair, strand 1
	//   chainbreak
	//  3 == new base pair, strand 2
	//  4 == starting base pair, strand 2
	//   end of pose.
	/////////////////////////////////////////////
	/////////////////////////////////////////////

	std::string const outfile_prefix = option[ out::file::o ]();

	/////////////////////////////////////////////
	// set up starting pose.
	/////////////////////////////////////////////
	Pose pose;
	setup_two_base_pair_pose( pose );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	if ( option[ score::weights ].user() ) scorefxn = get_score_function();
	ScoreFunctionOP rep_scorefxn = new ScoreFunction;
	rep_scorefxn->set_weight( fa_rep, 0.12 );


	std::string silent_file = "";
	if ( option[ out::file::silent ].user() )	silent_file = option[ out::file::silent ]();

	// Need to iterate through suite conformations for residue 1.
	utility::vector1< utility::vector1< Real > > torsion_list1, torsion_list2;
	utility::vector1< Matrix > M_list1, M_list2;
	utility::vector1< Vector > v_list1, v_list2;

	//////////////////////////////////////////////////////////
	// either sample explicitly -- or compute once and read from disk?
	if ( option[ input_base1_torsion_M_v_lists ].user() ){
		input_torsion_M_v_lists( torsion_list1, M_list1, v_list1,  option[ input_base1_torsion_M_v_lists ]() );
	} else {
		sample_new_base_in_two_base_pair_pose( pose, 1 /*second base*/, torsion_list1, M_list1, v_list1, scorefxn );
	}
	if ( option[ output_base1_torsion_M_v_lists ].user() ) output_torsion_M_v_lists( torsion_list1, M_list1, v_list1, option[ output_base1_torsion_M_v_lists ] );

	//////////////////////////////////////////////////////////
	if ( option[ input_base2_torsion_M_v_lists ].user() ){
		input_torsion_M_v_lists( torsion_list2, M_list2, v_list2,  option[ input_base2_torsion_M_v_lists ]() );
	} else {
		sample_new_base_in_two_base_pair_pose( pose, 2 /*second base*/, torsion_list2, M_list2, v_list2, scorefxn );
	}
	if ( option[ output_base2_torsion_M_v_lists ].user() ) output_torsion_M_v_lists( torsion_list2, M_list2, v_list2, option[ output_base2_torsion_M_v_lists ] );

	///////////////////////////////////////////
	// Now sample base in first strand
	///////////////////////////////////////////
	// As test, put it into A-form conformation.
	utility::vector1< Real > strand1_torsion_set = get_suite_ideal_A_form_torsions();
	apply_suite_torsions( strand1_torsion_set, pose, 1 );
	apply_suite_torsions( strand1_torsion_set, pose, 3 );

	pose.dump_pdb( "ideal.pdb" );
	Pose ideal_pose = pose;

	//////////////////////////////////////////////
	utility::vector1< utility::vector1< Real > > reference_rigid_body_settings;
	std::string const infile_reference = option[ reference_rigid_body_samples_new_pair ]();
	std::cout << "Reading rigid body settings from " << infile_reference << std::endl;
	read_rigid_body_settings( infile_reference, reference_rigid_body_settings );

	utility::vector1< utility::vector1< utility::vector1< Real > > > strand1_strand2_info_for_each_cluster;
	utility::vector1<  utility::vector1< Real > > blank_vector;
	for ( Size n = 1; n <= reference_rigid_body_settings.size(); n++ ) strand1_strand2_info_for_each_cluster.push_back( blank_vector );

	////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Rigid body associated with strand 1's base involved in the second base pair (residue 2 in the pose).
	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	Real box_size, box_spacing;
	FArray3D< utility::vector1< Size > > grid_strand2_index;
	initialize_for_grid_matcher( box_size, box_spacing, grid_strand2_index, v_list2);

	if ( option[ test_ideal ]() ){
		// Assume that base 1 is in ideal A-form helix conformation -- don't sample it.
		// This was useful for testing.
		clock_t const time_start( clock() );

		//	brute_force_matcher( strand1_strand2_info_for_each_cluster,
		grid_matcher( strand1_strand2_info_for_each_cluster,
									pose,
									strand1_torsion_set,
									torsion_list2, M_list2, v_list2,
									reference_rigid_body_settings, ideal_pose, scorefxn, rep_scorefxn, silent_file,
									box_size, box_spacing, grid_strand2_index );
		std::cout << "Total time in matcher: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

	} else {

		Real total_time_in_matcher( 0.0 );

		Size num_torsion_list1_ = torsion_list1.size();
		if ( option[ num_torsion_list1 ].user() ) num_torsion_list1_ = option[ num_torsion_list1 ]();

		for( Size n = 1; n <= num_torsion_list1_; n++ ){
			std::cout << "strand torsion1: " << n << " of " <<  num_torsion_list1_ << std::endl;
			strand1_torsion_set = torsion_list1[n];
			apply_suite_torsions( strand1_torsion_set , pose, 1 );
			clock_t const time_start( clock() );
			grid_matcher( strand1_strand2_info_for_each_cluster,
										pose,
										strand1_torsion_set,
										torsion_list2, M_list2, v_list2,
										reference_rigid_body_settings, ideal_pose, scorefxn, rep_scorefxn, silent_file,
										box_size, box_spacing, grid_strand2_index );
			Real const time_in_matcher = static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC;
			std::cout << "Time in matcher: " <<  time_in_matcher << std::endl;
			total_time_in_matcher += time_in_matcher;
		}
		std::cout << "TOTAL time in matcher: " <<  total_time_in_matcher << std::endl;

	}


	// output: strand1_strand2_info_for_each_cluster
	output_strand1_strand2_info( outfile_prefix, strand1_strand2_info_for_each_cluster );

}


////////////////////////////////////////////////////////////////////////
void
setup_two_base_pair_pose_with_chainbreak( pose::Pose & pose,
																					Size const chainbreak_suite = 3){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace io::silent;
	using namespace id;
	using namespace protocols::swa;
	using namespace utility::io;
	using namespace chemical::rna;
	using namespace protocols::stepwise::sampling::rna;

	/////////////////////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	std::string sequence = "ccgg";
	if ( option[ seq ].user() ) sequence = option[ seq ]();
	make_pose_from_sequence( pose, sequence, *rsd_set );

	FoldTree f( 4 );
	f.new_jump( 1, 4, chainbreak_suite );
	f.set_jump_atoms( 1,
										chemical::rna::chi1_torsion_atom( pose.residue( 1 ) ),
										chemical::rna::chi1_torsion_atom( pose.residue( 4 ) )   );
	f.new_jump( 2, 3, 2 );
	f.set_jump_atoms( 2,
										chemical::rna::chi1_torsion_atom( pose.residue( 2 ) ),
										chemical::rna::chi1_torsion_atom( pose.residue( 3 ) )   );
	pose.fold_tree( f );

	add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );
	add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 3 );
	add_variant_type_to_pose_residue( pose, "CUTPOINT_LOWER", chainbreak_suite );
	add_variant_type_to_pose_residue( pose, "CUTPOINT_UPPER", chainbreak_suite+1 );

	// For now assume main chain torsions are at ideal values.
	utility::vector1< Real > strand1_torsion_set = get_suite_ideal_A_form_torsions();
	apply_suite_torsions( strand1_torsion_set, pose, 1 );
	apply_suite_torsions( strand1_torsion_set, pose, 3 );

	// For now assume delta, chi are at ideal values.
	RNA_FittedTorsionInfo rna_fitted_torsion_info;
	StepWiseRNA_BaseSugarRotamerOP base_sugar_rotamer = new StepWiseRNA_BaseSugarRotamer( ANTI, NORTH, rna_fitted_torsion_info, 20.0, 3 );
	std::cout << "NEXT ROTAMER: " << base_sugar_rotamer->get_next_rotamer() << std::endl;
	std::cout << "NEXT ROTAMER: " << base_sugar_rotamer->get_next_rotamer() << std::endl;
	for ( Size i = 1; i <= pose.total_residue(); i++ ){
		pose.set_torsion( TorsionID( i, id::BB, DELTA ), base_sugar_rotamer->delta() );
		pose.set_torsion( TorsionID( i, id::CHI, NU2 - NUM_RNA_MAINCHAIN_TORSIONS ), base_sugar_rotamer->nu2() );
		pose.set_torsion( TorsionID( i, id::CHI, NU1 - NUM_RNA_MAINCHAIN_TORSIONS ), base_sugar_rotamer->nu1() );
	}


	//pose.dump_pdb( "init.pdb" );

}

using protocols::stepwise::PoseList;

//////////////////////////////////////////////////////////////////////////
void
minimize_poses( pose::Pose & pose,
								PoseList & minimize_pose_list,
								core::io::silent::SilentFileData & minimize_silent_file_data,
								core::scoring::ScoreFunctionOP scorefxn )
{

	using namespace core::optimization;
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::chemical::rna;
	using namespace core::kinematics;
	using namespace core::io::silent;
	using namespace options;
	using namespace options::OptionKeys;

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.00000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true );

	// Define movemap.
	MoveMap mm;
	mm.set_bb( false );
	mm.set_jump( false );
	mm.set_chi( false );

	mm.set( TorsionID( 1, id::BB, EPSILON ), true );
	mm.set( TorsionID( 1, id::BB, ZETA ),    true );
	mm.set( TorsionID( 2, id::BB, ALPHA ),   true );
	mm.set( TorsionID( 2, id::BB, BETA ),    true );
	mm.set( TorsionID( 2, id::BB, GAMMA ),   true );

	mm.set( TorsionID( 3, id::BB, EPSILON ), true );
	mm.set( TorsionID( 3, id::BB, ZETA ),    true );
	mm.set( TorsionID( 4, id::BB, ALPHA ),   true );
	mm.set( TorsionID( 4, id::BB, BETA ),    true );
	mm.set( TorsionID( 4, id::BB, GAMMA ),   true );

	utility::vector1< Size > const minimize_sidechain_set = option[ minimize_sidechain_res ]();
	for (Size n = 1; n <= minimize_sidechain_set.size(); n++ ){
		std::cout << "ENABLING sidechain movement for " << minimize_sidechain_set[n] << std::endl;
		mm.set_chi( minimize_sidechain_set[n], true );
	}


	if ( option[ minimize_jump ]() ) mm.set_jump( true );

	for ( PoseList::iterator iter = minimize_pose_list.begin(); iter != minimize_pose_list.end(); iter++ ) {

		std::string const tag = iter->first;
		pose = *(iter->second);

		minimizer.run( pose, mm, *(scorefxn), options );

		BinarySilentStruct s( pose, tag );
		minimize_silent_file_data.add_structure( s );

	}

}



////////////////////////////////////////////////////////////////////////////
bool check_filter_base_stack( pose::Pose const & pose ){

	using namespace chemical::rna;
	using namespace kinematics;
	static RNA_CentroidInfo rna_centroid_info;

	// Which way is base1 pointing?
	Vector centroid1 = rna_centroid_info.get_base_centroid( pose.residue(1) );
	Stub s1 = rna_centroid_info.get_base_coordinate_system( pose.residue(1), centroid1 );

	// Which way is base2 pointing?
	Vector centroid2 = rna_centroid_info.get_base_centroid( pose.residue(2) );
	Stub s2 = rna_centroid_info.get_base_coordinate_system( pose.residue(2), centroid2 );

	if( dot( s1.M.col_z(),  centroid2 - centroid1 ) < 0.0 ) return false;
	if( dot( s2.M.col_z(),  centroid1 - centroid2 ) > 0.0 ) return false;
	return true;
}

//////////////////////////////////////////////////////////////////////////
void
assign_stack_faces( core::io::silent::SilentStructOP & s ){

	using namespace pose;
	using namespace chemical::rna;
	using namespace kinematics;
	static RNA_CentroidInfo rna_centroid_info;

	Pose pose;
	s->fill_pose( pose );

	// Which way is base1 pointing?
	Vector centroid1 = rna_centroid_info.get_base_centroid( pose.residue(1) );
	Stub s1 = rna_centroid_info.get_base_coordinate_system( pose.residue(1), centroid1 );

	// Which way is base2 pointing?
	Vector centroid2 = rna_centroid_info.get_base_centroid( pose.residue(2) );
	Stub s2 = rna_centroid_info.get_base_coordinate_system( pose.residue(2), centroid2 );

	// this is kind of coarse...
	// could actually impose a filter, e.g., on maximal z and orientation (see, e.g., RNA_LowResPotential)
	// then "zero" could mean any kind of stack is OK. This may be important when
	// modeling the unfolded state from dinucleotides.
	int stack_face1 = ( dot( s1.M.col_z(),  centroid2 - centroid1 )  > 0 ?   1 : -1 );
	int stack_face2 = ( dot( s2.M.col_z(),  centroid1 - centroid2 )  > 0 ?   1 : -1 );

	s->add_energy( "stackface1", stack_face1 );
	s->add_energy( "stackface2", stack_face2 );

}

//////////////////////////////////////////////////////////////////////////
void
apply_filter_base_stack_direction( protocols::stepwise::PoseList & minimize_pose_list ){

	using namespace protocols::swa;
	using namespace pose;


	PoseList new_pose_list;
	for ( PoseList::iterator iter = minimize_pose_list.begin(); iter != minimize_pose_list.end(); iter++ ){
		Pose & pose = *(iter->second);
		if (!check_filter_base_stack( pose ) ) continue;
		new_pose_list[ iter->first ] = iter->second;
	}

	minimize_pose_list = new_pose_list;
}

//////////////////////////////////////////////////////////////////////////
void
apply_filter_base_stack_direction( core::io::silent::SilentFileDataOP &  sfd ){

	using namespace chemical::rna;
	using namespace protocols::swa;
	using namespace pose;
	using namespace kinematics;
	using namespace core::io::silent;

	SilentFileDataOP sfd_new = new SilentFileData;
	for ( core::io::silent::SilentFileData::iterator iter = sfd->begin(),
					end = sfd->end(); iter != end; ++iter ) {
		Pose pose;
		iter->fill_pose( pose );
		if (!check_filter_base_stack( pose ) ) continue;
		sfd_new->add_structure( *iter );
	}

	std::cout << "AFTER STACK DIRECTION FILTER, retained " << sfd_new->size() << " of " << sfd->size() << " decoys" << std::endl;

	sfd = sfd_new;

}

//////////////////////////////////////////////////////////////////////////
void
base_pair_to_base_pair_test(){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace io::silent;
	using namespace protocols::swa;
	using namespace protocols::stepwise::sampling::rna;
	using namespace utility::io;

	/////////////////////////////////////////////
	// set up starting pose.
	/////////////////////////////////////////////
	Pose pose;
	setup_two_base_pair_pose_with_chainbreak( pose );
	Pose start_pose = pose;
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	std::string silent_file = "";
	if ( option[ out::file::silent ].user() )	silent_file = option[ out::file::silent ]();

	////////////////////////////////////////////////////////////////////
	utility::vector1< InputStreamWithResidueInfoOP > input_streams;

	protocols::stepwise::initialize_input_streams( input_streams );
	ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	input_streams[1]->set_rsd_set( rsd_set );
	input_streams[2]->set_rsd_set( rsd_set );

	StepWisePoseSampleGeneratorOP sample_generator = new StepWisePoseCombineSampleGenerator( input_streams );
	Size count( 0 );

	Size moving_suite( 1 ), chainbreak_suite( 3 );
	RNA_LoopCloseSampler rna_loop_close_sampler( moving_suite, chainbreak_suite );
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	if ( option[ score::weights ].user() ) scorefxn = get_score_function();
	// This is intrusive, but needs to be done!

	ScoreFunctionOP minimize_scorefxn = new ScoreFunction;
	*minimize_scorefxn = *scorefxn;
	minimize_scorefxn->set_weight( linear_chainbreak, 5.0 );

	rna_loop_close_sampler.set_scorefxn( scorefxn );
	rna_loop_close_sampler.set_bin_size( get_bin_size() );
	rna_loop_close_sampler.set_rep_cutoff ( 3 * option[ rep_cutoff ]() );
	//	rna_loop_close_sampler.set_silent_file( "raw_"+silent_file );


	while( sample_generator->has_another_sample() ){

		sample_generator->get_next_sample( pose );

		count++;
		std::string const tag = "S_"+lead_zero_string_of( count,5 );
		//pose.dump_pdb( tag+".pdb" );

		//sample epsilon,zeta,alpha,beta,gamma of two suites.
		rna_loop_close_sampler.apply( pose );

		// did we get anything out?
		SilentFileDataOP silent_file_data = rna_loop_close_sampler.silent_file_data();

		if ( silent_file_data->size() > 0 ){

			//			if ( option[ filter_base_stack_direction ]() ) {
			//				apply_filter_base_stack_direction( silent_file_data );
			//			}

			//cluster
			protocols::stepwise::StepWiseLegacyClusterer stepwise_clusterer( silent_file_data );
			Size max_decoys( 1000 );
			if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
			stepwise_clusterer.set_max_decoys( max_decoys );
			stepwise_clusterer.set_cluster_by_all_atom_rmsd( true );
			stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );
			Real cluster_radius( 0.25 );
			if ( option[ OptionKeys::cluster::radius ].user() ) cluster_radius = option[ OptionKeys::cluster::radius ]();
			stepwise_clusterer.set_cluster_radius( cluster_radius	);

			stepwise_clusterer.cluster();
			//		stepwise_clusterer.output_silent_file( silent_file );

			//minimize. Probably need to specify movemap explicitly!
			PoseList minimize_pose_list = stepwise_clusterer.clustered_pose_list();

			//if ( option[ filter_base_stack_direction ]() ) apply_filter_base_stack_direction( minimize_pose_list );

			SilentFileDataOP minimize_silent_file_data = new SilentFileData;
			minimize_poses( pose, minimize_pose_list, *minimize_silent_file_data, minimize_scorefxn );

			//			if ( minimize_silent_file_data->size() == 0 ) { //need a garbage pose -- must output something.
			//				(*minimize_scorefxn)( start_pose );
			//				BinarySilentStruct s( start_pose, "GARBAGE" );
			//				minimize_silent_file_data->add_structure( s );
			//			}

			//cluster again. This isn't really necessary.
			//save lowest energy state --> note that we need some kind of tag of: parent outfiles, parent tags.
			stepwise_clusterer.set_silent_file_data( minimize_silent_file_data );
			stepwise_clusterer.cluster();

			// save lowest energy pose.
			SilentStructOP best_model = stepwise_clusterer.silent_struct_output_list()[1];
			best_model->set_decoy_tag( tag );
			assign_stack_faces( best_model );
			silent_file_data->write_silent_struct( *best_model, silent_file, false /*write score only*/ );

			if (option[ output_all ]() ) stepwise_clusterer.output_silent_file( "ALL_"+silent_file );

		} else {
			// this is just a 'garbage' pose.
			std::cout << "GARBAGE POSE!!" << std::endl;
			(*minimize_scorefxn)( pose );
			SilentStructOP s = new BinarySilentStruct( pose, tag );
			assign_stack_faces( s );
			silent_file_data->write_silent_struct( *s, silent_file, false );
		}

	}

}

///////////////////////////////////////////////////////////////
void
rna_close_loop_test(){

	using namespace pose;
	using namespace protocols::swa;
	using namespace protocols::stepwise::sampling::rna;
	using namespace options;
	using namespace options::OptionKeys;

	// set up two-base pair pose.
	Pose pose;
	setup_two_base_pair_pose_with_chainbreak( pose );

	//Need to setup starting base pair with user-input rigid body setting.
	utility::vector1< utility::vector1< Real > > reference_rigid_body_settings_fixed_pair, reference_rigid_body_settings_new_pair;
	read_rigid_body_settings( option[ reference_rigid_body_samples_fixed_pair ](), reference_rigid_body_settings_fixed_pair );
	read_rigid_body_settings( option[ reference_rigid_body_samples_new_pair ](), reference_rigid_body_settings_new_pair );
	//apply rigid body settings to top and bottom
	utility::vector1< Real > const & rbs_fixed_pair = reference_rigid_body_settings_fixed_pair[ option[ fixed_pair_state_number ]() ];
	utility::vector1< Real > const & rbs_new_pair   = reference_rigid_body_settings_new_pair[ option[ new_pair_state_number ]() ];

	apply_rigid_body_settings( pose, rbs_fixed_pair, 1, 4 );
	apply_rigid_body_settings( pose, rbs_new_pair, 2, 3 );

	/////////////////////////////////////////////////////////////////
	// NOTE!!!!! Following no longer works.
	RNA_AnalyticLoopCloser rna_analytic_loop_closer( 1, 3 );
	rna_analytic_loop_closer.apply( pose );
	pose.dump_pdb( "closed.pdb" );

}

///////////////////////////////////////////////////////////////
void
sample_state_to_state(
											pose::Pose & pose,
											Size const moving_suite,
											Size const chainbreak_suite,
											utility::vector1< Real > const & rbs_new_pair /*is this necessary?*/,
											utility::vector1< utility::vector1< Real > > & all_torsion_info ){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace pose;
	using namespace scoring;
	using namespace io::silent;
	using namespace chemical;
	using namespace id;
	using namespace chemical::rna;
	using namespace protocols::swa;
	using namespace protocols::stepwise::sampling::rna;

	RNA_LoopCloseSampler rna_loop_close_sampler( moving_suite, chainbreak_suite );

	/////////////////////////////////////////////////
	// Basic setup for clash checks and scoring.
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	if ( option[ score::weights ].user() ) scorefxn = get_score_function();
	rna_loop_close_sampler.set_scorefxn( scorefxn );

	/////////////////////////////////////////////////
	std::string silent_file = "";
	if ( option[ out::file::silent ].user() )	silent_file = option[ out::file::silent ]();
	rna_loop_close_sampler.set_silent_file( silent_file );

	/////////////////////////////////////////////////
	PoseOP native_pose;
	if ( option[ in::file::native ].user() )	{
		native_pose = new Pose;
		ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
		io::pdb::pose_from_pdb( *native_pose, *rsd_set, option[ in::file::native ]() );
		rna_loop_close_sampler.set_native_pose( native_pose );
	}
	/////////////////////////////////////////////////

	rna_loop_close_sampler.set_bin_size( get_bin_size() );
	rna_loop_close_sampler.set_center_around_native( option[ center_around_native ]() );
	rna_loop_close_sampler.set_save_torsion_info( true );
	rna_loop_close_sampler.set_rbs_new_pair( rbs_new_pair );
	rna_loop_close_sampler.set_rep_cutoff ( 3 * option[ rep_cutoff ]() );
	rna_loop_close_sampler.set_torsion_range( option[ torsion_range ]() );
	rna_loop_close_sampler.set_torsion_increment( option[ torsion_increment ]() );
	rna_loop_close_sampler.set_just_output_score( option[ just_output_score ]() );

	rna_loop_close_sampler.apply( pose );

	all_torsion_info = rna_loop_close_sampler.all_torsion_info();

}



///////////////////////////////////////////////////////////////
void
two_base_pairs_via_loop_closure_test(
																		 utility::vector1< utility::vector1< Real > > const & reference_rigid_body_settings_fixed_pair,
																		 utility::vector1< utility::vector1< Real > > const & reference_rigid_body_settings_new_pair,
																		 utility::vector1< Size > const & which_rigid_body_settings_fixed_pair,
																		 utility::vector1< Size > const & which_rigid_body_settings_new_pair) {

	using namespace options;
	using namespace options::OptionKeys;
	using namespace pose;
	using namespace protocols::swa;
	using namespace protocols::stepwise::sampling::rna;
	using namespace utility::io;

	// set up two-base pair pose.
	Pose pose;
	Size moving_suite( 1 ), chainbreak_suite( 3 );
	if ( option[ switch_chainbreak ]() ) { //This option does not work.
		moving_suite = 3;
		chainbreak_suite = 1;
	}
	setup_two_base_pair_pose_with_chainbreak( pose, chainbreak_suite );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	std::string const outfile = option[ out::file::o ]();

	utility::vector1< std::pair<Size,Size> > which_state;
	utility::vector1< utility::vector1< Real > > all_torsion_info;

	for ( Size m = 1; m <= which_rigid_body_settings_fixed_pair.size(); m++ ){

		Size const & fixed_pair_setting = which_rigid_body_settings_fixed_pair[ m ];
		utility::vector1< Real > rbs_fixed_pair = reference_rigid_body_settings_fixed_pair[ fixed_pair_setting ];
		std::cout << "Fixed pair --> state: " << fixed_pair_setting << std::endl;

		for ( Size n = 1; n <= which_rigid_body_settings_new_pair.size(); n++ ){

			Size const & new_pair_setting = which_rigid_body_settings_new_pair[ n ];
			utility::vector1< Real > rbs_new_pair = reference_rigid_body_settings_new_pair[ which_rigid_body_settings_new_pair[ n ]  ];
			std::cout << "New pair --> state: " << new_pair_setting << std::endl;

			//apply rigid body settings to top and bottom
			apply_rigid_body_settings( pose, rbs_fixed_pair, 1, 4 );
			apply_rigid_body_settings( pose, rbs_new_pair, 2, 3 );

			utility::vector1< utility::vector1< Real > > torsion_info;

			sample_state_to_state( pose,
														 moving_suite, chainbreak_suite,
														 rbs_new_pair,
														 torsion_info );

			for ( Size i = 1; i <= torsion_info.size(); i++ ){
				which_state.push_back( std::make_pair( fixed_pair_setting, new_pair_setting ) );
				all_torsion_info.push_back( torsion_info[ i ] );
			}

		}
	}

	ozstream out( outfile );
	Size const MAX_TRIES = 100;
	for ( Size n = 1; n <= MAX_TRIES; n++ ){
		if ( out.good() ) break;
		sleep( 5 );
		out.open( outfile );
	}
	for ( Size i = 1; i <= all_torsion_info.size(); i++ ){
		out << which_state[i].first << ' ' << which_state[i].second;
		utility::vector1< Real > const & info = all_torsion_info[i];
		for ( Size j = 1; j <= info.size(); j++ ) out << ' ' << info[j];
		out << std::endl;
	}
	out.close();



}


///////////////////////////////////////////////////////////////
void
two_base_pairs_via_loop_closure_test(){

	using namespace options;
	using namespace options::OptionKeys;

	utility::vector1< utility::vector1< Real > > reference_rigid_body_settings_fixed_pair, reference_rigid_body_settings_new_pair;
	read_rigid_body_settings( option[ reference_rigid_body_samples_fixed_pair ](), reference_rigid_body_settings_fixed_pair );
	read_rigid_body_settings( option[ reference_rigid_body_samples_new_pair ](), reference_rigid_body_settings_new_pair );

	utility::vector1< Size > which_rigid_body_settings_fixed_pair, which_rigid_body_settings_new_pair;

	if ( !option[ fixed_pair_state_number ].user()  ) utility_exit_with_message( "Must supply -fixed_pair_state_number" );
	which_rigid_body_settings_fixed_pair.push_back(  option[ fixed_pair_state_number ]() );

	if ( option[ all_new_pair_states ] ) {
		for ( Size n = 1; n <= reference_rigid_body_settings_new_pair.size(); n++ ) which_rigid_body_settings_new_pair.push_back( n );
	} else if ( option[ num_new_pair_states ].user() ) {
		for ( int n = 1; n <= option[ num_new_pair_states ](); n++ ) which_rigid_body_settings_new_pair.push_back( n );
	} else{
		if ( !option[ new_pair_state_number ]()  ) utility_exit_with_message( "Must supply -new_pair_state_number or -all_new_pair_states" );
		which_rigid_body_settings_new_pair.push_back(  option[ new_pair_state_number ]() );
	}

	two_base_pairs_via_loop_closure_test( reference_rigid_body_settings_fixed_pair,
																				reference_rigid_body_settings_new_pair,
																				which_rigid_body_settings_fixed_pair,
																				which_rigid_body_settings_new_pair);

}


///////////////////////////////////////////////////////////////////////////////
void
check_determinant_test(){

	//	using namespace arma;

	// mat m;
	// m.zeros( 6, 6 );
	// for ( Size i = 0; i < 6; i++ ) m( i, i ) = 1.0;
	// m( 3, 4) = -0.55;
	// m( 4, 3) = -0.55;
	// m( 1, 5) =  0.55;
	// m( 5, 1) =  0.55;
	// m( 2, 4) =  0.20;
	// m( 4, 2) = -0.20;

	// utility::vector1< utility::vector1< Real > > m_rosetta;
	// for ( Size i = 0; i < 6; i++ ) {
	// 	utility::vector1< Real > v;
	// 	for ( Size j = 0; j < 6; j++ ) {
	// 		std::cout << ' ' << m( i, j );
	// 		v.push_back( m( i, j ) );
	// 	}
	// 	std::cout << std::endl;
	// 	m_rosetta.push_back( v );
	// }

	// std::cout << "TAKING DETERMINANT: " << det( m ) << std::endl;

	// std::cout << "BRUTE FORCE       : " << get_determinant( m_rosetta ) << std::endl;

}

utility::vector1< Real >
reverse( utility::vector1< Real > const rbs ){

	using namespace protocols::swa;

	utility::vector1< Real > rbs_reverse;
	Real const & alpha = rbs[1];
	Real const & beta  = rbs[2];
	Real const & gamma = rbs[3];

	//Real const alpha_new = -rbs[3];
	//	Real const beta_new  = -rbs[2];
	//	Real const gamma_new = -rbs[1];
	Real const alpha_new = 180.0-rbs[3];
	Real const beta_new  = rbs[2];
	Real const gamma_new = 180.0-rbs[1];

	Matrix M, M_new;
	create_euler_rotation( M, alpha, beta, gamma );
	create_euler_rotation( M_new, alpha_new, beta_new, gamma_new );

	// Matrix M_check = M * M_new;
	// for ( Size i = 1; i <= 3; i++ ){
	// 	for ( Size j = 1; j <= 3; j++ ){
	// 		std::cout << ' ' << M_check( i, j );
	// 	}
	// 	std::cout << std::endl;
	// }
	// std::cout << std::endl;

	Vector v_reverse = -1.0 * M.transposed() * Vector( rbs[4], rbs[5], rbs[6] );

	rbs_reverse.push_back( alpha_new );
	rbs_reverse.push_back( beta_new );
	rbs_reverse.push_back( gamma_new );
	rbs_reverse.push_back( v_reverse(1) );
	rbs_reverse.push_back( v_reverse(2) );
	rbs_reverse.push_back( v_reverse(3) );
	for ( Size n = 7; n <= rbs.size(); n++ ) {
		rbs_reverse.push_back( rbs[ n ] );
	}

	return rbs_reverse;

}

///////////////////////////////////////////////////////////////////////////////
void
reverse_rbs_test(){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace utility::io;

	utility::vector1< utility::vector1< Real > > rigid_body_settings, rigid_body_settings_reverse;
	read_rigid_body_settings( option[ rigid_body_samples ](), rigid_body_settings );

	for ( Size n = 1; n <= rigid_body_settings.size(); n++ ) {
		rigid_body_settings_reverse.push_back( reverse( rigid_body_settings[ n ] ) );
	}

	std::string const outfile = option[ out::file::o ]();
	ozstream out( outfile );

	for ( Size n = 1; n <= rigid_body_settings.size(); n++ ) {
		utility::vector1< Real > const & rbs = rigid_body_settings_reverse[ n ];
		for ( Size n = 1; n <= 8; n++ )	out << ' ' << rbs[ n ];
		out << std::endl;
	}

	out.close();

}


///////////////////////////////////////////////////////////////
void
setup_dinucleotide_pose( pose::Pose & pose ){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace kinematics;
	using namespace pose;
	using namespace io::silent;
	using namespace id;
	using namespace protocols::swa;
	using namespace protocols::stepwise::sampling::rna;
	using namespace chemical::rna;

	ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	std::string sequence = "cc";
	if ( option[ seq ].user() ) sequence = option[ seq ]();
	make_pose_from_sequence( pose, sequence, *rsd_set );

	add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );

	// For now assume delta, chi are at ideal values.
	RNA_FittedTorsionInfo rna_fitted_torsion_info;
	StepWiseRNA_BaseSugarRotamerOP base_sugar_rotamer = new StepWiseRNA_BaseSugarRotamer( ANTI, NORTH, rna_fitted_torsion_info, 20, 3 );

	// This is kind of ugly. Need to know that the second BaseSugarRotamer  is the "ideal" one.
	// Anyway, stick with this for now, then consult with Parin on better design.
	std::cout << "NEXT ROTAMER: " << base_sugar_rotamer->get_next_rotamer() << std::endl;
	std::cout << "NEXT ROTAMER: " << base_sugar_rotamer->get_next_rotamer() << std::endl;

	for ( Size i = 1; i <= pose.total_residue(); i++ ){
		pose.set_torsion( TorsionID( i, id::BB, DELTA ), base_sugar_rotamer->delta() );
		pose.set_torsion( TorsionID( i, id::CHI, NU2 - NUM_RNA_MAINCHAIN_TORSIONS ), base_sugar_rotamer->nu2() );
		pose.set_torsion( TorsionID( i, id::CHI, NU1 - NUM_RNA_MAINCHAIN_TORSIONS ), base_sugar_rotamer->nu1() );
	}

	utility::vector1< Real > strand1_torsion_set = get_suite_ideal_A_form_torsions();
	apply_suite_torsions( strand1_torsion_set, pose, 1 );

	apply_south_syn_to_dinucleotide_pose( pose );


}


///////////////////////////////////////////////////////////////
void
output_torsion_list( std::string const & outfile,
										 utility::vector1< utility::vector1< Real > > const & torsion_info,
										 utility::vector1< Real > const & scores,
										 Real const bin_size ){

	using namespace utility::io;

	//	std::string const outfile = outfile_prefix + ".cluster" + ObjexxFCL::lead_zero_string_of( n, 5 );
	ozstream out( outfile );

	Real bin_size_in_radians = radians( bin_size );
	Real volume_element = 1.0;
	for( Size i = 1; i <= 5; i++ )  volume_element *= bin_size_in_radians;

	for ( Size n = 1; n <= torsion_info.size(); n++ ){

		utility::vector1< Real > const & info = torsion_info[ n ];

		for ( Size i = 1; i <= info.size(); i++ ){
			out << ' ' << info[i];
		}

		out << ' ' << scores[n];
		out << ' ' << volume_element << std::endl;
	}

	out.close();


}


///////////////////////////////////////////////////////////////
void
dinucleotide_test(){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace pack;
	using namespace pack::task;
	using namespace io::silent;
	using namespace id;
	using namespace protocols::swa;
	using namespace utility::io;
	using namespace chemical::rna;
	using namespace protocols::stepwise::sampling::rna;

	std::string const outfile_prefix = option[ out::file::o ]();

	/////////////////////////////////////////////
	// set up starting pose.
	/////////////////////////////////////////////
	Pose pose;
	setup_dinucleotide_pose( pose );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );
	Pose ideal_pose = pose;
	ideal_pose.dump_pdb( "ideal.pdb" );

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	if ( option[ score::weights ].user() ) scorefxn = get_score_function();
	ScoreFunctionOP rep_scorefxn = new ScoreFunction;
	rep_scorefxn->set_weight( fa_rep, 0.12 );

	std::string silent_file = "";
	if ( option[ out::file::silent ].user() )	silent_file = option[ out::file::silent ]();
	bool const just_score = option[ just_output_score ]();

	// Need to iterate through suite conformations for residue 1.
	utility::vector1< utility::vector1< Real > > torsion_list1;

	Size const moving_suite = 1;
	utility::vector1< Size > suite_res_list = make_vector1( moving_suite );

	Real const bin_size = get_bin_size();
	StepWiseRNA_RotamerGeneratorWrapperOP rotamer_generator = new StepWiseRNA_RotamerGeneratorWrapper( pose,
																																																			 suite_res_list,
																																																			 false /*sample_sugar_and_base1*/,
																																																			 false /*sample_sugar_and_base2*/,
																																																			 bin_size );
	Real const fa_rep_score_baseline = initialize_fa_rep( pose, make_vector1( moving_suite ), rep_scorefxn );
	Real const rep_cutoff_ = option[ rep_cutoff ]();

	// initialize for o2prime rotamer trials. Probably should stuff this into an initialization function.
	PackerTaskOP o2prime_pack_task =  pack::task::TaskFactory::create_packer_task( pose );
	for (Size i = 1; i <= pose.total_residue(); i++) {
		o2prime_pack_task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
		o2prime_pack_task->nonconst_residue_task(i).or_ex4( true ); //extra rotamers?? Parin S. Jan 28, 2010
		o2prime_pack_task->nonconst_residue_task(i).or_include_current( true );
	}
	ScoreFunctionOP o2prime_pack_scorefxn = new ScoreFunction;
	// Each of the following terms have been pretty optimized for the packer (trie, etc.)
	o2prime_pack_scorefxn->set_weight( fa_atr, scorefxn->get_weight( fa_atr ) );
	o2prime_pack_scorefxn->set_weight( fa_rep, scorefxn->get_weight( fa_rep ) );
	o2prime_pack_scorefxn->set_weight( hbond_lr_bb_sc, scorefxn->get_weight( hbond_lr_bb_sc ) );
	o2prime_pack_scorefxn->set_weight( hbond_sr_bb_sc, scorefxn->get_weight( hbond_sr_bb_sc ) );
	o2prime_pack_scorefxn->set_weight( hbond_sc, scorefxn->get_weight( hbond_sc ) );
	o2prime_pack_scorefxn->set_energy_method_options( scorefxn->energy_method_options() );
	// note that geom_sol is not optimized well --> replace with lk_sol for now.
	o2prime_pack_scorefxn->set_weight( fa_sol, scorefxn->get_weight( lk_nonpolar ) );


	// Get ready for main loop.
	Size count( 0 );
	Size count_saved( 0 );
	clock_t const time_start( clock() );

	utility::vector1< utility::vector1< Real > > torsion_list;
	utility::vector1< Real > score_list;

	while( rotamer_generator->has_another_rotamer() ){

		count++;

		utility::vector1< Torsion_Info > current_rotamer = rotamer_generator->get_next_rotamer();
		apply_rotamer( pose, current_rotamer);

		// Disallow steric clashes.
		if ( !check_clash( pose, fa_rep_score_baseline, rep_cutoff_, rep_scorefxn  ) ) continue;

		if ( option[ o2prime_trials ]() ) pack::rotamer_trials( pose, *o2prime_pack_scorefxn, o2prime_pack_task );

		Real const score = (*scorefxn)( pose );
		save_torsions( pose, moving_suite, torsion_list );
		score_list.push_back( score );

		// Save torsion angles, centroid, and base coordinate system.
		count_saved++;

		if ( silent_file.size() > 0  ){

			std::string const tag = "S_" + ObjexxFCL::lead_zero_string_of( count, 6 );
			BinarySilentStruct s( pose, tag ); // this is RNA-centric -- could make it OK for proteins.
			Real const rmsd_to_ideal = all_atom_rmsd( pose, ideal_pose );
			s.add_energy( "all_rms", rmsd_to_ideal );
			s.add_energy( "log_vol", log( radians(bin_size) * radians(bin_size) * radians(bin_size) * radians(bin_size) * radians(bin_size) ) );

			static SilentFileData sfd;
			sfd.write_silent_struct( s, silent_file, just_score );
		}

	}
	std::cout << "Total number of rotamers applied: " << count << std::endl;
	std::cout << "Total number that passed cuts:    " << count_saved << std::endl;

	Real const time_in_dinucleotide_test = static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC;
	std::cout << "Time in dinucleotide sampler: " <<  time_in_dinucleotide_test << std::endl;

	output_torsion_list( option[ out::file::o ], torsion_list, score_list, bin_size );

}

///////////////////////////////////////////////////////////////
void
delta_chi_correction_test(){

	using namespace options;
	using namespace options::OptionKeys;
	using namespace chemical;
	using namespace scoring;
	using namespace chemical::rna;
	using namespace kinematics;
	using namespace optimization;
	using namespace pose;
	using namespace io::silent;
	using namespace id;
	using namespace protocols::swa;
	using namespace utility::io;
	using namespace chemical::rna;
	using namespace protocols::stepwise::sampling::rna;

	ResidueTypeSetCAP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_HIRES_WTS );
	if ( option[ score::weights ].user() ) scorefxn = get_score_function();

	utility::vector1< std::string >  nts;
	nts.push_back("a");
	nts.push_back("c");
	nts.push_back("g");
	nts.push_back("u");

	std::string const outfile = option[ out::file::o ];
	ozstream out( outfile );

	for ( Size n = 1; n <= nts.size(); n++ ){

		Pose pose;

		std::string const & seq = nts[n];
		make_pose_from_sequence( pose, seq,	*rsd_set );

		std::cout << "--------------------------------------" << std::endl;
		std::cout << "--------------------------------------" << std::endl;
		std::cout << "Doing nucleotide: " << seq << std::endl;
		std::cout << "--------------------------------------" << std::endl;
		std::cout << "--------------------------------------" << std::endl;

		apply_ideal_A_form_torsions( pose );

		std::cout << "------- STARTING CONFORMATION -------" << std::endl;
		(*scorefxn)( pose );
		scorefxn->show( std::cout, pose );

		add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );

		std::cout << "------- VIRTUAL_PHOSPHATE ADDED -------" << std::endl;
		Real const start_score = (*scorefxn)( pose );
		scorefxn->show( std::cout, pose );

		remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );
		add_variant_type_to_pose_residue( pose, "VIRTUAL_BACKBONE_EXCEPT_C1PRIME", 1 );

		std::cout << "------- VIRTUAL_BACKBONE_EXCEPT_C1PRIME ADDED -------" << std::endl;
		Real const end_score = (*scorefxn)( pose );
		scorefxn->show( std::cout, pose );

		Real const diff_score = start_score - end_score;
		out << n
			<< ' ' << pose.torsion( TorsionID( 1, id::BB, DELTA ) )
			<< ' ' << pose.torsion( TorsionID( 1, id::CHI, chemical::rna::CHI - NUM_RNA_MAINCHAIN_TORSIONS) )
			<< ' ' << start_score << ' ' << end_score << ' ' << diff_score << std::endl;

		remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE_AND_", 1 );
		add_variant_type_to_pose_residue( pose, "VIRTUAL_BACKBONE_EXCEPT_C1PRIME", 1 );
		std::cout << "------- VIRTUAL_BACKBONE_EXCEPT_C1PRIME ADDED -------" << std::endl;
		(*scorefxn)( pose );
		scorefxn->show( std::cout, pose );

		remove_variant_type_from_pose_residue( pose, "VIRTUAL_BACKBONE_EXCEPT_C1PRIME", 1 );
		add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESIDUE", 1 );

		std::cout << "------- VIRTUAL RNA RESIDUE -------" << std::endl;
		(*scorefxn)( pose );
		scorefxn->show( std::cout, pose );

	}

	std::cout << "Output results to: " << outfile << std::endl;
	out.close();

}

/////////////////////////////////////////////////////////////////////
void
reverse_doublet_test(){

	using namespace core::io::silent;
	using namespace core::kinematics;
	using namespace core::pose;
	using namespace options;
	using namespace options::OptionKeys;


	SilentFileData silent_file_data;
	silent_file_data.read_file( option[ in::file::silent ]()[ 1 ] );

	std::string silent_file_out = option[ out::file::silent ]();
	Pose pose;

	for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(),
					end = silent_file_data.end(); iter != end; ++iter ) {

		std::string tag = iter->decoy_tag();
		iter->fill_pose( pose );

		// reversal
		FoldTree f( pose.fold_tree() ), f_new( f );
		Pose new_pose;
		new_pose.append_residue_by_bond( pose.residue( 2 ) );
		new_pose.append_residue_by_jump( pose.residue( 1 ), 1  );
		f_new.set_jump_atoms( 1, f.downstream_atom(1),f.upstream_atom(1) );
		new_pose.fold_tree( f_new );

		if ( tag == "S_0" ) {
			pose.dump_pdb( "START.pdb");
			new_pose.dump_pdb( "REVERSE.pdb");
		}

		BinarySilentStruct s( new_pose, tag ); // will this copy in scores?
		s.silent_energies(  iter->get_silent_energies() );

		silent_file_data.write_silent_struct( s, silent_file_out, false /*just score*/ );

	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace options;

	if ( option[ base_doublet_rmsd ]() ){
		base_doublet_rmsd_test();
	} else if ( option[ cluster_rigid_body_settings ]() ){
		cluster_rigid_body_settings_test();
	} else if ( option[ assign_to_clusters ]() ){
		assign_rigid_body_settings_to_clusters_test();
	} else if ( option[ two_base_pairs ]() ){
		two_base_pairs_test();
	} else if ( option[ two_base_pairs_via_loop_closure ]() ){
		two_base_pairs_via_loop_closure_test();
	} else if ( option[ close_loop_test ]() ){
	  rna_close_loop_test();
	} else if ( option[ check_determinant ]() ){
	  check_determinant_test();
	} else if ( option[ reverse_rbs ]() ){
	  reverse_rbs_test();
	} else if ( option[ dinucleotide ]() ){
	  dinucleotide_test();
	} else if ( option[ delta_chi_correction ]() ){
	  delta_chi_correction_test();
	} else if ( option[ finely_sample_base_pair ]() ){
	  finely_sample_base_pair_test();
	} else if ( option[ base_pair_to_base_pair ]() ){
	  base_pair_to_base_pair_test();
	} else if ( option[ reverse_doublet ]() ){
	  reverse_doublet_test();
	} else {
		define_states_test();
	}

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace options;

	utility::vector1< Size > blank_size_vector;

	NEW_OPT( n_sample, "number of samples per torsion angle", 18 );
	NEW_OPT( n_sample_beta, "number of samples in tilt angle beta", 18 );
	NEW_OPT( xyz_sample, "spacing in xyz search, in Angstroms", 1.0 );
	NEW_OPT( box_radius, "spacing in xyz search, in Angstroms", 10.0 );
	NEW_OPT( score_cutoff, "Scoring cutoff", 10.0 );
	NEW_OPT( temperature, "Temperature", 3.0 );
	NEW_OPT( contact_cutoff, "how close atoms need to be to define contact", 4.5 );
	NEW_OPT( steric_dist_cutoff, "how close heavy atoms need to be to define clash", 2.5 );
	NEW_OPT( bin_sample, "Torsion bin size (in degrees)", 20 );
	NEW_OPT( RBangle_range, "rigid body angle range for alpha, gamma (in degrees)", 20.0 );
	NEW_OPT( RBangle_increment, "rigid body angle increment for alpha, gamma (in degrees)", 2.0 );
	NEW_OPT( torsion_range, "Torsion range for alpha, gamma (in degrees)", 20.0 );
	NEW_OPT( torsion_increment, "Torsion increment for alpha, gamma (in degrees)", 2.0 );
	NEW_OPT( min_contacts, "minimum number of contacts", 1 );
	NEW_OPT( min_hbonds, "minimum number of H-bonds", 0 );
	NEW_OPT( fa_rep_cutoff, "fa rep cutoff for rigid bdy sampling", 0.0 );
	NEW_OPT( fixed_pair_state_number, "from reference file, which rigid body setting to use", 1 );
	NEW_OPT( new_pair_state_number, "from reference file, which rigid body setting to use", 1 );
	NEW_OPT( only_positive_Z, "only allow positive contributions to partition function", false );
	NEW_OPT( o2prime_trials, "in dinucleotide test, do rotamer trials", false );
	NEW_OPT( cycle_axes, "different coordinate system", false );
	NEW_OPT( do_not_rotate_base2, "consistency test -- leave base 2 in arbitrary rotation", false );
	NEW_OPT( cluster_poses, "Cluster after sampling", false );
	NEW_OPT( filter_base_stack_direction, "filter base stack directionality in base_pair_to_base_pair", false );
	NEW_OPT( expand_chi, "Expand chi after state definition", false );
	NEW_OPT( base_doublet_rmsd, "Testing base doublet rmsd", false );
	NEW_OPT( cluster_rigid_body_settings, "Cluster list of rigid body settings", false );
	NEW_OPT( two_base_pairs, "Build next base pair on existing base pair", false );
	NEW_OPT( two_base_pairs_via_loop_closure, "Build next base pair on existing base pair, using loop closure for speed", false );
	NEW_OPT( dinucleotide, "Build two bases", false );
	NEW_OPT( test_ideal, "Build next base pair on existing base pair", false );
	NEW_OPT( fine_torsions, "Fine torsions", false );
	NEW_OPT( super_fine_torsions, "Super fine torsions", false );
	NEW_OPT( gzip_out, "Gzip outfile", false );
	NEW_OPT( close_loop_test, "RNA loop close test", false );
	NEW_OPT( check_determinant, "Check determinant test", false );
	NEW_OPT( assign_to_clusters, "Assign list of rigid body settings to clusters defined by another list of rigid body settings", false );
	NEW_OPT( all_new_pair_states, "In two base pair run, iterate over all states of new pair", false );
	NEW_OPT( num_new_pair_states, "In two base pair run, iterate over all states of new pair", 0 );
	NEW_OPT( reverse_rbs, "Take a list of rigid body settings and reverse it", false );
	NEW_OPT( switch_chainbreak, "switch chainbreak", false );
	NEW_OPT( delta_chi_correction, "delta_chi_correction", false );
	NEW_OPT( finely_sample_base_pair, "finely sample base pair", false );
	NEW_OPT( base_pair_to_base_pair, "connect base pair to next base pair", false );
	NEW_OPT( just_output_score, "just output score", false );
	NEW_OPT( force_antiparallel_bases, "force antiparallel bases", false );
	NEW_OPT( force_parallel_bases, "force parallel bases", false );
	NEW_OPT( center_around_native, "center around native", false );
	NEW_OPT( ignore_o2prime_hbonds_in_filter, "Ignore O2' hbonds in cutoff of rigid body sampler", false );
	NEW_OPT( assign_WC_edges, "Ignore O2' hbonds in cutoff of rigid body sampler", false );
	NEW_OPT( virtualize_phosphate,   "virtualize phosphate instead of backbone", false );
	NEW_OPT( superimpose_over_all_res,   "during clustering, calculate rms over all residues", false );
	NEW_OPT( south1,   "", false );
	NEW_OPT( south2,   "", false );
	NEW_OPT( syn_chi1, "", false );
	NEW_OPT( syn_chi2, "", false );
	NEW_OPT( chi1, "", 0.0 );
	NEW_OPT( chi2, "", 0.0 );
	NEW_OPT( rigid_body_samples, "input file with alpha, beta, gamma, x, y, and z", "" );
	NEW_OPT( reference_rigid_body_samples, "input file with alpha, beta, gamma, x, y, and z", "rbs_cluster.txt" );
	NEW_OPT( reference_rigid_body_samples_fixed_pair, "input file with alpha, beta, gamma, x, y, and z", "rbs_cluster.txt" );
	NEW_OPT( reference_rigid_body_samples_new_pair  , "input file with alpha, beta, gamma, x, y, and z", "rbs_cluster.txt" );
	NEW_OPT( x_min, "input parameter", 0.0 );
	NEW_OPT( x_max, "input parameter", 0.0 );
	NEW_OPT( x_increment, "input parameter", 1.0 );
	NEW_OPT( y_min, "input parameter", 0.0 );
	NEW_OPT( y_max, "input parameter", 0.0 );
	NEW_OPT( y_increment, "input parameter", 1.0 );
	NEW_OPT( z_min, "input parameter", 0.0 );
	NEW_OPT( z_max, "input parameter", 0.0 );
	NEW_OPT( z_increment, "input parameter", 1.0);
	NEW_OPT( alpha_min, "input parameter", 0.0 );
	NEW_OPT( alpha_max, "input parameter", 0.0 );
	NEW_OPT( alpha_increment, "input parameter", 1.0 );
	NEW_OPT( cosbeta_min, "input parameter", 0.0 );
	NEW_OPT( cosbeta_max, "input parameter", 0.0 );
	NEW_OPT( cosbeta_increment, "input parameter", 1.0 );
	NEW_OPT( gamma_min, "input parameter", 0.0 );
	NEW_OPT( gamma_max, "input parameter", 0.0 );
	NEW_OPT( gamma_increment, "input parameter", 1.0 );
	NEW_OPT( seq, "sequence to model", "cg" );
	NEW_OPT( rmsd_cutoff, "Cutoff for rigid body clustering", 3.0 );
	NEW_OPT( rep_cutoff, "Cutoff in steric check on fa_rep", 2.0 );
	NEW_OPT( input_base1_torsion_M_v_lists, "input list for base1 in two_base_pair", "blah.txt" );
	NEW_OPT( output_base1_torsion_M_v_lists, "input list for base1 in two_base_pair", "blah.txt" );
	NEW_OPT( input_base2_torsion_M_v_lists, "input list for base2 in two_base_pair", "blah.txt" );
	NEW_OPT( output_base2_torsion_M_v_lists, "input list for base2 in two_base_pair", "blah.txt" );
	NEW_OPT( num_torsion_list1, "for test runs on two base pairs, number of strand1 torsions to sample", 400 );
	NEW_OPT( score_diff_cut, "score difference cut for clustering", 10.0 );
	NEW_OPT( minimize_sidechain_res, "in base_pair_to_base_pair, also minimize chi,pucker of these side-chains", blank_size_vector );
	NEW_OPT( minimize_jump, "in base_pair_to_base_pair, also minimize jump", false );
	NEW_OPT( reverse_doublet, "take a base doublet silent file and siwtch residues 1 and 2", false );
	NEW_OPT( output_all, "output all decoys from base pair to base pair", false );
	NEW_OPT( split_silent_files, "output clustered states in separate silent files", false );

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
