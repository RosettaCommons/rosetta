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

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <basic/database/open.hh>
#include <core/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/swa.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
///////////////////////////////////////////////////
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <protocols/idealize/idealize.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/viewer/viewers.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh> //for EnergyMap
#include <core/scoring/EnergyMap.fwd.hh> //for EnergyMap
#include <core/scoring/constraints/ConstraintSet.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>


#include <string>
#include <map>

//////////////////////////////////////////////////////////
#include <protocols/swa/rna/RNA_LoopCloseSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_Clusterer.hh>
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/StepWiseUtil.hh>
#include <protocols/swa/StepWisePoseSetup.hh>
#include <protocols/swa/InputStreamWithResidueInfo.hh>
#include <protocols/swa/StepWiseJobParameters.hh>
#include <protocols/swa/StepWiseClusterer.hh>


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

//Added by Parin
//#include <core/scoring/ScoreType.hh>
#include <list>
#include <stdio.h>
#include <math.h>

using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

OPT_KEY( Boolean, parin_favorite_output )
OPT_KEY( Boolean, close_chainbreak )
OPT_KEY( IntegerVector, sample_res )
OPT_KEY( IntegerVector, cutpoint_open )
OPT_KEY( IntegerVector, cutpoint_closed )
OPT_KEY( IntegerVector, fixed_res )
OPT_KEY( IntegerVector, working_res )
OPT_KEY( IntegerVector, superimpose_res )
OPT_KEY( IntegerVector, calc_rms_res )
OPT_KEY( IntegerVector, virtual_res ) //For testing purposes. Use to make a residue dissapear. Parin S, Jan 29, 2010
OPT_KEY( IntegerVector, bulge_res )
OPT_KEY( IntegerVector, terminal_res ) //Sampling res cannot stack with terminal_res Parin S, 2010 (same as my no_stack_res option)
OPT_KEY( IntegerVector, jump_res )
OPT_KEY( Boolean, center_around_native )
OPT_KEY( Boolean, center_around_Aform )
OPT_KEY( Boolean, prepend )
OPT_KEY( Boolean, no_o2star_screen )
OPT_KEY( Boolean, o2star_screen ) //o2star_screen option doesn't seem to be use in the code???? Parin S, 2010
OPT_KEY( Boolean, fast )
OPT_KEY( Boolean, allow_bulge_at_chainbreak )
OPT_KEY( Boolean, sampler_verbose )
OPT_KEY( Boolean, sampler_native_rmsd_screen )
OPT_KEY( Real, sampler_native_rmsd_screen_cutoff )
OPT_KEY( Boolean, auto_tune )
OPT_KEY( Boolean, just_combine )
OPT_KEY( Boolean, skip_minimize )
OPT_KEY( String,  params_file )
OPT_KEY( String,  input_stream_list )
OPT_KEY( String, cst_file )
OPT_KEY( Real, score_diff_cut )
OPT_KEY( Real, torsion_increment )
OPT_KEY( Real, loop_closer_rep_cutoff )
OPT_KEY( String, 	algorithm)
OPT_KEY( String, 	cluster_type)
OPT_KEY( Integer, num_pose_kept)
OPT_KEY( Integer, combo)
OPT_KEY( Integer, chainbreak_res )
OPT_KEY( Integer, minimize_rounds )

void
figure_out_tags_for_combo( utility::vector1< std::string > const & silent_files_in,
													 utility::vector1< std::string > & pdb_tags ); // see below.


void
initialize_input_stream_from_list(  	utility::vector1< protocols::swa::InputStreamWithResidueInfoOP > & input_streams, std::string const input_stream_list ); // see below

void
convert_sfd_to_pose_data_list( core::io::silent::SilentFileData & silent_file_data,
															 utility::vector1< protocols::swa::rna::pose_data_struct2 > & pose_data_list );

core::io::silent::SilentFileDataOP
convert_pose_data_list_to_sfd( utility::vector1< protocols::swa::rna::pose_data_struct2 > & pose_data_list );

utility::vector1< protocols::swa::rna::pose_data_struct2 >
convert_pose_list_to_pose_data_list( protocols::swa::PoseList const & pose_list );


//////////////////////////////////////////////////////////////////////////////////////
void
rna_resample_test()
{

  using namespace core::pose;
  using namespace core::chemical;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace core::kinematics;
  using namespace core::scoring;
	using namespace protocols::swa::rna;
	using namespace protocols::swa;
	using namespace core::io::silent;

	///////////////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	///////////////////////////////
	// Read in sequence.
	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
	std::string const desired_sequence = fasta_sequence->sequence();

	///////////////////////////////
	// Read in native_pose.
	PoseOP native_pose;
	if (option[ in::file::native ].user() ) {
		native_pose = PoseOP( new Pose );
		import_pose::pose_from_pdb( *native_pose, *rsd_set, option[ in::file::native ]() );
	}

	//	std::cout << native_pose->fold_tree() << " " << native_pose->annotated_sequence( true ) << std::endl;
  std::string const silent_file = option[ out::file::silent  ]();
	utility::vector1< InputStreamWithResidueInfoOP > input_streams;
	initialize_input_streams( input_streams );

	if ( option[ input_stream_list ].user() ) initialize_input_stream_from_list( input_streams, option[ input_stream_list]() ); //later move this into InputStreamWithResidueInfo.cc

  Pose pose;
	//if ( !option[ sample_res ].user() ) utility_exit_with_message( "Must supply sample_res!" );
  utility::vector1< core::Size > const moving_res_list = option[ sample_res ]();
	StepWisePoseSetup stepwise_pose_setup( moving_res_list, /*the first element of moving_res_list is the sampling_res*/
																				 desired_sequence,
																				 input_streams,
																				 option[ cutpoint_open ](),
																				 option[ cutpoint_closed ]() );
	stepwise_pose_setup.set_native_pose( native_pose );
	stepwise_pose_setup.set_fixed_res( option[ fixed_res ]() );
	stepwise_pose_setup.set_superimpose_res( option[ superimpose_res ]() );
	stepwise_pose_setup.set_calc_rms_res( option[ calc_rms_res ]() );
	stepwise_pose_setup.set_bulge_res( option[ bulge_res ]() );
	stepwise_pose_setup.set_terminal_res( option[ terminal_res ]() );
	stepwise_pose_setup.set_virtual_res( option[ virtual_res ]() );
	stepwise_pose_setup.set_jump_res( option[ jump_res ]() );
	stepwise_pose_setup.set_parin_favorite_output( option[ parin_favorite_output ]());
	stepwise_pose_setup.set_rsd_set( rsd_set ); // forces RNA!
	if ( option[ cst_file ].user() ) stepwise_pose_setup.set_cst_file( option[ cst_file ]() );
	stepwise_pose_setup.apply( pose );


  //////////////////////////////////////////////////////////////////
	StepWiseJobParametersOP job_parameters_to_change( stepwise_pose_setup.job_parameters() );
	StepWiseJobParametersCOP job_parameters( job_parameters_to_change );
	StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener = new StepWiseRNA_BaseCentroidScreener( pose, job_parameters );

  //////////////////////////////////////////////////////////////////
	core::scoring::ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "single_strand_benchmark" );
	if ( option[ score::weights ].user() ) scorefxn = getScoreFunction();

  protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	SilentFileDataOP sfd = new SilentFileData;
	Real cluster_radius( 0.25 );
	if ( option[ OptionKeys::cluster::radius ].user() ) cluster_radius = option[ OptionKeys::cluster::radius ]();

	if ( option[ just_combine ]() ){
		// ugly hack -- move to its own function? Time to refactor this code to match stepwise_protein_test()?
		// nothing to do, right? Basically like an IdentityPoseSampleGenerator...
		(*scorefxn)( pose );
		std::string const tag( "START" );
		BinaryRNASilentStruct s( pose, tag );
		pose.dump_pdb( "combine.pdb" );
		//		sfd->write_silent_struct( s, silent_file, false );
		sfd->add_structure( s );

	} else if ( option[ close_chainbreak ]() || option[ chainbreak_res ].user()  ) {

		Size chainbreak_res_ = option[ chainbreak_res ]();
		if ( !option[ chainbreak_res ].user() )			chainbreak_res_ = option[ cutpoint_closed ]()[1] ;

		Size const sample_res = job_parameters->working_moving_res_list()[ 1 ];
		RNA_LoopCloseSampler rna_loop_close_sampler( sample_res /*moving_suite*/,
																								 job_parameters->apply_full_to_sub_mapping( chainbreak_res_ ) );
		rna_loop_close_sampler.set_scorefxn( scorefxn );
		rna_loop_close_sampler.set_bin_size( 20 );
		rna_loop_close_sampler.set_rep_cutoff ( option[ loop_closer_rep_cutoff ]() );
		rna_loop_close_sampler.set_torsion_increment( option[ torsion_increment ]() );
		if ( native_pose ) rna_loop_close_sampler.set_native_pose( native_pose );
		rna_loop_close_sampler.set_center_around_native( option[center_around_native]() );
		rna_loop_close_sampler.set_center_around_Aform( option[center_around_Aform]() );
		rna_loop_close_sampler.apply( pose );
		sfd = rna_loop_close_sampler.silent_file_data();

	} else {

		StepWiseRNA_ResidueSampler stepwise_rna_residue_sampler( job_parameters );
		stepwise_rna_residue_sampler.set_silent_file( "sample_"+silent_file );
		stepwise_rna_residue_sampler.set_scorefxn( scorefxn );
		stepwise_rna_residue_sampler.set_num_pose_kept( option[ num_pose_kept ]() );
		stepwise_rna_residue_sampler.set_fast( option[ fast ]() );
		stepwise_rna_residue_sampler.set_native_rmsd_screen( option[ sampler_native_rmsd_screen ]() );
		stepwise_rna_residue_sampler.set_native_rmsd_screen_cutoff( option[ sampler_native_rmsd_screen_cutoff ]() );
		stepwise_rna_residue_sampler.set_allow_bulge_at_chainbreak( option[ allow_bulge_at_chainbreak ]() );
		stepwise_rna_residue_sampler.set_o2star_screen( !option[ no_o2star_screen ]() );
		stepwise_rna_residue_sampler.set_verbose( option[ sampler_verbose ]() );
		stepwise_rna_residue_sampler.set_base_centroid_screener( base_centroid_screener );
		stepwise_rna_residue_sampler.set_cluster_rmsd(	cluster_radius	);
		stepwise_rna_residue_sampler.set_parin_favorite_output( option[ parin_favorite_output ]());
		stepwise_rna_residue_sampler.apply( pose );

		sfd = convert_pose_data_list_to_sfd( stepwise_rna_residue_sampler.get_pose_data_list() ); //silly conversion -- need to use only one type!

	}


	/////////////////////////////
	// Cluster...
	/////////////////////////////
	// Have an option to read structure back in from disk?  may be useful for checkpointing.
	// For now just take in silent structs prepared by the stepwise_residue_sampler.
	protocols::swa::StepWiseClusterer stepwise_clusterer( sfd );
	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer.set_max_decoys( max_decoys );
	stepwise_clusterer.set_cluster_by_all_atom_rmsd( true ); // false by default
	stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );
	stepwise_clusterer.set_calc_rms_res( job_parameters->working_moving_res_list() );
	if (option[ calc_rms_res ].user() )		stepwise_clusterer.set_calc_rms_res( job_parameters->working_calc_rms_res() ); // used for calculating rmsds to native.
	stepwise_clusterer.set_score_diff_cut( option[ score_diff_cut ] );
	stepwise_clusterer.set_cluster_radius( cluster_radius	);
	stepwise_clusterer.cluster();

  ////////////////////////////////////////////////////////////////
	if ( !option[ skip_minimize ]() ){
		StepWiseRNA_Minimizer stepwise_rna_minimizer( convert_pose_list_to_pose_data_list( stepwise_clusterer.clustered_pose_list() ),
																									job_parameters );
		stepwise_rna_minimizer.set_silent_file( silent_file );
		stepwise_rna_minimizer.set_scorefxn( scorefxn );
		stepwise_rna_minimizer.set_base_centroid_screener( base_centroid_screener );
		stepwise_rna_minimizer.set_minimize_rounds( option[ minimize_rounds ]() );
		// Following needs to be fixed/unified.
		//		stepwise_rna_minimizer.set_calc_rms_res( job_parameters->working_calc_rms_res() ); // used for calculating rmsds to native.
		stepwise_rna_minimizer.apply( pose );

		stepwise_clusterer.set_silent_file_data( stepwise_rna_minimizer.silent_file_data() );
		stepwise_clusterer.cluster();
	}

	stepwise_clusterer.output_silent_file( silent_file );

}

/////////////////////////////////////////////////////////////////////////////
utility::vector1< std::pair< core::Real, std::string > >
read_scores_and_tags( std::string const silent_file ){

	using namespace core::io::silent;

	SilentFileData silent_file_data;
	silent_file_data.read_file( silent_file );

	utility::vector1< std::pair< core::Real, std::string > > score_tags_list;
	for ( SilentFileData::iterator iter = silent_file_data.begin();
				iter != silent_file_data.end(); ++iter ) {
		Real const & silent_score = (*iter)->get_energy( "score" );
		std::string const tag = (*iter)->decoy_tag();
		score_tags_list.push_back( std::make_pair( silent_score, tag ) );
	}
	return score_tags_list;

}

/////////////////////////////////////////////////////////////////////////////
void
figure_out_tags_for_combo( utility::vector1< std::string > const & silent_files_in,
													 utility::vector1< std::string > & pdb_tags ){ // see below.

	using namespace std;

	if ( !option[ combo ].user() ||
			 silent_files_in.size() < 2 ) utility_exit_with_message( "silent file, but no tags or combo option?" );

	Size const which_combo( option[ combo ]() );
	utility::vector1< pair< Real, string > > scores_tags1 = read_scores_and_tags( silent_files_in[ 1 ] );
	utility::vector1< pair< Real, string > > scores_tags2 = read_scores_and_tags( silent_files_in[ 2 ] );

	list< pair< Real, pair< string,string> >  > sum_scores_and_two_tags;
	for ( Size n = 1; n <= scores_tags1.size(); n++ ) {
		for ( Size m = 1; m <= scores_tags2.size(); m++ ) {
			Real const sum_scores = scores_tags1[ n ].first + scores_tags2[ m ].first;
			sum_scores_and_two_tags.push_back(
			 make_pair( sum_scores,	make_pair(  scores_tags1[ n ].second, scores_tags2[ m ].second) ) );
		}
	}

	sum_scores_and_two_tags.sort();
	Size n( 0 );
	for ( list< pair< Real, pair< string,string> >  >::const_iterator iter = sum_scores_and_two_tags.begin();
				iter != sum_scores_and_two_tags.end(); iter++ ) {
		n++;
		if ( n > which_combo ) {
			pdb_tags.clear();
			pdb_tags.push_back( (iter->second).first );
			pdb_tags.push_back( (iter->second).second );
			std::cout << "Found " << which_combo << "-th lowest energy tag combination: " <<
				pdb_tags[ 1] << " " << pdb_tags[ 2 ] << std::endl;
			break;
		}
	}

}

/////////////////////////////////////////////////////////////////////////////
// Do we need this here anymore?
void
cluster_test()
{

	using namespace protocols::swa::rna;

	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
	std::string const full_fasta_sequence = fasta_sequence->sequence();

	std::cout << "full_fasta_sequence= " << full_fasta_sequence << std::endl;

	//Test for building res 3 or 4 of the 3bp helix
	utility::vector1<std::string> input_silent_file_list;
	input_silent_file_list.push_back("blah.out"); //output from rna_resample_test();

	utility::vector1<Size> cluster_res_seq_num_list;
	std::map< core::Size, core::Size > res_map;
	std::map< core::Size, bool > Is_prepend_map;


	cluster_res_seq_num_list.push_back(Size(3));
	cluster_res_seq_num_list.push_back(Size(4));

	Is_prepend_map[3]=false; //append
	Is_prepend_map[4]=true; //prepend

	res_map[1]=1;	//full_pose to rebuild_pose res_map
	res_map[2]=2;
	res_map[3]=3;
	res_map[4]=4;
	res_map[5]=5;
	res_map[6]=6;


	StepWiseRNA_Clusterer stepwise_rna_clusterer(input_silent_file_list);
	stepwise_rna_clusterer.set_output_silent_file("cluster.out");
	stepwise_rna_clusterer.set_cluster_mode(std::string(option[cluster_type]));
	stepwise_rna_clusterer.create_cluster_residue_list(cluster_res_seq_num_list, full_fasta_sequence);
	stepwise_rna_clusterer.set_res_map(res_map);
	stepwise_rna_clusterer.set_is_prepend_map(Is_prepend_map);

	stepwise_rna_clusterer.cluster();

}

///////////////////////////////////////////////////////////////
void
cluster_outfile_test_OLD(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::swa;

	utility::vector1< std::string > const silent_files_in( option[ in::file::silent ]() );

	// Replace this with Parin's clusterer when he checks it in!
	protocols::swa::StepWiseClusterer stepwise_clusterer( silent_files_in );

	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer.set_max_decoys( max_decoys );

	utility::vector1< Size > const working_calc_rms_res = convert_to_working_res( option[ calc_rms_res ](), option[ working_res ]() );
	stepwise_clusterer.set_calc_rms_res( working_calc_rms_res );

	// hey do we want this??????????? protein version has it? but supposed we are doing de novo building?
	//stepwise_clusterer.set_force_align( true );

	Real cluster_radius( 0.25 );
	if ( option[ OptionKeys::cluster::radius ].user() ) cluster_radius = option[ OptionKeys::cluster::radius ]();
	stepwise_clusterer.set_cluster_radius( cluster_radius	);

	stepwise_clusterer.set_cluster_by_all_atom_rmsd( true );

	stepwise_clusterer.set_score_diff_cut( option[ score_diff_cut ] );

	stepwise_clusterer.set_auto_tune( option[ auto_tune ] );

	stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );

	stepwise_clusterer.cluster();

	std::string const silent_file_out( option[ out::file::silent  ]() );
	stepwise_clusterer.output_silent_file( silent_file_out );

}

///////////////////////////////////////////////////////////////
void
pucker_torsion_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::id;
	using namespace core::kinematics;
	using namespace core::optimization;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	//create a pose with a single residue.
	pose::Pose pose, start_pose, save_pose;
	std::string sequence = "g";
	core::pose::make_pose_from_sequence( pose,sequence,*rsd_set );

	//Score these suckers.
	ScoreFunctionOP scorefxn = getScoreFunction();
	scorefxn->show( std::cout, pose );

	core::scoring::rna::RNA_FittedTorsionInfo rna_fitted_torsion_info;
	start_pose = pose;

	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.0000025);
	bool const use_nblist( true );
	bool deriv_check_( false );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, deriv_check_, deriv_check_ );
	options.nblist_auto_update( true );
	kinematics::MoveMap mm;
	//	mm.set( TorsionID(1, BB, 4), true ); //delta
	//	mm.set( TorsionID( 1, CHI, 2), true ); //nu2
	//	mm.set( TorsionID( 1, CHI, 3), true ); //nu1
	mm.set_bb ( true );
	mm.set_chi( true  );

	// apply "ideal" north sugar pucker -- score it.
	pose = start_pose;
	pose.set_torsion( TorsionID( 1, BB, 4 ),  rna_fitted_torsion_info.gaussian_parameter_set_delta_north()[1].center );
	pose.set_torsion( TorsionID( 1, BB, 5 ),  rna_fitted_torsion_info.gaussian_parameter_set_epsilon_north()[1].center );
	pose.set_torsion( TorsionID( 1, CHI, 1 ),  rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[1].center );
	pose.set_torsion( TorsionID( 1, CHI, 2 ),  rna_fitted_torsion_info.gaussian_parameter_set_nu2_north()[1].center );
	pose.set_torsion( TorsionID( 1, CHI, 3 ),  rna_fitted_torsion_info.gaussian_parameter_set_nu1_north()[1].center );
	std::cout << std::endl << "APPLY NORTH ANGLES " << std::endl;
	scorefxn->show( std::cout, pose );
	pose.dump_pdb( "north.pdb");
	save_pose = pose;

	std::cout << std::endl << "MINIMIZE NORTH " << std::endl;
	minimizer.run( pose, mm, *scorefxn, options );
	scorefxn->show( std::cout, pose );
	pose.dump_pdb( "north_min.pdb" );
	std::cout << save_pose.torsion( TorsionID( 1, BB, 4 ) ) << "  " << pose.torsion( TorsionID( 1, BB, 4 ) ) << std::endl;
	std::cout << save_pose.torsion( TorsionID( 1, CHI, 2 ) ) << "  " << pose.torsion( TorsionID( 1, CHI, 2 ) ) << std::endl;
	std::cout << save_pose.torsion( TorsionID( 1, CHI, 3 ) ) << "  " << pose.torsion( TorsionID( 1, CHI, 3 ) ) << std::endl;

	// apply "ideal" south sugar pucker -- score it.
	pose = start_pose;
	pose.set_torsion( TorsionID( 1, BB, 4 ),  rna_fitted_torsion_info.gaussian_parameter_set_delta_south()[1].center );
	pose.set_torsion( TorsionID( 1, BB, 5 ),  rna_fitted_torsion_info.gaussian_parameter_set_epsilon_south()[1].center );
	pose.set_torsion( TorsionID( 1, CHI, 1 ),  rna_fitted_torsion_info.gaussian_parameter_set_chi_south()[1].center );
	pose.set_torsion( TorsionID( 1, CHI, 2 ),  rna_fitted_torsion_info.gaussian_parameter_set_nu2_south()[1].center );
	pose.set_torsion( TorsionID( 1, CHI, 3 ),  rna_fitted_torsion_info.gaussian_parameter_set_nu1_south()[1].center );
	std::cout << std::endl << "APPLY SOUTH ANGLES " << std::endl;
	scorefxn->show( std::cout, pose );
	pose.dump_pdb( "south.pdb");
	save_pose = pose;

	std::cout << std::endl << "MINIMIZE SOUTH " << std::endl;
	minimizer.run( pose, mm, *scorefxn, options );
	scorefxn->show( std::cout, pose );
	pose.dump_pdb( "south_min.pdb" );
	std::cout << save_pose.torsion( TorsionID( 1, BB, 4 ) ) << "  " << pose.torsion( TorsionID( 1, BB, 4 ) ) << std::endl;
	std::cout << save_pose.torsion( TorsionID( 1, CHI, 2 ) ) << "  " << pose.torsion( TorsionID( 1, CHI, 2 ) ) << std::endl;
	std::cout << save_pose.torsion( TorsionID( 1, CHI, 3 ) ) << "  " << pose.torsion( TorsionID( 1, CHI, 3 ) ) << std::endl;


}


///////////////////////////////////////////////////////////////
// Move this to InputStreamWithResidueInfo.cc!!!
void
initialize_input_stream_from_list(  	utility::vector1< protocols::swa::InputStreamWithResidueInfoOP > & input_streams, std::string const input_stream_list ){

	using namespace protocols::swa;
	using namespace core::import_pose::pose_stream;

	input_streams.clear();

	utility::io::izstream data( input_stream_list.c_str() );

	std::string line;
	while( getline(data, line) ) {

		std::istringstream is( line );

		bool is_silent( false );
		utility::vector1< std::string > tags, silent_files;
		std::string input_file, tag;
		is >> input_file;

		if ( input_file == "-silent" ){
			is_silent = true;
			is >> input_file;
			silent_files.push_back( input_file );
			is >> tag;
			tags.push_back( tag );
		}

		std::cout << "INPUT_RES! " << input_file;


		utility::vector1< Size > input_res, slice_res;
		Size res( 0 );

		while ( is.good() ){
			is >> res;
			std::cout << ' ' << res;
			input_res.push_back( res );
		}
		std::cout << std::endl;

		InputStreamWithResidueInfoOP stream;
		if (is_silent ) {
			stream = new InputStreamWithResidueInfo( new SilentFilePoseInputStream( silent_files, tags ),
																																						input_res,
																																						slice_res  );
		} else {
			stream = new InputStreamWithResidueInfo( new PDBPoseInputStream( input_file ),
																																						input_res,
																																						slice_res  );
		}

		input_streams.push_back( stream );

	}


}


///////////////////////////////////////////////////////////////
void
convert_sfd_to_pose_data_list( core::io::silent::SilentFileData & silent_file_data,
															 utility::vector1< protocols::swa::rna::pose_data_struct2 > & pose_data_list ){

	using namespace core::io::silent;
  using namespace core::pose;
  using namespace protocols::swa::rna;

	Pose pose;
	pose_data_list.clear();

	for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(),
					end = silent_file_data.end(); iter != end; ++iter ) {

		std::string tag = iter->decoy_tag();
		//std::cout << tag << std::endl;

		pose_data_struct2 current_pose_data;
		current_pose_data.pose_OP=new pose::Pose;
		iter->fill_pose( (*current_pose_data.pose_OP) );
		current_pose_data.score = iter->get_energy( "score" );
		current_pose_data.tag = tag;

		pose_data_list.push_back( current_pose_data );

	}

}


///////////////////////////////////////////////////////////////
core::io::silent::SilentFileDataOP
convert_pose_data_list_to_sfd( utility::vector1< protocols::swa::rna::pose_data_struct2 > & pose_data_list ){

	using namespace core::io::silent;
  using namespace core::pose;
  using namespace protocols::swa::rna;

	SilentFileDataOP sfd = new SilentFileData;

	for ( Size n = 1; n <= pose_data_list.size(); n++ ) {
		BinaryRNASilentStruct s( *pose_data_list[n].pose_OP,  pose_data_list[n].tag );
		sfd->add_structure( s );
	}

	return sfd;

}

////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< protocols::swa::rna::pose_data_struct2 >
convert_pose_list_to_pose_data_list( protocols::swa::PoseList const & pose_list ){

	using namespace core::io::silent;
  using namespace core::pose;
  using namespace protocols::swa::rna;
  using namespace protocols::swa;

	utility::vector1< protocols::swa::rna::pose_data_struct2 > pose_data_list;

	for ( PoseList::const_iterator iter = pose_list.begin(),
					end = pose_list.end(); iter != end; ++iter ) {

		std::string tag = iter->first;
		PoseOP pose_op = iter->second;
		Real score = pose_op->energies().total_energy();

		pose_data_struct2 current_pose_data;
		current_pose_data.score = score;
		current_pose_data.tag = tag;
		current_pose_data.pose_OP = pose_op;

		pose_data_list.push_back( current_pose_data );
	}

	return pose_data_list;
}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

  using namespace basic::options;

	std::string algorithm_input = option[ algorithm];

	//	system(std::string("mkdir pose/").c_str()); //Output all the poses generated by the code in here. Parin S Jan 28, 2010

	if (algorithm_input=="cluster_test"){
		cluster_test();
	} else if (algorithm_input=="cluster_old"){
		cluster_outfile_test_OLD();
	} else if (algorithm_input=="pucker"){
		pucker_torsion_test();
	}	else {
	  rna_resample_test();
		//	} else {
		//		std::cout << "Error no algorithm selected" << std::endl;
	}

	protocols::viewer::clear_conformation_viewers();
  exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
  using namespace basic::options;

	utility::vector1< Size > blank_size_vector;

	NEW_OPT( parin_favorite_output , " parin_favorite_output ", false );
	NEW_OPT( sample_res, "residues to build, the first element is the actual sample res while the other are the bulge residues", blank_size_vector ); //I am here.
  NEW_OPT( prepend , "prepend ", true );
	NEW_OPT( params_file, "Input file for pairings", "default.prm" );
	NEW_OPT( input_stream_list, "Input file for input streams", "" );
	// Note that this could be specified in a PDB INFO file -- but how about for silent files?
	NEW_OPT( cluster_type, "cluster_type", "all_atom" );
	NEW_OPT( cutpoint_open, "optional: chainbreak in full sequence", blank_size_vector );
	NEW_OPT( cutpoint_closed, "optional: cutpoint at which to apply chain closure", blank_size_vector );
	NEW_OPT( fixed_res, "optional: residues to be held fixed in minimizer", blank_size_vector );
	NEW_OPT( working_res, "optional: residues in input pose in numbering of full-length pose. Used in clustering.", blank_size_vector );
	NEW_OPT( virtual_res, "optional: residues to be made virtual", blank_size_vector ); //For testing purposes. Use to make a residue dissapear. Parin S, Jan 29, 2010
	NEW_OPT( jump_res, "optional: residues for defining jumps -- please supply in pairs", blank_size_vector );
	NEW_OPT( superimpose_res, "optional: residues fixed in space which are superimposable in all poses", blank_size_vector );
	NEW_OPT( calc_rms_res, "optional: residues to calculate rmsds over", blank_size_vector );
	NEW_OPT( chainbreak_res, "must specify along with close_chainbreaks", 0 );
	NEW_OPT( minimize_rounds, "number of minimize rounds", 1 );
	NEW_OPT( num_pose_kept, "optional: set_num_pose_kept by ResidueSampler", 108 );
	NEW_OPT( terminal_res, "optional: residues that are not allowed to stack during sampling", blank_size_vector );
	NEW_OPT( bulge_res, "optional: residues to be turned into a bulge variant", blank_size_vector );
	NEW_OPT( no_o2star_screen, "Turn off O2* hydrogen sample during ResidueSampler", false );
	NEW_OPT( o2star_screen, "Turn on O2* hydrogen sample during ResidueSampler -- NOW ON BY DEFAULT!!!!", false );
	NEW_OPT( fast, "quick runthrough for debugging", false );
	NEW_OPT( allow_bulge_at_chainbreak, "Allow sampler to replace chainbreak residue with bulge variant if it looks bulged.", false );
	NEW_OPT( sampler_verbose, "verbose ResidueSampler", false );
	NEW_OPT( sampler_native_rmsd_screen, "native_rmsd_screen ResidueSampler", false );
	NEW_OPT( sampler_native_rmsd_screen_cutoff, "native_rmsd_screen ResidueSampler", 2.0 );
	NEW_OPT( auto_tune, "autotune rmsd for clustering between 0.1A up to 2.0A", false );
	NEW_OPT( skip_minimize, "no minimize step in rna_swa residue sampling", false );
	NEW_OPT( score_diff_cut, "score difference cut for clustering", 1000000.0 );
	NEW_OPT( torsion_increment, "torsion increment in center_around_native mode of loop-closing", 10.0 );
	NEW_OPT( loop_closer_rep_cutoff, "max rep. cutoff for loop-closing", 6.0 );
	NEW_OPT( algorithm, "Specify algorithm to execute", "");
	NEW_OPT( combo, "Sort combinations of tags from two silent files by energy, return n-th lowest energy combo", 0 );
	NEW_OPT( just_combine, "just_combine", false );
	NEW_OPT( close_chainbreak, "close chainbreak", false );
	NEW_OPT( center_around_native, "center around native", false );
	NEW_OPT( center_around_Aform, "center around A-form in loop closing", false );
	NEW_OPT( cst_file, "Input file for constraints", "default.constraints" );


  ////////////////////////////////////////////////////////////////////////////
  // setup
  ////////////////////////////////////////////////////////////////////////////
  core::init(argc, argv);


  ////////////////////////////////////////////////////////////////////////////
  // end of setup
  ////////////////////////////////////////////////////////////////////////////

  protocols::viewer::viewer_main( my_main );

}



