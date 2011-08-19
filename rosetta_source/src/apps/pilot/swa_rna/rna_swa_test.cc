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

#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <basic/database/open.hh>
#include <core/init.hh>
#include <core/io/pdb/pose_io.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
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

#include <protocols/viewer/viewers.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh> //for EnergyMap
#include <core/scoring/EnergyMap.fwd.hh> //for EnergyMap
#include <basic/basic.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>


#include <string>
#include <map>

//////////////////////////////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Minimizer.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueSampler.hh>
#include <protocols/swa/rna/StepWiseRNA_PoseSetup.hh>
#include <protocols/swa/rna/StepWiseRNA_Clusterer.hh>
//#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
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

//silly using/typedef
#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


//#include <basic/tracer.hh>
using basic::T;
using basic::Error;
using basic::Warning;
using namespace core;
using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

OPT_KEY( Integer, sample_res )
OPT_KEY( IntegerVector, input_res )
OPT_KEY( IntegerVector, input_res2 )
OPT_KEY( IntegerVector, cutpoint_open )
OPT_KEY( Integer, cutpoint_closed )
OPT_KEY( IntegerVector, fixed_res )
OPT_KEY( IntegerVector, bulge_res )
OPT_KEY( IntegerVector, terminal_res )
OPT_KEY( Boolean, prepend )
OPT_KEY( Boolean, no_o2star_screen )
OPT_KEY( Boolean, o2star_screen )
OPT_KEY( Boolean, fast )
OPT_KEY( Boolean, allow_bulge_at_chainbreak )
OPT_KEY( Boolean, sampler_verbose )
OPT_KEY( String,  params_file )
OPT_KEY( Real, score_diff_cut )
OPT_KEY( String, 	algorithm)
OPT_KEY( String, 	cluster_type)


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
		core::import_pose::pose_from_pdb( *native_pose, *rsd_set, option[ in::file::native ]() );
	}

  Pose pose;

	////////////////////////////////////////////////////////////////////
	// StepWisePoseSetup should create the starting pose.
	// This class might eventually be united with the protein StepWisePoseSetup.
	utility::vector1< std::string > pdb_tags, silent_files_in;

	// First read in any information on pdb read in from silent files.
	if ( option[ in::file::silent ].user() ) {
		silent_files_in = option[ in::file::silent ]();
		pdb_tags = option[ in::file::tags ]();
	}

	if ( option[ in::file::s ].user() ) {
		// Then any pdbs that need to be read in from disk.
		utility::vector1< std::string > const	pdb_tags_from_disk( option[ in::file::s ]() );
		for ( Size n = 1; n <= pdb_tags_from_disk.size(); n++ ) pdb_tags.push_back( pdb_tags_from_disk[ n ] );
	}

	assert( pdb_tags.size() > 0 );

	if ( !option[ sample_res ].user() ) utility_exit_with_message( "Must supply sample_res!" );
  Size const moving_res = option[ sample_res ];
	StepWiseRNA_PoseSetup stepwise_rna_pose_setup( moving_res /*later generalize to more than one residue*/,
																								 desired_sequence,
																								 pdb_tags,
																								 silent_files_in,
																								 option[ input_res ](),
																								 option[ input_res2 ](),
																								 option[ cutpoint_open ](),
																								 option[ cutpoint_closed ]() );
	stepwise_rna_pose_setup.set_native_pose( native_pose );
	stepwise_rna_pose_setup.set_fixed_res( option[ fixed_res ]() );
	stepwise_rna_pose_setup.set_bulge_res( option[ bulge_res ]() );
	stepwise_rna_pose_setup.set_terminal_res( option[ terminal_res ]() );
	stepwise_rna_pose_setup.apply( pose );

  //////////////////////////////////////////////////////////////////
	StepWiseRNA_JobParametersCOP job_parameters( stepwise_rna_pose_setup.job_parameters() );
	StepWiseRNA_BaseCentroidScreenerOP base_centroid_screener = new StepWiseRNA_BaseCentroidScreener( pose, job_parameters );

  //////////////////////////////////////////////////////////////////
	core::scoring::ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( "single_strand_benchmark" );
	if ( option[ score::weights ].user() ) scorefxn = getScoreFunction();

  protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	StepWiseRNA_ResidueSampler stepwise_rna_residue_sampler( job_parameters );
  std::string const silent_file = option[ out::file::silent  ]();
	stepwise_rna_residue_sampler.set_silent_file( "sample_"+silent_file );
	stepwise_rna_residue_sampler.set_scorefxn( scorefxn );
	stepwise_rna_residue_sampler.set_fast( option[ fast ]() );
	stepwise_rna_residue_sampler.set_allow_bulge_at_chainbreak( option[ allow_bulge_at_chainbreak ]() );
	stepwise_rna_residue_sampler.set_o2star_screen( !option[ no_o2star_screen ]() );
	stepwise_rna_residue_sampler.set_verbose( option[ sampler_verbose ]() );
	stepwise_rna_residue_sampler.set_base_centroid_screener( base_centroid_screener );
	stepwise_rna_residue_sampler.apply( pose );

	// let's output the final pose_data_list, just to have a look
	stepwise_rna_residue_sampler.output_pose_data_list( "final_sample_"+silent_file );

  ////////////////////////////////////////////////////////////////
	StepWiseRNA_Minimizer stepwise_rna_minimizer( stepwise_rna_residue_sampler.get_pose_data_list(),
																								job_parameters );
	stepwise_rna_minimizer.set_silent_file( silent_file );
	stepwise_rna_minimizer.set_scorefxn( scorefxn );
	stepwise_rna_minimizer.set_base_centroid_screener( base_centroid_screener );
	stepwise_rna_minimizer.apply( pose );

}

/////////////////////////////////////////////////////////////////////////////
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

	utility::vector1< std::string > const silent_files_in( option[ in::file::silent ]() );

	// Replace this with Parin's clusterer when he checks it in!
	protocols::swa::StepWiseClusterer stepwise_clusterer( silent_files_in );

	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer.set_max_decoys( max_decoys );

	stepwise_clusterer.set_cluster_radius(	option[ OptionKeys::cluster::radius ]()	);
	stepwise_clusterer.set_cluster_by_all_atom_rmsd( true );
	stepwise_clusterer.set_score_diff_cut( option[ score_diff_cut ] );
	stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );

	stepwise_clusterer.cluster();

	std::string const silent_file_out( option[ out::file::silent  ]() );
	stepwise_clusterer.output_silent_file( silent_file_out );

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

  using namespace basic::options;

	std::string algorithm_input = option[ algorithm];

	if (algorithm_input=="cluster_test"){
		cluster_test();
	} else if (algorithm_input=="cluster_old"){
		cluster_outfile_test_OLD();
	}	else {
	  rna_resample_test();
		//	} else {
		//		std::cout << "Error no algorithm selected" << std::endl;
	}

  exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
  using namespace basic::options;

	utility::vector1< Size > blank_size_vector;

  NEW_OPT( sample_res, "residue to sample", 5 );
  NEW_OPT( prepend , "prepend ", true );
	NEW_OPT( params_file, "Input file for pairings", "default.prm" );
	// Note that this could be specified in a PDB INFO file -- but how about for silent files?
	NEW_OPT( cluster_type, "cluster_type", "all_atom" );
	NEW_OPT( input_res, "Residues already present in starting file", blank_size_vector );
	NEW_OPT( input_res2, "Residues already present in starting file2", blank_size_vector );
	NEW_OPT( cutpoint_open, "optional: chainbreak in full sequence", blank_size_vector );
	NEW_OPT( cutpoint_closed, "optional: cutpoint at which to apply chain closure", 0 );
	NEW_OPT( fixed_res, "optional: residues to be held fixed in minimizer", blank_size_vector );
	NEW_OPT( terminal_res, "optional: residues that are not allowed to stack during sampling", blank_size_vector );
	NEW_OPT( bulge_res, "optional: residues to be turned into a bulge variant", blank_size_vector );
	NEW_OPT( no_o2star_screen, "Turn off O2* hydrogen sample during ResidueSampler", false );
	NEW_OPT( o2star_screen, "Turn on O2* hydrogen sample during ResidueSampler -- NOW ON BY DEFAULT!!!!", false );
	NEW_OPT( fast, "quick runthrough for debugging", false );
	NEW_OPT( allow_bulge_at_chainbreak, "Allow sampler to replace chainbreak residue with bulge variant if it looks bulged.", false );
	NEW_OPT( sampler_verbose, "verbose ResidueSampler", false );
	NEW_OPT( score_diff_cut, "score difference cut for clustering", 1000000.0 );
	NEW_OPT( algorithm, "Specify algorithm to execute", "");


  ////////////////////////////////////////////////////////////////////////////
  // setup
  ////////////////////////////////////////////////////////////////////////////
  core::init(argc, argv);


  ////////////////////////////////////////////////////////////////////////////
  // end of setup
  ////////////////////////////////////////////////////////////////////////////

  protocols::viewer::viewer_main( my_main );

}



