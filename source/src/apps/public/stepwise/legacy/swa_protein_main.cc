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
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Conformation.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/RT.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/viewer/viewers.hh>

//StepWiseProtein!
#include <protocols/stepwise/modeler/StepWiseModeler.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/align/StepWiseLegacyClustererSilentBased.hh>
#include <protocols/stepwise/legacy/modeler/protein/StepWiseProteinPoseSetup.hh>
#include <protocols/stepwise/legacy/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/protein/InputStreamWithResidueInfo.hh>
#include <protocols/stepwise/modeler/file_util.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/sampler/NoOpStepWiseSampler.hh>
#include <protocols/stepwise/sampler/protein/util.hh>

//clustering
#include <protocols/cluster/cluster.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/optimizeH.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/util.hh>
#include <basic/options/option_macros.hh>
#include <core/pose/Pose.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/import_pose/pose_stream/ExtendedPoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <fstream>
#include <iostream>


#include <basic/Tracer.hh>
using basic::T;

static thread_local basic::Tracer TR( "swa_protein_main" );

// option key includes

#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>


using basic::Error;
using basic::Warning;

using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace protocols::stepwise::legacy::modeler;
using namespace protocols::stepwise::legacy::modeler::protein;

using utility::vector1;

using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;
//typedef std::map< std::string, core::pose::PoseOP > PoseList;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( Boolean, big_bins )
OPT_KEY( Boolean, repack )
OPT_KEY( Boolean, move_jumps_between_chains )
OPT_KEY( Boolean, rebuild )
OPT_KEY( Boolean, cluster_test )
OPT_KEY( Boolean, calc_rms )
OPT_KEY( Boolean, rename_tags )
OPT_KEY( Boolean, n_terminus )
OPT_KEY( Boolean, c_terminus )
OPT_KEY( Boolean, add_peptide_plane )
OPT_KEY( Boolean, minimize_test )
OPT_KEY( Boolean, deriv_check )
OPT_KEY( Boolean, no_sample_junction )
OPT_KEY( Boolean, generate_beta_database )
OPT_KEY( Boolean, combine_loops )
OPT_KEY( Boolean, add_virt_res )
OPT_KEY( Integer, sample_residue )
OPT_KEY( Integer, start_res )
OPT_KEY( Integer, end_res )
OPT_KEY( String, cst_file )
OPT_KEY( String, centroid_weights )
OPT_KEY( String, start_pdb )
OPT_KEY( String, secstruct )
OPT_KEY( IntegerVector, sample_res )
OPT_KEY( IntegerVector, superimpose_res )
OPT_KEY( IntegerVector, loop_bounds )
OPT_KEY( IntegerVector, working_loop_bounds )
OPT_KEY( Boolean, auto_tune )
OPT_KEY( Boolean, centroid )

using namespace protocols::stepwise::modeler;
using namespace protocols::stepwise::modeler::protein;

void
initialize_native_pose( core::pose::PoseOP & native_pose, core::chemical::ResidueTypeSetCAP & rsd_set );


///////////////////////////////////////////////////////////////////////
void
rebuild_test(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::kinematics;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::io::silent;
	using namespace core::pose;
	using namespace core::pack;
	using namespace protocols::stepwise::modeler;
	using namespace align;

	// A lot of the following might be better handled by a JobDistributor!?

	////////////////////////////////////////////////////
	//Read in sequence information and native
	////////////////////////////////////////////////////
	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	bool const centroid_mode = option[ centroid ]();
	if ( centroid_mode ) rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( CENTROID );

	//Read in desired fasta.
	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
	std::string const desired_sequence = fasta_sequence->sequence();

	// Read in native pose.
	PoseOP native_pose;
	initialize_native_pose( native_pose, rsd_set );

	////////////////////////////////////////////////////
	// Actual read in of any starting poses.
	////////////////////////////////////////////////////
	Pose pose;

	////////////////////////////////////////////////////////////////////
	utility::vector1< InputStreamWithResidueInfoOP > input_streams;
	initialize_input_streams( input_streams );

  utility::vector1< core::Size > const moving_res_list = option[ sample_res ]();

	// move following to its own function?
	StepWiseProteinPoseSetupOP stepwise_pose_setup( new StepWiseProteinPoseSetup( moving_res_list, /*the first element of moving_res_list is the modeler_res*/
																	desired_sequence,
																	input_streams,
																	option[ OptionKeys::full_model::cutpoint_open ](),
																	option[ OptionKeys::full_model::cutpoint_closed ]() ) );
	stepwise_pose_setup->set_native_pose( native_pose );
	// it would be better to have reasonable defaults for the following...
	stepwise_pose_setup->set_fixed_res( option[ OptionKeys::stepwise::fixed_res ]() );
	if ( option[ superimpose_res ].user() )	stepwise_pose_setup->set_superimpose_res( option[ superimpose_res ]() );
	else stepwise_pose_setup->set_superimpose_res( option[ OptionKeys::stepwise::fixed_res ]() );
	stepwise_pose_setup->set_calc_rms_res( option[ OptionKeys::full_model::calc_rms_res ]() );
	stepwise_pose_setup->set_jump_res( option[ OptionKeys::full_model::jump_res ]() );
	stepwise_pose_setup->set_virtual_res( option[ OptionKeys::full_model::virtual_res ]() );
	stepwise_pose_setup->set_bridge_res( option[ OptionKeys::stepwise::protein::bridge_res ]() );
	stepwise_pose_setup->set_add_peptide_plane_variants( option[ add_peptide_plane ]() );
	stepwise_pose_setup->set_align_file( option[ OptionKeys::stepwise::align_pdb ] );
	stepwise_pose_setup->set_dump( option[ OptionKeys::stepwise::dump ] );
	stepwise_pose_setup->set_rsd_set( rsd_set );
	stepwise_pose_setup->set_secstruct( option[ secstruct ] );
	stepwise_pose_setup->set_cst_file( option[ cst_file ]() );
	stepwise_pose_setup->set_disulfide_file( option[ OptionKeys::stepwise::protein::disulfide_file ]() );
	stepwise_pose_setup->set_add_virt_res( option[ add_virt_res ]() || option[ edensity::mapfile ].user() );

	stepwise_pose_setup->apply( pose );

	// NOTE: These are new options that have been explored in StepWise MonteCarlo and are set to true...
	//  when we go back to stepwise enumeration, let's make sure user is forced to recognize this.
	runtime_assert( option[ OptionKeys::stepwise::protein::protein_prepack ].user() );
	runtime_assert( option[ OptionKeys::stepwise::atr_rep_screen ].user() );
	runtime_assert( option[ OptionKeys::stepwise::protein::allow_virtual_side_chains ].user() );

	working_parameters::StepWiseWorkingParametersOP & working_parameters = stepwise_pose_setup->working_parameters();

	Vector center_vector = ( native_pose != 0 ) ? get_center_of_mass( *native_pose ) : Vector( 0.0 );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400, false, ( native_pose != 0 ), center_vector );

	options::StepWiseModelerOptionsOP stepwise_options( new options::StepWiseModelerOptions );
	stepwise_options->initialize_from_command_line();
	stepwise_options->set_output_minimized_pose_list( true );
	stepwise_options->set_silent_file( option[ out::file::silent ]() );
	stepwise_options->set_disallow_realign( true );
	if ( stepwise_options->pack_weights().size() == 0 ) stepwise_options->set_pack_weights( "stepwise/protein/pack_no_hb_env_dep.wts" );
	if ( stepwise_options->dump()	) pose.dump_pdb( "after_setup.pdb" );
	remove_silent_file_if_it_exists( option[ out::file::silent] );

	ScoreFunctionOP scorefxn;
	if ( ! option[ score::weights ].user() ) scorefxn = ScoreFunctionFactory::create_score_function( "score12_no_hb_env_dep.wts"  );
	else scorefxn = core::scoring::get_score_function();

	StepWiseModeler stepwise_modeler( scorefxn );
	stepwise_modeler.set_native_pose( native_pose );
	stepwise_modeler.set_moving_res_list( working_parameters->working_moving_res_list() );
	stepwise_modeler.set_input_streams( input_streams );
	stepwise_modeler.set_options( stepwise_options );
	if ( working_parameters->working_moving_res_list().size() == 0 ) stepwise_modeler.set_working_prepack_res( get_all_residues( pose ) );
	if ( !option[ OptionKeys::stepwise::test_encapsulation ]() ) stepwise_modeler.set_working_parameters( working_parameters );

	stepwise_modeler.apply( pose );

}


////////////////////////////////////////////////////////////////
void
initialize_native_pose( core::pose::PoseOP & native_pose, core::chemical::ResidueTypeSetCAP & rsd_set ){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;

	if ( !option[ in::file::native ].user() ) return;

	native_pose = core::pose::PoseOP( new Pose );

	std::string native_pdb_file  = option[ in::file::native ];
	import_pose::pose_from_pdb( *native_pose, *rsd_set.lock(), native_pdb_file );

	native_pose->conformation().detect_disulfides();
	if (!option[ OptionKeys::stepwise::protein::disulfide_file ].user() ){
		for (Size n = 1; n <= native_pose->total_residue(); n++){
			if ( native_pose->residue_type( n ).has_variant_type( core::chemical::DISULFIDE ) ) utility_exit_with_message( "native pose has disulfides -- you should probable specify disulfides with -disulfide_file" );
		}
	}

	// this is weird, but I'm trying to reduce memory footprint by not saving the big 2-body energy arrays that get allocated
	// when native_poses with missing atoms are 'packed' during pose_from_pdb().
	//	native_pose->set_new_energies_object( 0 );
	if ( option[ basic::options::OptionKeys::stepwise::dump ]() ) native_pose->dump_pdb("full_native.pdb");

}

///////////////////////////////////////////////////////////////
void
cluster_outfile_test(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::stepwise::modeler;
	using namespace align;

	utility::vector1< std::string > const silent_files_in( option[ in::file::silent ]() );
	StepWiseLegacyClustererSilentBased stepwise_clusterer( silent_files_in );

	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer.set_max_decoys( max_decoys );

	utility::vector1< Size > const working_calc_rms_res = convert_to_working_res( option[ OptionKeys::full_model::calc_rms_res ](),
																																								option[ basic::options::OptionKeys::full_model::working_res ]() );
	stepwise_clusterer.set_calc_rms_res( working_calc_rms_res );

	stepwise_clusterer.set_force_align( true );

	Real cluster_radius( 0.25 );
	if ( option[ OptionKeys::cluster::radius ].user() ) cluster_radius = option[ OptionKeys::cluster::radius ]();
	stepwise_clusterer.set_cluster_radius( cluster_radius	);

	stepwise_clusterer.set_cluster_by_all_atom_rmsd( option[ basic::options::OptionKeys::stepwise::protein::cluster_by_all_atom_rmsd ] );

	stepwise_clusterer.set_score_diff_cut( option[ basic::options::OptionKeys::stepwise::protein::score_diff_cut ] );

	stepwise_clusterer.set_auto_tune( option[ auto_tune ] );

	stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );

	// Do it!
	stepwise_clusterer.cluster();

	std::string const silent_file_out( option[ out::file::silent  ]() );
	stepwise_clusterer.output_silent_file( silent_file_out );

}


///////////////////////////////////////////////////////////////
void
calc_rms_test(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose;
	using namespace protocols::stepwise::modeler;
	using namespace core::chemical;

	PoseOP pose_op,native_pose;
	ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	utility::vector1< std::string > const silent_files_in( option[ in::file::silent ]() );
	SilentFilePoseInputStreamOP input( new SilentFilePoseInputStream( silent_files_in ) );

	native_pose = PoseOP( new Pose );
	std::string native_pdb_file  = option[ in::file::native ];
	import_pose::pose_from_pdb( *native_pose, *rsd_set.lock(), native_pdb_file );

	std::string const silent_file_out( option[ out::file::silent  ]() );
	core::io::silent::SilentFileDataOP sfd_dummy;
	utility::vector1< Size > calc_rms_res_ = option[ OptionKeys::full_model::calc_rms_res ]();

	// Do we need to slice up native pose?
	if ( option[  basic::options::OptionKeys::full_model::working_res ].user() ) {
		pdbslice( *native_pose, option[  basic::options::OptionKeys::full_model::working_res ]() );
		calc_rms_res_ = convert_to_working_res( calc_rms_res_, option[  basic::options::OptionKeys::full_model::working_res ]() );
	}

	while ( input->has_another_pose() ) {
		PoseOP pose_op( new Pose );
		core::io::silent::SilentStructOP silent_struct( input->next_struct() );
		silent_struct->fill_pose( *pose_op );
		output_silent_struct( *pose_op, native_pose, silent_file_out, silent_struct->decoy_tag(),
													sfd_dummy, calc_rms_res_	);
	}

}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace basic::options;

	if ( option[ cluster_test ] ){
		cluster_outfile_test(); // probably should unify with 'PoseSelection' soon.
	} else if ( option[ generate_beta_database ] ){
		protocols::stepwise::sampler::protein::generate_beta_database_test(); // in StepWiseBetaAntiParallelUtil.cc
	} else if ( option[ combine_loops ] ){
		utility_exit_with_message( "-combine_loops is deprecated." );
	} else if ( option[ big_bins ] ){
		utility_exit_with_message( "-big_bins is deprecated." );
	} else if ( option[ calc_rms ] ){
		calc_rms_test();
	}	else {
		rebuild_test();
	}

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );

}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

	using namespace basic::options;

	utility::vector1< Size > blank_size_vector;
	utility::vector1< std::string > blank_string_vector;

	option.add_relevant( OptionKeys::in::file::frag_files );
	option.add_relevant( OptionKeys::full_model::cutpoint_open );
	option.add_relevant( OptionKeys::full_model::cutpoint_closed );
	option.add_relevant( OptionKeys::full_model::virtual_res );
	option.add_relevant( OptionKeys::full_model::jump_res );
	option.add_relevant( OptionKeys::full_model::working_res );
	option.add_relevant( OptionKeys::stepwise::skip_minimize );
	option.add_relevant( OptionKeys::stepwise::min_type );
	option.add_relevant( OptionKeys::stepwise::min_tolerance );
	option.add_relevant( OptionKeys::stepwise::dump );
	option.add_relevant( OptionKeys::stepwise::rmsd_screen );
	option.add_relevant( OptionKeys::stepwise::atr_rep_screen );
	option.add_relevant( OptionKeys::stepwise::use_green_packer );
	option.add_relevant( OptionKeys::stepwise::fixed_res );
	option.add_relevant( OptionKeys::stepwise::align_pdb );
	option.add_relevant( OptionKeys::stepwise::protein::centroid_output );
	option.add_relevant( OptionKeys::stepwise::protein::n_sample );
	option.add_relevant( OptionKeys::stepwise::protein::score_diff_cut );
	option.add_relevant( OptionKeys::stepwise::protein::filter_native_big_bins );
	option.add_relevant( OptionKeys::stepwise::protein::centroid_screen );
	option.add_relevant( OptionKeys::stepwise::protein::centroid_weights );
	option.add_relevant( OptionKeys::stepwise::protein::ghost_loops );
	option.add_relevant( OptionKeys::stepwise::protein::sample_beta );
	option.add_relevant( OptionKeys::stepwise::protein::nstruct_centroid );
	option.add_relevant( OptionKeys::stepwise::protein::global_optimize );
	option.add_relevant( OptionKeys::stepwise::protein::disable_sampling_of_loop_takeoff );
	option.add_relevant( OptionKeys::stepwise::protein::use_packer_instead_of_rotamer_trials );
	option.add_relevant( OptionKeys::stepwise::protein::ccd_close );
	option.add_relevant( OptionKeys::stepwise::protein::move_jumps_between_chains );
	option.add_relevant( OptionKeys::stepwise::protein::bridge_res );
	option.add_relevant( OptionKeys::stepwise::protein::cart_min );
	option.add_relevant( OptionKeys::stepwise::protein::cluster_by_all_atom_rmsd );
	option.add_relevant( OptionKeys::stepwise::protein::protein_prepack );
	option.add_relevant( OptionKeys::stepwise::protein::allow_virtual_side_chains );
	option.add_relevant( OptionKeys::stepwise::protein::disulfide_file );
	NEW_OPT( rebuild, "rebuild", false );
	NEW_OPT( cluster_test, "cluster", false );
	NEW_OPT( calc_rms, "calculate rms for input silent file", false );
	NEW_OPT( rename_tags, "After clustering, rename tags S_0, S_1, etc.", false );
	NEW_OPT( n_terminus, "build N terminus", false );
	NEW_OPT( c_terminus, "build C terminus", false );
	NEW_OPT( cst_file, "Input file for constraints", "" );
	NEW_OPT( start_pdb, "For combine_loops, parent pdb", "" );
	NEW_OPT( no_sample_junction, "disable modeler of residue at junction inherited from start pose", false );
	NEW_OPT( add_peptide_plane, "Include N-acetylation and C-methylamidation caps at termini", false );
	NEW_OPT( generate_beta_database, "generate_beta_database", false );
	NEW_OPT( combine_loops, "take a bunch of pdbs with loops remodeled and merge them into the parent pose", false );
	NEW_OPT( sample_res, "residues to build, the first element is the actual sample res while the other are the bulge residues", blank_size_vector ); //I am here.
	NEW_OPT( superimpose_res, "optional: residues fixed in space which are superimposable in all poses", blank_size_vector );
	NEW_OPT( loop_bounds, "residue numbers for beginning and end of loops", blank_size_vector );
	NEW_OPT( working_loop_bounds, "residue numbers for beginning and end of loops", blank_size_vector );
	NEW_OPT( secstruct, "desired secondary structure for full length pose", "" );
	NEW_OPT( centroid, "Use centroid representation", false );
	NEW_OPT( auto_tune, "autotune rmsd for clustering between 0.1A up to 2.0A", false );
	NEW_OPT( big_bins, "Check out big bin assignment for an input pdb", false );
	NEW_OPT( add_virt_res, "Add virtual residue as root for edensity", false );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);

	if ( option[ OptionKeys::stepwise::protein::allow_virtual_side_chains ]() ) {
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_SIDE_CHAIN" );
	}
	option[ OptionKeys::chemical::patch_selectors ].push_back( "PEPTIDE_CAP" ); // N_acetylated.txt and C_methylamidated.txt

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
