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

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/Sequence.fwd.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/Ramachandran.hh>
#include <protocols/farna/RNA_ProtocolUtil.hh>

#include <protocols/viewer/viewers.hh>

//Mmmm.. constraints.
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/RT.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/DOF_ID.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/rigid/RigidBodyMover.hh>

//StepWiseProtein!
#include <protocols/stepwise/sampling/general/StepWiseClusterer.hh>
#include <protocols/stepwise/sampling/protein/StepWisePoseSetup.hh>
#include <protocols/stepwise/sampling/protein/StepWiseJobParameters.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <protocols/stepwise/sampling/protein/sample_generators/StepWiseDoNothingSampleGenerator.hh>
#include <protocols/stepwise/sampling/protein/sample_generators/StepWiseCombineSampleGenerator.hh>
#include <protocols/stepwise/sampling/protein/sample_generators/StepWiseIdentitySampleGenerator.hh>
#include <protocols/stepwise/sampling/protein/sample_generators/StepWisePoseCombineSampleGenerator.hh>
#include <protocols/stepwise/sampling/protein/PoseFilter.hh>
#include <protocols/stepwise/sampling/protein/PoseFilter_RMSD_Screen.hh>
#include <protocols/stepwise/sampling/protein/StepWiseBetaAntiParallelJumpSampleGenerator.hh>
#include <protocols/stepwise/sampling/protein/StepWiseBetaAntiParallelUtil.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinFilterer.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinLoopBridger.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinPoseMinimizer.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinScreener.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinCCD_Closer.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinUtil.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinPacker.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinFragmentSampleGenerator.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinJumpSampleGenerator.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinMainChainSampleGenerator.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinPacker.hh>
#include <protocols/stepwise/sampling/protein/MainChainTorsionSet.hh>
#include <protocols/stepwise/sampling/protein/InputStreamWithResidueInfo.hh>

//clustering
#include <protocols/cluster/cluster.hh>

//GreenPacker
#include <protocols/simple_moves/GreenPacker.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>

#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinaryProteinSilentStruct.hh>
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
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray3D.hh>

//Job dsitributor
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>
#include <deque>
#include <vector>

//silly using/typedef

#include <basic/Tracer.hh>
using basic::T;

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>


using basic::Error;
using basic::Warning;

using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;

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
OPT_KEY( Boolean, cluster_by_all_atom_rmsd )
OPT_KEY( Boolean, calc_rms )
OPT_KEY( Boolean, filter_native_big_bins )
OPT_KEY( Boolean, rename_tags )
OPT_KEY( Boolean, n_terminus )
OPT_KEY( Boolean, c_terminus )
OPT_KEY( Boolean, add_peptide_plane )
OPT_KEY( Boolean, minimize_test )
OPT_KEY( Boolean, deriv_check )
OPT_KEY( Boolean, centroid_screen )
OPT_KEY( Boolean, no_sample_junction )
OPT_KEY( Boolean, centroid_output )
OPT_KEY( Boolean, ghost_loops )
OPT_KEY( Boolean, dump )
OPT_KEY( Boolean, sample_beta )
OPT_KEY( Boolean, generate_beta_database )
OPT_KEY( Boolean, combine_loops )
OPT_KEY( Boolean, skip_minimize )
OPT_KEY( Boolean, rescore_only )
OPT_KEY( Boolean, add_virt_res )
OPT_KEY( Boolean, cart_min )
OPT_KEY( Real, score_diff_cut )
OPT_KEY( Real, centroid_score_diff_cut )
OPT_KEY( Real, rmsd_screen )
OPT_KEY( Integer, sample_residue )
OPT_KEY( Integer, start_res )
OPT_KEY( Integer, end_res )
OPT_KEY( Integer, n_sample )
OPT_KEY( Integer, nstruct_centroid )
OPT_KEY( String, pack_weights )
OPT_KEY( String, cst_file )
OPT_KEY( String, disulfide_file )
OPT_KEY( String, align_pdb )
OPT_KEY( String, centroid_weights )
OPT_KEY( String, min_type )
OPT_KEY( String, start_pdb )
OPT_KEY( String, secstruct )
OPT_KEY( Real, min_tolerance )
OPT_KEY( IntegerVector, sample_res )
OPT_KEY( IntegerVector, virtual_res )
OPT_KEY( IntegerVector, cutpoint_open )
OPT_KEY( IntegerVector, cutpoint_closed )
OPT_KEY( IntegerVector, fixed_res )
OPT_KEY( IntegerVector, working_res )
OPT_KEY( IntegerVector, superimpose_res )
OPT_KEY( IntegerVector, calc_rms_res )
OPT_KEY( IntegerVector, jump_res )
OPT_KEY( IntegerVector, loop_bounds )
OPT_KEY( IntegerVector, working_loop_bounds )
OPT_KEY( Boolean, use_green_packer )
OPT_KEY( Boolean, use_packer_instead_of_rotamer_trials )
OPT_KEY( Boolean, auto_tune )
OPT_KEY( Boolean, centroid )
OPT_KEY( Boolean, global_optimize )
OPT_KEY( Boolean, disallow_backbone_sampling )
OPT_KEY( Boolean, disable_sampling_of_loop_takeoff )
OPT_KEY( Boolean, ccd_close )
OPT_KEY( Integer, ccd_close_res )
OPT_KEY( IntegerVector, bridge_res )

void
initialize_native_pose( core::pose::PoseOP & native_pose, core::chemical::ResidueTypeSetCAP & rsd_set );

void
generate_samples_and_cluster( core::pose::Pose & pose,
															protocols::stepwise::sampling::protein::StepWiseJobParametersOP & job_parameters,
															protocols::stepwise::sampling::protein::StepWisePoseSetupOP & stepwise_pose_setup,
															utility::vector1< protocols::stepwise::sampling::protein::InputStreamWithResidueInfoOP > & input_streams,
															utility::vector1 < Size > & moving_residues,
															protocols::stepwise::sampling::general::StepWiseClustererOP & stepwise_clusterer,
															std::string const & silent_file );

void
enable_sampling_of_loop_takeoff( protocols::stepwise::sampling::protein::sample_generators::StepWisePoseSampleGeneratorOP & sample_generator,
																 protocols::stepwise::sampling::protein::StepWiseJobParametersOP job_parameters,
																 pose::Pose & pose );


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
	using namespace protocols::stepwise;
	using namespace protocols::stepwise::sampling::general;
	using namespace protocols::stepwise::sampling::protein;

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
	protocols::stepwise::sampling::protein::initialize_input_streams( input_streams );

  utility::vector1< core::Size > const moving_res_list = option[ sample_res ]();

	// move following to its own function?
	StepWisePoseSetupOP stepwise_pose_setup =
		new StepWisePoseSetup( moving_res_list, /*the first element of moving_res_list is the sampling_res*/
													 desired_sequence,
													 input_streams,
													 option[ cutpoint_open ](),
													 option[ cutpoint_closed ]() );
	stepwise_pose_setup->set_native_pose( native_pose );
	// it would be better to have reasonable defaults for the following...
	stepwise_pose_setup->set_fixed_res( option[ fixed_res ]() );
	if ( option[ superimpose_res ].user() )	stepwise_pose_setup->set_superimpose_res( option[ superimpose_res ]() );
	else stepwise_pose_setup->set_superimpose_res( option[ fixed_res ]() );
	stepwise_pose_setup->set_calc_rms_res( option[ calc_rms_res ]() );
	stepwise_pose_setup->set_jump_res( option[ jump_res ]() );
	stepwise_pose_setup->set_virtual_res( option[ virtual_res ]() );
	stepwise_pose_setup->set_bridge_res( option[ bridge_res ]() );
	stepwise_pose_setup->set_add_peptide_plane_variants( option[ add_peptide_plane ]() );
	stepwise_pose_setup->set_align_file( option[ align_pdb ] );
	stepwise_pose_setup->set_dump( option[ dump ] );
	stepwise_pose_setup->set_rsd_set( rsd_set );
	stepwise_pose_setup->set_secstruct( option[ secstruct ] );
	stepwise_pose_setup->set_cst_file( option[ cst_file ]() );
	stepwise_pose_setup->set_disulfide_file( option[ disulfide_file ]() );
	stepwise_pose_setup->set_add_virt_res( option[ add_virt_res ]() || option[ edensity::mapfile ].user() );

	stepwise_pose_setup->apply( pose );

	StepWiseJobParametersOP & job_parameters = stepwise_pose_setup->job_parameters();
	utility::vector1 < Size > moving_residues = job_parameters->working_moving_res_list();

	if ( option[ dump ]()	) pose.dump_pdb( "after_setup.pdb" );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	std::string const silent_file = option[ out::file::silent  ]();
	std::string const silent_file_minimize = get_file_name( silent_file, "_minimize" );

	////////////////////////////////////////////////////////////////////
	StepWiseClustererOP stepwise_clusterer;
	generate_samples_and_cluster( pose, job_parameters, stepwise_pose_setup, input_streams,
																moving_residues, stepwise_clusterer, silent_file );


	////////////////////////////////////////////////////////////////////
	// move following to its own function?
	// Minimize...
	//	PoseList minimize_pose_list = stepwise_clusterer.clustered_pose_list();
	StepWiseProteinPoseMinimizer stepwise_pose_minimizer( stepwise_clusterer->silent_file_data(), moving_residues );

	ScoreFunctionOP minimize_scorefxn;
	if ( ! option[ score::weights ].user() ) minimize_scorefxn =ScoreFunctionFactory::create_score_function( "score12_no_hb_env_dep.wts"  );
	else minimize_scorefxn = core::scoring::getScoreFunction();
	if (minimize_scorefxn->get_weight( atom_pair_constraint ) == 0.0) minimize_scorefxn->set_weight( atom_pair_constraint, 1.0 ); //go ahead and turn these on
	if (minimize_scorefxn->get_weight( coordinate_constraint) == 0.0) minimize_scorefxn->set_weight( coordinate_constraint, 1.0 ); // go ahead and turn these on
	check_scorefxn_has_constraint_terms_if_pose_has_constraints( pose, minimize_scorefxn );
	minimize_scorefxn->set_weight( linear_chainbreak, 150.0 );
	if ( option[cart_min]() && ( minimize_scorefxn->get_weight( cart_bonded ) == 0.0 ) ) minimize_scorefxn->set_weight( cart_bonded, 1.0 );
	if ( option[edensity::mapfile].user() && minimize_scorefxn->get_weight( elec_dens_atomwise ) == 0.0 ) minimize_scorefxn->set_weight( elec_dens_atomwise, 10.0 );
	stepwise_pose_minimizer.set_scorefxn( minimize_scorefxn );

	if( option[ dump ] ) stepwise_pose_minimizer.set_silent_file( silent_file_minimize );
	stepwise_pose_minimizer.set_move_jumps_between_chains( option[ move_jumps_between_chains ]() );
	//stepwise_pose_minimizer.set_constraint_set( cst_set );
	stepwise_pose_minimizer.set_native_pose( job_parameters->working_native_pose() );
	stepwise_pose_minimizer.set_calc_rms_res( job_parameters->working_calc_rms_res() ); // used for calculating rmsds to native.
	stepwise_pose_minimizer.set_fixed_res( job_parameters->working_fixed_res() );
	stepwise_pose_minimizer.set_move_takeoff_torsions( !option[ disable_sampling_of_loop_takeoff ]() );
	stepwise_pose_minimizer.set_rescore_only( option[ rescore_only ]() );
	stepwise_pose_minimizer.set_cartesian( option[ cart_min ]() );
	if ( option[ min_type ].user() )	stepwise_pose_minimizer.set_min_type( option[ min_type ]() );
	if ( option[ min_tolerance ].user() ) stepwise_pose_minimizer.set_min_tolerance( option[ min_tolerance ]() );

	if ( !option[ skip_minimize ]() ){

		stepwise_pose_minimizer.apply( pose );
		stepwise_clusterer->set_silent_file_data( stepwise_pose_minimizer.silent_file_data() );
		stepwise_clusterer->cluster();

	}

	stepwise_clusterer->output_silent_file( silent_file );

}

////////////////////////////////////////////////////////////////
void
initialize_native_pose( core::pose::PoseOP & native_pose, core::chemical::ResidueTypeSetCAP & rsd_set ){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;

	if ( !option[ in::file::native ].user() ) return;

	native_pose = new Pose;

	std::string native_pdb_file  = option[ in::file::native ];
	import_pose::pose_from_pdb( *native_pose, *rsd_set, native_pdb_file );

	native_pose->conformation().detect_disulfides();
	if (!option[ disulfide_file ].user() ){
		for (Size n = 1; n <= native_pose->total_residue(); n++){
			if ( native_pose->residue_type( n ).has_variant_type( chemical::DISULFIDE ) ) utility_exit_with_message( "native pose has disulfides -- you should probable specify disulfides with -disulfide_file" );
		}
	}

	// this is weird, but I'm trying to reduce memory footprint by not saving the big 2-body energy arrays that get allocated
	// when native_poses with missing atoms are 'packed' during pose_from_pdb().
	//	native_pose->set_new_energies_object( 0 );
	if ( option[ dump ]() ) native_pose->dump_pdb("full_native.pdb");


}

////////////////////////////////////////////////////////////////////////////////////////////
void
generate_samples_and_cluster( core::pose::Pose & pose,
															protocols::stepwise::sampling::protein::StepWiseJobParametersOP & job_parameters,
															protocols::stepwise::sampling::protein::StepWisePoseSetupOP & stepwise_pose_setup,
															utility::vector1< protocols::stepwise::sampling::protein::InputStreamWithResidueInfoOP > & input_streams,
															utility::vector1 < Size > & moving_residues,
															protocols::stepwise::sampling::general::StepWiseClustererOP & stepwise_clusterer,
															std::string const & silent_file ){

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
	using namespace protocols::stepwise;
	using namespace protocols::stepwise::sampling::protein;
	using namespace protocols::stepwise::sampling::protein::sample_generators;

	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	// What are we going to do? SampleGenerator encodes this information...
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	// Usually move_all is specified by sample_res...
	//bool optimize_all( false );  // unused ~Labonte
	sample_generators::StepWisePoseSampleGeneratorOP sample_generator;
	bool full_optimize = option[ global_optimize ]();
	bool disallow_backbone_sampling_ = option[ disallow_backbone_sampling ]();

	std::string const silent_file_sample   = get_file_name( silent_file, "_pack" );
	std::string const silent_file_centroid = get_file_name( silent_file, "_centroid" );

	if ( option[ in::file::frag_files ].user() ){

		std::string const frag_file  = option[ in::file::frag_files ]()[ 1 ];
		utility::vector1< Size > const & slice_res= job_parameters->working_res_list();
		sample_generator = new StepWiseProteinFragmentSampleGenerator( frag_file, slice_res, moving_residues );
		std::cout << "Using StepWiseProteinFragmentSampleGenerator" << std::endl;

		if ( input_streams.size() == 1 ){
			// Could follow this up with an Aligner. [That would remove the "stepwise_pose_setup" silliness].
			// Maybe the identity generator could be called an aligner?
			sample_generators::StepWisePoseSampleGeneratorOP sample_generator_identity = new StepWiseIdentitySampleGenerator( input_streams[1],
																																																		 stepwise_pose_setup /*to align poses. I'm not a big fan of this!*/ );
			sample_generator = new StepWiseCombineSampleGenerator( sample_generator_identity, sample_generator /*fragment generator above*/);
		}

	} else if ( input_streams.size() == 2 ){
		// assume that we want to "combine" two streams of poses...
		// This would be the mode if we have a bunch of templates from which we will graft chunks.
		// Or if we have SWA-based little fragments that we want to paste in.
		sample_generator = new StepWisePoseCombineSampleGenerator( input_streams,
																															 stepwise_pose_setup /*to align poses. I'm not a big fan of this!*/ );
		std::cout << "Using StepWisePoseCombineSampleGenerator" << std::endl;

	} else	if ( option[ sample_beta ]() ) {

		if ( moving_residues.size() !=  1 ) utility_exit_with_message( "Sample beta only works for adding one residue to a beta sheet...");
		sample_generator = new StepWiseBetaAntiParallelJumpSampleGenerator( pose, moving_residues[1] );
		std::cout << "Using StepWiseBetaAntiParallelJumpSampleGenerator" << std::endl;

	} else if ( moving_residues.size() > 0 && (!disallow_backbone_sampling_) ){

		//////////////////////////////////////////////////////////////////////
		//  DEFAULT -- enumeration of conformations for moving residues.
		// Screen that predefines (phi, psi, omega)  for moving residues
		// --> input for later mover carries out all the green packer moves.
		//////////////////////////////////////////////////////////////////////
		StepWiseProteinScreener stepwise_screener( job_parameters );
		stepwise_screener.set_n_sample( option[ n_sample ]() );
		stepwise_screener.set_native_pose( job_parameters->working_native_pose() );

		///////////////////////////////////////////////////////////////
		// Following has not been tested in a while, may not work:
		stepwise_screener.set_filter_native_big_bins( option[ filter_native_big_bins ]  );
		if ( option[ centroid_output ] ) stepwise_screener.set_silent_file( silent_file_centroid );
		if ( option[ centroid_screen] ) stepwise_screener.setup_centroid_screen( option[ centroid_score_diff_cut ](), option[ centroid_weights ](), option[ nstruct_centroid ](), option[ ghost_loops ]() );
		// Could also put loose chainbreak closure check here.
		stepwise_screener.apply( pose );
		/////////////////////////////////////////////////////////////////////////////////////
		// This is not ideal -- easy hack to get things working,
		//  but there's no reason to run the screener first.
		//  should probably get rid of the StepWiseProteinMainChainSampleGenerator...
		/////////////////////////////////////////////////////////////////////////////////////
		sample_generator = new StepWiseProteinMainChainSampleGenerator( stepwise_screener.which_torsions(),
																																		stepwise_screener.main_chain_torsion_set_lists_real() );
		std::cout << "Using StepWiseProteinMainChainSampleGenerator. Num poses: " << stepwise_screener.main_chain_torsion_set_lists_real().size() << std::endl;

	} else { // Just give the poses one at a time... useful in prepacks.

		sample_generator = new StepWiseIdentitySampleGenerator( input_streams[1],
																														stepwise_pose_setup /*to align poses. I'm not a big fan of this!*/ );
		std::cout << "Using StepWiseIdentitySampleGenerator." << std::endl;
		full_optimize = true;

	}

	if ( full_optimize ) {
		// This should force clustering, minimizing based on RMSD over all residues.
		moving_residues.clear();
		for ( Size i = 1; i <= pose.total_residue(); i++ ) moving_residues.push_back( i );
	}


	///////////////////////////////////////////////////////////////////////////
	// Loop closure.
	///////////////////////////////////////////////////////////////////////////
	bool const close_loops =  option[ ccd_close ]()  || option[ bridge_res ].user();
	if (  close_loops ){
		// move following to its own function!!!
		if ( !option[ cutpoint_closed ].user() ) utility_exit_with_message( "must specify -cutpoint_closed with loop closure!" );

		// close loops here, figuring out what torsions are compatible with input poses.
		// Will go through all combinations of poses in the sample generator, and close loops.
		utility::vector1< id::TorsionID > which_torsions;  // will include entire loop.
		utility::vector1<  utility::vector1< Real > > main_chain_torsion_set_lists; // torsions that correspond to closed loops.

		if ( option[ dump ] ) pose.dump_pdb("before_loop_close.pdb");

		if ( option[ ccd_close ]() ) {

			// CCD closure -- heuristic closer but will accept fewer than 6 torsions (and does not require
			//  the torsions to be triaxial like kinematic closer).
			// Involves 5 residues: takeoff - bridge_res1 - bridge_res2 -bridge_res3 - landing
			if ( job_parameters->working_bridge_res().size() >= 3 ) utility_exit_with_message( "cannot specify more than 2 bridge_res for CCD loop closure" );
			StepWiseProteinCCD_Closer stepwise_ccd_closer( sample_generator, job_parameters );
			stepwise_ccd_closer.set_ccd_close_res( option[ ccd_close_res ]() );
			stepwise_ccd_closer.apply( pose );

			which_torsions = stepwise_ccd_closer.which_torsions();
			main_chain_torsion_set_lists = stepwise_ccd_closer.main_chain_torsion_set_lists();

		} else if ( option[ bridge_res ].user() ) {

			// Kinematic Inversion (a.k.a., 'analytical' or 'triaxial') loop closure
			// Involves 5 residues: takeoff - bridge_res1 - bridge_res2 -bridge_res3 - landing
			if ( job_parameters->working_bridge_res().size() != 3 ) utility_exit_with_message( "must specify exactly 3 bridge_res for kinematic loop closure" );

			// sample N-terminal-phi of takeoff and C-terminal-psi of landing.
			if  ( !option[ disable_sampling_of_loop_takeoff ] ) enable_sampling_of_loop_takeoff( sample_generator, job_parameters, pose );
			// stepwiseproteinloopbridger figures out loop residues as those positions that are not 'fixed'.
			StepWiseProteinLoopBridger stepwise_loop_bridger( sample_generator, job_parameters );
			stepwise_loop_bridger.apply( pose );

			which_torsions = stepwise_loop_bridger.which_torsions();
			main_chain_torsion_set_lists = stepwise_loop_bridger.main_chain_torsion_set_lists();

		}

		sample_generator = new StepWiseProteinMainChainSampleGenerator( which_torsions, main_chain_torsion_set_lists );
		moving_residues = merge_vectors( moving_residues, job_parameters->working_bridge_res() );
		std::cout << "Using StepWiseProteinMainChainSampleGenerator with loops" << std::endl;

	}

	//////////////////////////////////////////////////////////////////////////
	PoseFilterOP pose_filter; // used in packer below to, e.g., only keep poses where rebuilt residues are within specific rmsd from target.
	if ( option[ rmsd_screen ].user() ){
		if ( ! job_parameters->working_native_pose() ) utility_exit_with_message( "must specify native pose if using -rmsd_screen!" );
		pose_filter = new PoseFilter_RMSD_Screen( moving_residues, job_parameters->working_native_pose(), option[ rmsd_screen ]() );
	}



	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	// StepWiseProteinPacker -- iterative and enumerative sampling
	//   of backbone degrees of freedom.
	//////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////
	ScoreFunctionOP pack_scorefxn = ScoreFunctionFactory::create_score_function( option[pack_weights] );
	check_scorefxn_has_constraint_terms_if_pose_has_constraints( pose, pack_scorefxn );
	pack_scorefxn->set_weight( linear_chainbreak, 0.2 /*arbitrary*/ );
	if ( option[edensity::mapfile].user() && pack_scorefxn->get_weight( elec_dens_atomwise ) == 0.0 ) pack_scorefxn->set_weight( elec_dens_atomwise, 10.0 );
	StepWiseProteinPacker stepwise_packer( moving_residues, sample_generator );
	stepwise_packer.set_native_pose( job_parameters->working_native_pose() );
	stepwise_packer.set_scorefxn( pack_scorefxn );
	if( option[ dump ] ) 	stepwise_packer.set_silent_file( silent_file_sample /*useful for checkpointing*/ );
	if( pose_filter ) stepwise_packer.set_pose_filter( pose_filter );
	stepwise_packer.set_use_green_packer( option[ use_green_packer ]() );
	stepwise_packer.set_use_packer_instead_of_rotamer_trials( option[ use_packer_instead_of_rotamer_trials ]() );
	stepwise_packer.set_calc_rms_res( job_parameters->working_calc_rms_res() ); // used for calculating rmsds to native.
	stepwise_packer.set_rescore_only( option[ rescore_only ]() );

	stepwise_packer.apply( pose );


	/////////////////////////////
	// Cluster...
	/////////////////////////////
	// Have an option to read structure back in from disk?  may be useful for checkpointing.
	// For now just take in silent structs prepared by the stepwise_residue_sampler.
	stepwise_clusterer = new protocols::stepwise::sampling::general::StepWiseClusterer (  stepwise_packer.silent_file_data() );
	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer->set_max_decoys( max_decoys );
	stepwise_clusterer->set_cluster_by_all_atom_rmsd( option[ cluster_by_all_atom_rmsd ] ); // false by default
	stepwise_clusterer->set_rename_tags( true /*option[ rename_tags ]*/ );
	// this is important (trick from parin ) -- since most of the pose is fixed during
	//  sampling, only calculate rmsd over moving residues for this clustering!
	stepwise_clusterer->set_calc_rms_res( moving_residues );
	if (full_optimize) stepwise_clusterer->set_force_align( true );
	Real cluster_radius( 0.1 );
	if ( option[rescore_only]() ) cluster_radius = 0.0; // no clustering will actually happen
	if ( option[ OptionKeys::cluster::radius ].user() ) cluster_radius = option[ OptionKeys::cluster::radius ]();
	stepwise_clusterer->set_cluster_radius( cluster_radius	);

	stepwise_clusterer->cluster();

	// Perhaps we should output decoys into a silent file at this point -- for checkpointing.
	if ( option[ dump ]()	) stepwise_clusterer->output_silent_file( "CLUSTER_"+silent_file );
}

///////////////////////////////////////////////////////////////
void
cluster_outfile_test(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::stepwise;
	using namespace protocols::stepwise::sampling::general;
	using namespace protocols::stepwise::sampling::protein;

	utility::vector1< std::string > const silent_files_in( option[ in::file::silent ]() );
	StepWiseClusterer stepwise_clusterer( silent_files_in );

	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer.set_max_decoys( max_decoys );

	utility::vector1< Size > const working_calc_rms_res = convert_to_working_res( option[ calc_rms_res ](), option[ working_res ]() );
	stepwise_clusterer.set_calc_rms_res( working_calc_rms_res );

	stepwise_clusterer.set_force_align( true );

	Real cluster_radius( 0.25 );
	if ( option[ OptionKeys::cluster::radius ].user() ) cluster_radius = option[ OptionKeys::cluster::radius ]();
	stepwise_clusterer.set_cluster_radius( cluster_radius	);

	stepwise_clusterer.set_cluster_by_all_atom_rmsd( option[ cluster_by_all_atom_rmsd ] );

	stepwise_clusterer.set_score_diff_cut( option[ score_diff_cut ] );

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
	using namespace protocols::stepwise;
	using namespace protocols::stepwise::sampling::protein;
	using namespace core::chemical;

	PoseOP pose_op,native_pose;
	ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	utility::vector1< std::string > const silent_files_in( option[ in::file::silent ]() );
	SilentFilePoseInputStreamOP input = new SilentFilePoseInputStream( silent_files_in );

	native_pose = PoseOP( new Pose );
	std::string native_pdb_file  = option[ in::file::native ];
	import_pose::pose_from_pdb( *native_pose, *rsd_set, native_pdb_file );

	std::string const silent_file_out( option[ out::file::silent  ]() );
	core::io::silent::SilentFileDataOP sfd_dummy;
	utility::vector1< Size > calc_rms_res_ = option[ calc_rms_res ]();

	// Do we need to slice up native pose?
	if ( option[ working_res ].user() ) {
		pdbslice( *native_pose, option[ working_res ]() );
		calc_rms_res_ = convert_to_working_res( calc_rms_res_, option[ working_res ]() );
	}


	while ( input->has_another_pose() ) {

		PoseOP pose_op( new Pose );
		core::io::silent::SilentStructOP silent_struct( input->next_struct() );
		silent_struct->fill_pose( *pose_op );
		output_silent_struct( *pose_op, native_pose, silent_file_out, silent_struct->decoy_tag(),
													sfd_dummy, calc_rms_res_	);

	}


}



///////////////////////////////////////////////////////////////////////////////////////////////////////////
void
enable_sampling_of_loop_takeoff( protocols::stepwise::sampling::protein::sample_generators::StepWisePoseSampleGeneratorOP & sample_generator,
																 protocols::stepwise::sampling::protein::StepWiseJobParametersOP job_parameters,
																 pose::Pose & pose ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace protocols::stepwise;
	using namespace protocols::stepwise::sampling::protein;

	StepWiseProteinScreener stepwise_screener( job_parameters );
	stepwise_screener.set_n_sample( option[ n_sample ] );

	utility::vector1< Size > takeoff_res;
	Size const pre_loop_res  =  job_parameters->working_bridge_res()[1] - 1;
	Size const post_loop_res =  job_parameters->working_bridge_res()[3] + 1;
	takeoff_res.push_back( pre_loop_res  );
	takeoff_res.push_back( post_loop_res );
	stepwise_screener.set_moving_residues( takeoff_res );

	utility::vector1< Size > fixed_res_for_screener = takeoff_res;
	if ( pre_loop_res  > 1                    ) fixed_res_for_screener.push_back( pre_loop_res  - 1 );
	if ( post_loop_res < pose.total_residue() ) fixed_res_for_screener.push_back( post_loop_res + 1 );
	stepwise_screener.set_fixed_residues( fixed_res_for_screener ); //
	stepwise_screener.apply( pose );

	std::cout << "Going to sample this many takeoff psi/phi combinations: " <<  stepwise_screener.main_chain_torsion_set_lists_real().size() << std::endl;
	sample_generators::StepWisePoseSampleGeneratorOP sample_generator_for_takeoff_res = new StepWiseProteinMainChainSampleGenerator( stepwise_screener.which_torsions(),
																																																								stepwise_screener.main_chain_torsion_set_lists_real() );

	sample_generator = new sample_generators::StepWiseCombineSampleGenerator( sample_generator /*input pose sample generator from above*/, sample_generator_for_takeoff_res);

}


///////////////////////////////////////////////////////////////////////////////
void
combine_loops_test(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::scoring;
	using namespace core::import_pose::pose_stream;
	using namespace core::io::pdb;
	using namespace core::id;
	using namespace core::pose;
	using namespace protocols::stepwise;
	using namespace protocols::stepwise::sampling::protein;

	// Read in main pose.
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	std::string start_pdb_file  = option[ start_pdb ]();
	if ( start_pdb_file.size() == 0 ) utility_exit_with_message( "Must provide -start_pdb" );

	Pose pose;
	import_pose::pose_from_pdb( pose, *rsd_set, start_pdb_file );
	Size const nres = pose.total_residue();
	//	std::cout << pose.annotated_sequence();
	//	pose.dump_pdb( "input_pose.pdb" );


	// Read in loop res.
	utility::vector1<Size> loop_boundaries = option[ loop_bounds ]();
	utility::vector1<Size> working_loop_boundaries = option[ working_loop_bounds ]();
	if ( loop_boundaries.size() != working_loop_boundaries.size() ) utility_exit_with_message( "Must specify equal number of -working_loop_bounds and -loop_bounds" );

	utility::vector1< std::pair< Size, Size > > all_loop_boundaries, all_working_loop_boundaries;
	for ( Size n = 1; n <= (loop_boundaries.size()/2); n++ ){

		Size const loop_start = loop_boundaries[ 2*(n-1) + 1 ];
		Size const loop_end = loop_boundaries[ 2*(n-1) + 2 ];
		all_loop_boundaries.push_back(  std::make_pair( loop_start, loop_end ) );

		Size const working_loop_start = working_loop_boundaries[ 2*(n-1) + 1 ];
		Size const working_loop_end = working_loop_boundaries[ 2*(n-1) + 2 ];
		all_working_loop_boundaries.push_back(  std::make_pair( working_loop_start, working_loop_end ) );

		std::cout << working_loop_end << " " << working_loop_start << " " <<  (working_loop_end - working_loop_start )
							<< "    vs.   "
							<< loop_end << " " << loop_start << " " << ( loop_end - loop_start ) << std::endl;
		if ( (working_loop_end - working_loop_start) != ( loop_end - loop_start ) ){
			utility_exit_with_message( "Problem with loop and working_loop correspondence" );
		}
	}

	utility::vector1<Size> obligate_cuts = option[ cutpoint_open ]();

	utility::vector1< std::string > infiles = option[ in::file::s ]();
	if ( infiles.size() != all_loop_boundaries.size() || infiles.size() == 0 ){
		std::cout << infiles.size() << " " << all_loop_boundaries.size() << std::endl;
		utility_exit_with_message( "Must supply twice as many loop boundaries (-loop_bounds) as files (-s)." );
	}
	Size const num_loops = infiles.size();


	for ( Size n = 1; n <= num_loops; n++ ){

		// // copy_dofs. Cross your fingers!
		Pose import_pose;
		import_pose::pose_from_pdb( import_pose, *rsd_set, infiles[ n ] );
		Size const working_nres =  import_pose.total_residue();

		Size const loop_start = all_loop_boundaries[n].first;
		Size const loop_end   = all_loop_boundaries[n].second;
		Size const working_loop_start = all_working_loop_boundaries[n].first;
		Size const working_loop_end   = all_working_loop_boundaries[n].second;

		std::cout << "check this out: " << working_loop_end << " " << working_loop_start << std::endl;

		bool const loop_start_is_after_cutpoint = ( loop_start == 1  || working_loop_start == 1 || obligate_cuts.has_value( loop_start-1) );
		bool const loop_end_is_before_cutpoint  = ( loop_end == nres || working_loop_end == working_nres || obligate_cuts.has_value( loop_end) );
		assert( !( loop_start_is_after_cutpoint && loop_end_is_before_cutpoint ) );

		FoldTree f( nres );
		if ( loop_start > 1 && loop_end < nres  /*Internal loop requires a cut*/){
			Size cutpoint( loop_end );
			if ( loop_start_is_after_cutpoint ) cutpoint = loop_start-1;
			std::cout << "Putting cutpoint at: " << cutpoint << std::endl;
			f.new_jump( loop_start-1, loop_end+1, cutpoint );
			f.set_jump_atoms( 1, " CA ", " CA ", true /*bKeepStubInResidue. I think this is still protein-centric*/ );
			f.reassign_atoms_for_intra_residue_stubs(); // it seems silly that we need to do this separately.
		}
		pose.fold_tree( f );


		if (false ){
			utility::vector1< Size > input_res, slice_res;
			for ( Size i = 0; i <= ( loop_end - loop_start ); i++ ) {
				input_res.push_back( loop_start + i );
				slice_res.push_back( working_loop_start + i );
			}

			if ( !loop_start_is_after_cutpoint ) {
				input_res.push_back( loop_start - 1 );
				slice_res.push_back( working_loop_start - 1 );
			}
			if ( !loop_end_is_before_cutpoint ) {
				input_res.push_back( loop_end + 1 );
				slice_res.push_back( working_loop_end + 1 );
			}

			InputStreamWithResidueInfo input_stream( new PDBPoseInputStream( infiles[ n ] ), input_res, slice_res );
			input_stream.set_rsd_set( rsd_set );
			input_stream.copy_next_pose_segment( pose );
		}


		/////////////////////////////////////////////////////////////
		std::map < id::AtomID , id::AtomID > atom_id_map;

		for ( Size i = 0; i <= ( loop_end - loop_start ); i++ ) {

			for ( Size j = 1; j <= pose.residue_type( loop_start + i ).natoms(); j++ ){

				std::string const name = pose.residue_type( loop_start + i ).atom_name( j );

				if ( !import_pose.residue_type( working_loop_start + i).has( name ) ){
					std::cout  << "Problem with atom_name: " << name << " at residue " << i << std::endl;
					//					utility_exit_with_message( "PROBLEM!" );
					continue;
				}

				atom_id_map[ AtomID( j, loop_start + i )  ] =
					AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( name, working_loop_start + i ), import_pose  ));
			}
		}

		if ( !loop_start_is_after_cutpoint ){
			std::cout << "adding special atoms for " << loop_start - 1 << std::endl;
			atom_id_map[ AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " N  ", loop_start-1 ), pose))  ] =
				AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " N  ", working_loop_start-1 ), import_pose ));
			atom_id_map[ AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " O  ", loop_start-1 ), pose))  ] =
				AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " O  ", working_loop_start-1 ), import_pose ));
			atom_id_map[ AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " CA ", loop_start-1 ), pose))  ] =
				AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " CA ", working_loop_start-1 ), import_pose ));
			atom_id_map[ AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " C  ", loop_start-1 ), pose))  ] =
				AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " C  ", working_loop_start-1 ), import_pose ));
		}
		if ( !loop_end_is_before_cutpoint ){
			std::cout << "adding special atoms for " << loop_end+1 << std::endl;
			if ( pose.residue_type( loop_end+1).has( " H  " ) ){
				atom_id_map[ AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " H  ", loop_end+1 ), pose))  ] =
					AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " H  ", working_loop_end+1 ), import_pose ));
			}
			atom_id_map[ AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " N  ", loop_end+1 ), pose))  ] =
				AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " N  ", working_loop_end+1 ), import_pose ));
			atom_id_map[ AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " CA ", loop_end+1 ), pose))  ] =
				AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " CA ", working_loop_end+1 ), import_pose ));
			atom_id_map[ AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " C  ", loop_end+1 ), pose))  ] =
				AtomID( core::pose::named_atom_id_to_atom_id( NamedAtomID( " C  ", working_loop_end+1 ), import_pose ));
		}

		//		pose.dump_pdb( "before_copy_dofs.pdb" );
		copy_dofs( pose, import_pose, atom_id_map );
		//		pose.dump_pdb( "after_copy_dofs.pdb" );


	}

	if ( option[ pack_weights ].user() ){

		ScoreFunctionOP pack_scorefxn = ScoreFunctionFactory::create_score_function( option[pack_weights] );
		check_scorefxn_has_constraint_terms_if_pose_has_constraints( pose, pack_scorefxn );
		pack_scorefxn->set_weight( linear_chainbreak, 0.2 /*arbitrary*/ );

		utility::vector1< Size > moving_residues;
		for ( Size i = 1; i <= nres; i++ ) moving_residues.push_back( i );

		// This should be called an IdentityPoseSampleGenerator, and the IdentityPoseSampleGenerator should be called a PoseStreamSampleGenerator, or something like that.
		sample_generators::StepWisePoseSampleGeneratorOP sample_generator = new sample_generators::StepWiseDoNothingSampleGenerator();
		StepWiseProteinPacker stepwise_packer( moving_residues, sample_generator );
		stepwise_packer.set_scorefxn( pack_scorefxn );
		stepwise_packer.set_use_packer_instead_of_rotamer_trials( true );

		stepwise_packer.apply( pose );
	}

	// output the pdb.
	pose.dump_pdb( option[ out::file::o ]() );

	// minimize, and output

}

///////////////////////////////////////////////////////////////
void
big_bins_test(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::scoring;
	using namespace core::import_pose::pose_stream;
	using namespace core::io::pdb;
	using namespace core::id;
	using namespace core::pose;
	using namespace protocols::stepwise::sampling::protein;
	using namespace protocols::stepwise;

	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	std::string start_pdb_file  = option[ in::file::s ]()[1];
	Pose pose;
	import_pose::pose_from_pdb( pose, *rsd_set, start_pdb_file );

	StepWiseJobParametersOP dummy_parameters( new StepWiseJobParameters );
	StepWiseProteinScreener screener( dummy_parameters );

	std::cout << std::endl;
	std::cout << "BIG BINS: " << std::endl;
	for ( Size n = 1; n <= pose.total_residue(); n++ ){
		Size const big_bin = screener.get_big_bin( pose.phi(n), pose.psi(n) );
		std::cout << big_bin;
	}
	std::cout << std::endl;
	std::cout << std::endl;

	std::cout << "SECSTRUCT: " << std::endl;
	for ( Size n = 1; n <= pose.total_residue(); n++ ){
		Size const big_bin = screener.get_big_bin( pose.phi(n), pose.psi(n) );
		char secstruct = 'L';
		if ( big_bin == 1) secstruct = 'H';
		if ( big_bin == 2) secstruct = 'E';
		std::cout << secstruct;
	}
	std::cout << std::endl;

}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace basic::options;

	if ( option[ cluster_test ] ){
		cluster_outfile_test();
	} else if ( option[ generate_beta_database ] ){
		protocols::stepwise::sampling::protein::generate_beta_database_test(); // in StepWiseBetaAntiParallelUtil.cc
	} else if ( option[ combine_loops ] ){
		combine_loops_test();
	} else if ( option[ big_bins ] ){
		big_bins_test();
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

	//Uh, options?
	NEW_OPT( rebuild, "rebuild", false );
	NEW_OPT( cluster_test, "cluster", false );
	NEW_OPT( calc_rms, "calculate rms for input silent file", false );
	NEW_OPT( cluster_by_all_atom_rmsd, "cluster by all atom rmsd", false );
	NEW_OPT( centroid_output, "output centroid structure during screening", false );
	NEW_OPT( n_sample, "number of samples per torsion angle", 18 );
	//	NEW_OPT( filter_rmsd, "for fast sampling", -1.0 );
	NEW_OPT( score_diff_cut, "score difference cut for clustering", 10.0 );
	NEW_OPT( filter_native_big_bins, "Figure out various terms for score12", false );
	NEW_OPT( rename_tags, "After clustering, rename tags S_0, S_1, etc.", false );
	NEW_OPT( n_terminus, "build N terminus", false );
	NEW_OPT( c_terminus, "build C terminus", false );
	NEW_OPT( centroid_screen, "Centroid Screen", false );
	NEW_OPT( skip_minimize, "Skip minimize, e.g. in prepack step", false );
	NEW_OPT( pack_weights, "weights for green packing", "pack_no_hb_env_dep.wts" );
	NEW_OPT( centroid_weights, "weights for centroid filter", "score3.wts" );
	NEW_OPT( cst_file, "Input file for constraints", "" );
	NEW_OPT( disulfide_file, "Input file for disulfides", "" );
	NEW_OPT( align_pdb, "PDB to align to. Default will be native, or no alignment", "" );
	NEW_OPT( min_type, "Minimizer type", "dfpmin_armijo_nonmonotone" );
	NEW_OPT( min_tolerance, "Minimizer tolerance", 0.000025 );
	NEW_OPT( start_pdb, "For combine_loops, parent pdb", "" );
	NEW_OPT( rescore_only, "skip packing,clustering,minimizing -- just get scores & rmsds", false );
	NEW_OPT( no_sample_junction, "disable sampling of residue at junction inherited from start pose", false );
	NEW_OPT( add_peptide_plane, "Include N-acetylation and C-methylamidation caps at termini", false );
	NEW_OPT( ghost_loops, "Virtualize loops in centroid screening", false );
	NEW_OPT( dump, "Dump intermediate silent files", false );
	NEW_OPT( sample_beta, "sample beta strand pairing -- later need to specify parallel/antiparallel", false );
	NEW_OPT( generate_beta_database, "generate_beta_database", false );
	NEW_OPT( combine_loops, "take a bunch of pdbs with loops remodeled and merge them into the parent pose", false );
	NEW_OPT( nstruct_centroid, "Number of decoys to output from centroid screening", 0 );
	NEW_OPT( sample_res, "residues to build, the first element is the actual sample res while the other are the bulge residues", blank_size_vector ); //I am here.
	NEW_OPT( cutpoint_open, "optional: chainbreak in full sequence", blank_size_vector );
	NEW_OPT( cutpoint_closed, "optional: cutpoint at which to apply chain closure", blank_size_vector );
	NEW_OPT( fixed_res, "optional: residues to be held fixed in minimizer", blank_size_vector );
	NEW_OPT( working_res, "optional: residues in input pose in numbering of full-length pose. Used in clustering.", blank_size_vector );
	NEW_OPT( superimpose_res, "optional: residues fixed in space which are superimposable in all poses", blank_size_vector );
	NEW_OPT( calc_rms_res, "optional: residues over which rmsds will be calculated", blank_size_vector );
	NEW_OPT( jump_res, "optional: residues for defining jumps -- please supply in pairs", blank_size_vector );
	NEW_OPT( virtual_res, "optional: residues for defining virtual residues", blank_size_vector );
	NEW_OPT( loop_bounds, "residue numbers for beginning and end of loops", blank_size_vector );
	NEW_OPT( working_loop_bounds, "residue numbers for beginning and end of loops", blank_size_vector );
	NEW_OPT( secstruct, "desired secondary structure for full length pose", "" );
	NEW_OPT( centroid, "Use centroid representation", false );
	NEW_OPT( global_optimize, "In clustering, packing, minimizing, use all residues.", false );
	NEW_OPT( disallow_backbone_sampling, "Just pack and minimize!", false );
	NEW_OPT( disable_sampling_of_loop_takeoff, "For protein loop closure, disallow sampling of psi at N-terminus and phi at C-terminus takeoff residues", false );
	NEW_OPT( use_green_packer, "Use green packer instead of rotamer trials in residue sampling", false );
	NEW_OPT( use_packer_instead_of_rotamer_trials, "Use packer instead of rotamer trials in residue sampling", false );
	NEW_OPT( auto_tune, "autotune rmsd for clustering between 0.1A up to 2.0A", false );
	NEW_OPT( big_bins, "Check out big bin assignment for an input pdb", false );
	NEW_OPT( ccd_close, "Close loops with CCD", false );
	NEW_OPT( move_jumps_between_chains, "Move all jumps", false );
	NEW_OPT( ccd_close_res, "Position at which to close loops with CCD [optional if there is only one cutpoint_closed]", 0);
	NEW_OPT( bridge_res, "instead of enumerative sampling of backbone torsions, combine silent files that contains pieces of loops", blank_size_vector );
	NEW_OPT( rmsd_screen, "keep sampled residues within this rmsd from the native pose", 0.0 );
	NEW_OPT( add_virt_res, "Add virtual residue as root for edensity", false );
	NEW_OPT( cart_min, "Use cartesian minimizer", false );

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
	}
}
