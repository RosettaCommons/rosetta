// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file domain_assembly.cc
/// @brief App that I stole from Frank DiMaio that utilizes the DockPdbIntoDensityMover instead of DockFragmentsIntoDensityMover
///        that app has the capability to split most tasks. This can be used for debugging, or parallel execution
/// @author Danny Farrell

#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/testing.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>

#include <core/fragment/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>

#include <devel/init.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>

#include <protocols/electron_density/DockPDBIntoDensityMover.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <string>


OPT_KEY( Boolean, beta_conv )
OPT_KEY( Boolean, centroid_silent_out )
OPT_KEY( Boolean, cheat_native_com )
OPT_KEY( Boolean, cheat_native_mca )
OPT_KEY( Boolean, convolute_single_residue)
OPT_KEY( Boolean, dump_inter )
OPT_KEY( Boolean, dump_inter_silent )
OPT_KEY( Boolean, gaussian_blur)
OPT_KEY( Boolean, min )
OPT_KEY( Boolean, min_bb )
OPT_KEY( Boolean, rot_middle_ca )
OPT_KEY( Boolean, rot_seq_center )
OPT_KEY( Boolean, score_natives )
OPT_KEY( Integer, bw )
OPT_KEY( Integer, max_rot_per_trans )
OPT_KEY( Integer, movestep )
OPT_KEY( Integer, n_filtered )
OPT_KEY( Integer, n_output )
OPT_KEY( Integer, n_to_search )
OPT_KEY( Integer, ncyc )
OPT_KEY( Integer, refine_end )
OPT_KEY( Integer, refine_start )
OPT_KEY( Integer, search_end )
OPT_KEY( Integer, search_start )
OPT_KEY( Integer, searchsep )
OPT_KEY( String, start_model )

OPT_KEY( IntegerVector, core_idx )
OPT_KEY( Real, clust_radius )
OPT_KEY( Real, constrain_refinement )
OPT_KEY( Boolean, cart_ref )
OPT_KEY( Real, delR )
OPT_KEY( Real, frag_dens )
OPT_KEY( Real, laplacian_offset )
OPT_KEY( Real, point_radius )
OPT_KEY( String, mode )
OPT_KEY( String, man_pts )
OPT_KEY( String, points_to_search_fname )
OPT_KEY( String, points_to_search_pdb_fname )
OPT_KEY( String, point_search_results_fname )
OPT_KEY( String, combined_search_results_fname )
OPT_KEY( StringVector, local_result_files )
OPT_KEY( StringVector, multi_native )
OPT_KEY( String, final_chain )

using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TR( "apps.dgdp" );


void
check_inputs_for_errors() {
	// Check for multiple PDBs to dock, we don't support that
	utility::vector1< std::string > pose_filenames = option[ in::file::s ]();
	if ( pose_filenames.size() != 1 && option[ mode ]() != "cluster_silent" ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "Found more or less then one input pose, please only use one! Exiting...");
	}
	// Check for native flag, we don't support that
	if ( option[ in::file::native ].user() ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "This protocol requires that you use the multi_native flag, not in:file:native, it's okay even if you only have 1 native location!");
	}
	// Check for naming conventions, required for silent entry naming
	if ( option[ mode ]() == "cluster_silent" ) return; // TODO we have to make this smarter instead of just blanket adding if statements
	utility::vector1< std::string > fn = utility::string_split( pose_filenames[1], '/' );
	utility::vector1< std::string > split_filename = utility::string_split( fn.back(), '.' );
	std::string base_filename = split_filename[1];
	utility::vector1< std::string > split_base_filename = utility::string_split( base_filename, '_' );
	// if ( split_base_filename.size() != 2 ) throw std::runtime_error( "Sorry but we require you to name your input pdbs (-s flag) in a specific way. Please follow '{name}_{chain}.pdb >>" + base_filename + "<<" );
	if ( option[ mode ]() == "search_points" ) {
		if ( option[ point_search_results_fname ]() == "" ) {
			throw CREATE_EXCEPTION(utility::excn::BadInput, "The 'search_points' mode of this protocol requires that you set the '-point_search_results_fname' option");
		}
	}

	if ( option[ mode ]() == "combine_search" ) {
		if ( option[ combined_search_results_fname ]() == "" ) {
			throw CREATE_EXCEPTION(utility::excn::BadInput, "The 'search_points' mode of this protocol requires that you set the '-combined_search_results_fname' option");
		}
	}
}


void
set_basic_dock_options( protocols::electron_density::DockPDBIntoDensityMoverOP const & dock ) {
	check_inputs_for_errors();
	// search options
	dock->setB( option[ bw ] );
	dock->setBetaConv( option[ beta_conv ]() );
	dock->setConvoluteSingleR( option[ convolute_single_residue ]());
	dock->setDelR( option[ delR ]() );
	dock->setFragDens(option[ frag_dens ]());
	dock->setGaussianBlur( option[ gaussian_blur ]());
	dock->setGridStep( option[ movestep ] );
	dock->setLaplacianOffset( option[ laplacian_offset ]());
	dock->setMaxRotPerTrans( option[ max_rot_per_trans ]());
	dock->setNCyc(option[ ncyc ]());
	dock->setPointRadius( option[ point_radius ]());
	dock->setPointsToSearchFname( option[ points_to_search_fname ]() );
	dock->setPointsToSearchPDBFname( option[ points_to_search_pdb_fname ]() );
	dock->setPointSearchResultsFname( option[ point_search_results_fname ]() );
	dock->setPointSearchResultsFname( option[ combined_search_results_fname ]() );
	dock->setRotateMiddleCA( option[ rot_middle_ca ]() );
	dock->setRotateSeqCenter( option[ rot_seq_center ]() );
	dock->setTopN( option[ n_to_search ], option[ n_filtered ] , option[ n_output ] );
	dock->setStartModel( option[ start_model ]() );

	// Refinement options
	dock->setConstrainRefinement( option[ constrain_refinement ]());
	dock->setDoRefine(option[ min ]());
	dock->setMinBackbone(option[ min_bb ]());

	// General options
	dock->setClusterRadius(option[ clust_radius ]());
	dock->setOverwrite( option[ out::overwrite ]() );
	dock->setDumpInter( option[ dump_inter ]() );
	dock->setDumpInterSilent( option[ dump_inter_silent ]() );

	// Parallel options
	dock->setCoreIdx( option[ core_idx ]()); // TODO: check for more than 2 values
	dock->setRefineStart( option[refine_start]() );
	dock->setRefineEnd( option[refine_end]() );


	// Setup native entries
	utility::vector1< std::string > const pose_filenames = option[ in::file::s ]();
	utility::vector1< core::pose::PoseOP > all_native_poseOP;
	if ( option[ multi_native ]().size() != 0 ) {
		if ( option[ mode ]() == "cluster_silent" ) {
			// in this case we don't want to pre-trim our native since we will be working with a lot
			// of different poses
			for ( core::Size i = 1; i <= option[ multi_native ]().size(); ++i ) {
				TR << "Setting " << option[ multi_native ]()[i] << " to be indexed to native # " << i << std::endl;
				std::string filename = option[ multi_native ]()[i];
				core::pose::PoseOP current_native_poseOP = core::import_pose::pose_from_file( filename );
				protocols::simple_moves::SwitchResidueTypeSetMover to_all_atom( core::chemical::FA_STANDARD );
				to_all_atom.apply( *current_native_poseOP );
				all_native_poseOP.push_back( current_native_poseOP);
			}
		} else {
			TR << "Found native and setting it now!" << std::endl;
			for ( core::Size i = 1; i <= option[ multi_native ]().size(); ++i ) {
				TR << "Setting " << option[ multi_native ]()[i] << " to be indexed to native # " << i << std::endl;
				std::string filename = option[ multi_native ]()[i];
				core::pose::PoseOP current_native_poseOP = core::import_pose::pose_from_file( filename );
				protocols::simple_moves::SwitchResidueTypeSetMover to_all_atom( core::chemical::FA_STANDARD );
				to_all_atom.apply( *current_native_poseOP );
				all_native_poseOP.push_back( current_native_poseOP);
			}
			// all natives are set, but we want to compare native against our query model, not the best possible
			for ( core::Size i = 1; i <= all_native_poseOP.size(); ++i ) {
				core::pose::PoseOP modelOP = core::import_pose::pose_from_file( pose_filenames[1]);
				core::pose::PoseOP cut_nativeOP = all_native_poseOP[i]->clone();
				for ( core::Size j = cut_nativeOP->size(); j > 0 ; --j ) {
					if ( cut_nativeOP->residue(j).aa() == core::chemical::aa_vrt || cut_nativeOP->size() == 1 ) continue;
					int nat_resnum = cut_nativeOP->pdb_info()->number(j);
					char nat_chain = cut_nativeOP->pdb_info()->chain(j);
					core::Size model_posenum = modelOP->pdb_info()->pdb2pose( nat_chain, nat_resnum );
					//int r2_resnum = r2->pdb_info()->number(r2_posenum);
					// cut if no alignment
					if ( model_posenum == 0 ) {
						cut_nativeOP->delete_residue_slow(j);
						continue;
					}
					if ( cut_nativeOP->residue(j).name1() != modelOP->residue(model_posenum).name1() ) {
						TR.Fatal << "Note, residues are removed at this step, remember that when you're debugging" << std::endl;
						TR.Fatal << "nat seq: " << cut_nativeOP->sequence() << std::endl;
						TR.Fatal << "model seq: " << modelOP->sequence() << std::endl;
						TR.Fatal << "nat chain: " << cut_nativeOP->pdb_info()->chain(j) << " nat pdbinfo res: " << cut_nativeOP->pdb_info()->number(j) << std::endl;
						TR.Fatal << "model chain: " << modelOP->pdb_info()->chain(model_posenum) << " model pdbinfo res: " << modelOP->pdb_info()->number(model_posenum) << std::endl;
						TR.Fatal << "nat name: " << cut_nativeOP->residue(j).name3() << " model name: " << modelOP->residue(model_posenum).name3() << std::endl;
						TR.Fatal << "native pose num: " << j << " model pose num: " << model_posenum << std::endl;
						throw CREATE_EXCEPTION(utility::excn::BadInput, "Found that in compare_and_align_poses, poses did not align in residue numbering or sequence");
					}
				}
				if ( cut_nativeOP->size() <= 1 ) {
					core::Size nativeend = 0;
					core::Size modelend = 0;
					if ( all_native_poseOP[i]->residue(all_native_poseOP[i]->size()).aa() == core::chemical::aa_vrt ) nativeend = all_native_poseOP[i]->size() - 1;
					else nativeend = all_native_poseOP[i]->size();
					if ( modelOP->residue(modelOP->size()).aa() == core::chemical::aa_vrt ) modelend = modelOP->size() - 1;
					else modelend = modelOP->size();
					TR.Warning << "After aligning native and query pose, Found that there were sequence alignment between the two" << std::endl;
					TR.Warning << "native seq: " << all_native_poseOP[i]->sequence() << std::endl;
					TR.Warning << "native start: " << all_native_poseOP[i]->pdb_info()->number(1) << " end: " << all_native_poseOP[i]->pdb_info()->number(nativeend) << std::endl;
					TR.Warning << "query seq: " << modelOP->sequence() << std::endl;
					TR.Warning << "query start: " << modelOP->pdb_info()->number(1) << " end: " << modelOP->pdb_info()->number(modelend) << std::endl;
				}
				core::Real rms = dock->compare_and_align_poses(*modelOP, *all_native_poseOP[i]);
				TR << "Aligned rms to native IDX " << i << " is " << rms << std::endl;
				all_native_poseOP[i] = cut_nativeOP;
			}
		}
		dock->setMultiNative( all_native_poseOP );
		dock->setScoreNatives( option[ score_natives ]() );
	}


	// Predefine search on Native COMs ( cheat )
	if ( option[ cheat_native_com ]() == true ) {
		if ( option[ multi_native ]().size() == 0 ) {
			throw CREATE_EXCEPTION(utility::excn::BadInput, "To 'cheat_native_com' you must include a native using the multi_native flag" );
		} else {
			dock->setCheatNativeCOM( true );
		}
	}
	// Predefine search on Native Middle CAs xyz ( cheat ) ( CA closest to com )
	if ( option[ cheat_native_com ]() == true ) {
		if ( option[ multi_native ]().size() == 0 ) {
			throw CREATE_EXCEPTION(utility::excn::BadInput, "To 'cheat_native_com' you must include a native using the multi_native flag" );
		} else {
			dock->setCheatNativeMCA( true );
		}
	}

	// set silent prefix
	if ( option[ out::file::silent ].user() ) {
		std::string silent_fn = option[ out::file::silent ]();
		dock->setOutputSilent( silent_fn );
	}

	std::string filename;
	// don't use in::file::s if we are clustering ;
	if ( option[ mode ]() == "cluster_silent" ) {
		utility::vector1< std::string > complete_filename_split = utility::string_split( option[out::file::silent](), '/' );
		filename = complete_filename_split.back();
		utility::vector1< std::string > split_filename = utility::string_split( filename, '.' );
		dock->setTag( split_filename[1] );
		return;
	} else {
		// Set tag and check filename organization
		utility::vector1< std::string > complete_filename = option[ in::file::s ]();
		utility::vector1< std::string > complete_filename_split = utility::string_split( complete_filename[1], '/' );
		filename = complete_filename_split.back();
		// remove file extension
		utility::vector1< std::string > split_filename = utility::string_split( filename, '.' );
		// check for multiple '.'s.
		if ( split_filename.size() > 2 ) throw CREATE_EXCEPTION(utility::excn::BadInput, "Sorry we require your filename to have only 1 period ( '.' ) in it. \n Found ==> " + filename + " <== to have more than 1 period");
		dock->setTag( split_filename[1] );
	}
}

void
get_pts( protocols::electron_density::DockPDBIntoDensityMoverOP const & dock ) {
	utility::vector1< std::string > filenames = option[ in::file::s ]();
	if ( filenames.size() != 1 ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "Found more then one input pose, please only use one! exiting");
	}
	core::pose::PoseOP poseOP = core::import_pose::pose_from_file( filenames[1]);

	if ( option[ man_pts ].user() ) {
		core::pose::PoseOP points = core::import_pose::pose_from_file(option[ man_pts ]());
		utility::vector1< numeric::xyzVector< core::Real > > pts;
		for ( core::Size i = 1; i <= points->size(); ++i ) {
			pts.push_back( points->residue(i).nbr_atom_xyz() );
		}
		dock->predefine_search( pts );
	}

	dock->get_points_to_search( poseOP );
}

void
do_search( protocols::electron_density::DockPDBIntoDensityMoverOP const & dock ) { // TODO: set nRsteps in here
	utility::vector1< std::string > const filenames = option[ in::file::s ]();
	if ( filenames.size() != 1 ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "Found more then one input pose, please only use one! exiting");
	}
	core::pose::PoseOP poseOP = core::import_pose::pose_from_file( filenames[1]);

	// check for search parameters
	// if ( option[ search_start ]() > option[ search_end ]() ) {
	//  throw std::runtime_error("Found search start greater than search_end, please check your input flags");
	// }
	// if ( option[ search_start ]() == 0 ||  option[ search_end ]() == 0 ) {
	//  throw std::runtime_error("Found search start or end at 0!, indexing starts at 1 so please check your input flags");
	// }

	dock->apply_search( *poseOP, 10 );
}


void
combine_search( protocols::electron_density::DockPDBIntoDensityMoverOP const & dock ) {
	utility::vector1< std::string > const filenames = option[ in::file::s ]();
	utility::vector1< std::string > const local_result_filenames = option[ local_result_files ]();
	if ( local_result_filenames.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "Found zero local result filenames, please add a flag for those! (hint: *.sscores) ");
	}

	core::pose::PoseOP const poseOP = core::import_pose::pose_from_file( filenames[1]);

	dock->set_nRsteps_from_pose( *poseOP );
	dock->combine_search( local_result_filenames, poseOP );
}


void
results_to_pdb( protocols::electron_density::DockPDBIntoDensityMoverOP const & dock ) {
	utility::vector1< std::string > filenames = option[ in::file::s ]();
	if ( filenames.size() != 1 ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "Found more or less then one input pose, please only use one! exiting");
	}
	utility::vector1< std::string > local_result_filenames = option[ local_result_files ]();
	if ( local_result_filenames.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "Found zero local result filenames, please add a flag for those! (hint: *.sscores) ");
	}

	core::pose::PoseOP const poseOP = core::import_pose::pose_from_file( filenames[1]);

	dock->search_results_to_pdb( local_result_filenames, poseOP );
}


void
do_refinement( protocols::electron_density::DockPDBIntoDensityMoverOP const & dock ) {
	utility::vector1< std::string > const filenames = option[ in::file::s ]();
	if ( filenames.size() != 1 ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "Found more or less then one input pose, please only use one! exiting");
	}

	utility::vector1< std::string > local_result_filenames = option[ local_result_files ]();
	if ( local_result_filenames.size() == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "Found zero local result filenames, please add a flag for those! (hint: *.sscores) ");
	}

	core::pose::PoseOP const poseOP = core::import_pose::pose_from_file( filenames[1]);

	dock->apply_refinement( local_result_filenames, poseOP );
}

void
manual_refine( protocols::electron_density::DockPDBIntoDensityMoverOP const & dock ) {
	utility::vector1< std::string > const filenames = option[ in::file::s ]();
	if ( filenames.size() != 1 ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "Found more or less then one input pose, please only use one! exiting");
	}
	core::pose::PoseOP poseOP = core::import_pose::pose_from_file( filenames[1]);
	dock->manual_refine_pdb( poseOP );
}


void
combine_refinement( protocols::electron_density::DockPDBIntoDensityMoverOP const & dock ) {
	if ( !option[ in::file::silent ].user() ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "We need silent files to combine, please include the in:file:silent flag");
	}
	utility::vector1< std::string > silent_filenames = option[ in::file::silent ]();

	dock->combine_refinement( silent_filenames );
}


void
cluster_silent( protocols::electron_density::DockPDBIntoDensityMoverOP const & dock ) {
	if ( !option[ in::file::silent ].user() ) {
		throw CREATE_EXCEPTION(utility::excn::BadInput, "We need silent files to combine, please include the in:file:silent flag");
	}
	utility::vector1< std::string > silent_filenames = option[ in::file::silent ]();

	dock->setFinalChain( option[ final_chain ]() );
	TR.Warning << "Setting final chain: " << option[ final_chain ]() << std::endl;

	dock->cluster_silent( silent_filenames );
}



// main
int main(int argc, char* argv[]) {
	try {
		utility::vector1< std::string > blank_string_vector;
		utility::vector1< core::Size > blank_size_vector;

		// Options

		// Search options
		NEW_OPT( beta_conv, "Modify convolution to take into account shape of the domain (beta)", false );
		NEW_OPT( bw, "SPHARM bandwidth", 32 );
		NEW_OPT( clust_radius, "Cluster radius", 8.0 );
		NEW_OPT( convolute_single_residue, "Convolute only on middle reside", false);
		NEW_OPT( delR, "Sets the delR", 2 );
		NEW_OPT( frag_dens, "Fragment density", 0.9 );
		NEW_OPT( gaussian_blur, "Blur the map using a gaussian kernel instead of a rotational convolution", false);
		NEW_OPT( laplacian_offset, "Activates laplacian scoring, and sets laplacian filter offset distance ", 0 );
		NEW_OPT( max_rot_per_trans, "Maximum number of rotations we will take at any translational point during density search", 10 );
		NEW_OPT( movestep, "Grid spacing over which to search", 1 );
		NEW_OPT( man_pts, "pdb file of points that you have manually selected to search", "" );
		NEW_OPT( n_filtered,  "How many solutions to take to refinement", 100 );
		NEW_OPT( n_output, "How many solutions to output", 10 );
		NEW_OPT( n_to_search, "How many translations to search", 100 );
		NEW_OPT( ncyc, "Min cycles", 1 );
		NEW_OPT( point_radius, "Minimum translation point to point radius", 4 );
		NEW_OPT( points_to_search_fname, "File containing points to search (input/output)", "" );
		NEW_OPT( points_to_search_pdb_fname, "File of where to dump points to search as a pdb file (output)", "" );
		NEW_OPT( point_search_results_fname, "File of where to dump results of point search", "" );
		NEW_OPT( combined_search_results_fname, "File of where to dump results of point search", "" );
		NEW_OPT( rot_middle_ca, "In SHARM docking, will rotate pose on middle CAs (ca closest to COM)", false );
		NEW_OPT( rot_seq_center, "In SHARM docking, will rotate pose on CA of residue at the center of the sequence", false );
		NEW_OPT( searchsep, "Min distance between search points", 3 );
		NEW_OPT( start_model, "Aligned pdb file that we can use to remove density from the map before docking", "" );

		// Cheat options
		NEW_OPT( cheat_native_com, "Will only dock on native COMs, (point selection cheat)", false );
		NEW_OPT( cheat_native_mca, "Will only dock on native middle CAs (ca closest to COM), (point selection cheat)", false );

		// Refinement options
		NEW_OPT( min_bb, "Minimize backbone?", false );
		NEW_OPT( min, "RigidBody min?", true );
		NEW_OPT( constrain_refinement, "Constrain the pose xyz coordinates during refinement ( 0 is don't constrain! )", 0);
		NEW_OPT( cart_ref, "Use cartesian refinement during backbone minimization", false);

		// General options
		NEW_OPT( mode, "What mode to run docking?", "legacy");
		NEW_OPT( centroid_silent_out, "Do you want to output final structures as centroid models for assembly?", false);
		NEW_OPT( multi_native, "All pdbs that correspond to the native positions for the current -s pdb", blank_string_vector);
		NEW_OPT( score_natives, "Will dump scored native silent files during point seleciton", false );
		NEW_OPT( dump_inter, "Dump intermediate structures, WARNING!! THIS WILL DUMP A LOT OF PDBS", false );
		NEW_OPT( dump_inter_silent, "Dump intermediate structures as a silent file, WARNING!! THIS WILL DUMP A LOT OF SILENTFILES", false );

		NEW_OPT( final_chain, "The final chain to dump your pose with", "^" );

		// Parallel options
		NEW_OPT( core_idx, "Please enter the core index, and the total cores you're running, ex: '{3} {16}' (nothing implies serial running)", blank_size_vector);
		NEW_OPT( local_result_files, "All of the files that came from the density search step", blank_string_vector);
		NEW_OPT( refine_end, "Search result refinement range start", 0);
		NEW_OPT( refine_start, "Search result refinement range start", 0);
		NEW_OPT( search_end, "What point to end search on", 0);
		NEW_OPT( search_start, "What point to start search on", 0);

		devel::init(argc, argv);

		option[ out::nooutput ].value(true);
		std::chrono::time_point< std::chrono::system_clock > start;
		start = std::chrono::system_clock::now();

		protocols::electron_density::DockPDBIntoDensityMoverOP dock( utility::pointer::make_shared< protocols::electron_density::DockPDBIntoDensityMover >() );
		set_basic_dock_options( dock );
		TR << "Finished setting basic dock options" << std::endl;

		if ( option[ mode ]() == "find_pts" ) {
			get_pts( dock );
		} else if ( option[ mode ]() == "search_points" ) { // TODO: check for overwrite and kill bad shit
			do_search( dock );
		} else if ( option[ mode]() == "combine_search" ) {
			combine_search( dock );
		} else if ( option[ mode ]() == "search_results_to_pdb" ) { // serves as a debugging step
			results_to_pdb( dock );
		} else if ( option[ mode]() == "refine" ) {
			do_refinement( dock );
		} else if ( option[ mode ]() == "combine_refine" ) {
			combine_refinement( dock );
		} else if ( option[ mode ]() == "manual_refine" ) {
			manual_refine( dock );
		} else if ( option[ mode ]() == "cluster_silent" ) {
			cluster_silent( dock );
		} else {
			throw CREATE_EXCEPTION(utility::excn::BadInput, "Didn't recognize your option " + option[ mode ]() + " please change it to one of the following:\n" +
				"find_pts\nearch_points\ncombine_search\nearch_results_to_pdb\nrefine\ncombine_refine");
		}

		auto const end = std::chrono::system_clock::now();
		std::chrono::duration< double > const elapsed_seconds = end-start;
		if ( !option[ testing::INTEGRATION_TEST ].value() ) {
			std::cout << "Took " << elapsed_seconds.count() << " seconds to run " << option[ mode ]() << std::endl;
		}

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
