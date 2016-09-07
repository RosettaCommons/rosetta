// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Mike Tyka
/// @brief


// libRosetta headers
#include <protocols/simple_moves/ScoreMover.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/chemical/ChemicalManager.hh>

#include <protocols/cluster/cluster.hh>
#include <protocols/loops/Loops.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/util.hh>

#include <basic/options/option.hh>

#include <devel/init.hh>

// C++ headers
#include <iostream>
#include <string>
#include <deque>

// option key includes

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


using namespace core;
using namespace ObjexxFCL;
using namespace core::pose;
using namespace protocols;
using namespace basic::options;

int
main( int argc, char * argv [] ) {
	try {

		using namespace protocols;
		using namespace protocols::moves;
		using namespace core::scoring;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace utility::file;
		using namespace protocols::cluster;
		using namespace basic::options::OptionKeys::cluster;

		option.add_relevant( OptionKeys::cluster::input_score_filter           );
		option.add_relevant( OptionKeys::cluster::output_score_filter          );
		option.add_relevant( OptionKeys::cluster::exclude_res                  );
		option.add_relevant( OptionKeys::cluster::thinout_factor               );
		option.add_relevant( OptionKeys::cluster::radius                       );
		option.add_relevant( OptionKeys::cluster::limit_cluster_size           );
		option.add_relevant( OptionKeys::cluster::limit_cluster_size_percent   );
		option.add_relevant( OptionKeys::cluster::random_limit_cluster_size_percent   );
		option.add_relevant( OptionKeys::cluster::limit_clusters               );
		option.add_relevant( OptionKeys::cluster::limit_total_structures       );
		option.add_relevant( OptionKeys::cluster::sort_groups_by_energy        );
		option.add_relevant( OptionKeys::cluster::remove_highest_energy_member );
		option.add_relevant( OptionKeys::cluster::limit_dist_matrix            );
		option.add_relevant( OptionKeys::cluster::make_ensemble_cst            );
		simple_moves::ScoreMover::register_options();

		// initialize core
		devel::init(argc, argv);

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << " Rosetta Tool:  cluster - clustering tool for PDBs and or silent files " << std::endl;
		std::cout << " Usage:                                                                  " << std::endl;
		std::cout << "   PDB input:      -in:file:s *.pdb   or  " << std::endl;
		std::cout << "                   -in:file:l  list_of_pdbs  " << std::endl;
		std::cout << "                   -no_optH                                    Dont change positions of Hydrogen atoms! " << std::endl;
		std::cout << "   Silent input:   -in:file:silent silent.out                  silent input filesname " << std::endl;
		std::cout << "                   -in:file:s                                  specify specific tags to be extracted, if left out all will be taken " << std::endl;
		std::cout << "                   -in:file:fullatom                           for full atom structures " << std::endl;
		std::cout << "                   -in:file:silent_struct_type <type>          specify the input silent-file format " << std::endl;
		std::cout << "   Native:         -in:file:native                             native PDB if CaRMS is required " << std::endl;
		std::cout << "   Scorefunction:  -score:weights  weights                     weight set or weights file " << std::endl;
		std::cout << "                   -score:patch  patch                         patch set " << std::endl;
		std::cout << "                   -rescore:verbose                            display score breakdown " << std::endl;
		std::cout << "                   -rescore:output_only                        don't rescore " << std::endl;
		std::cout << "   Output:         -nooutput                                   don't print PDB structures " << std::endl;
		std::cout << "                   -out:prefix  myprefix                       prefix the output structures with a string " << std::endl;
		std::cout << "   Clustering:     -cluster:radius  <float>                    Cluster radius in A (for RMS clustering) or in inverse GDT_TS for GDT clustering. Use \"-1\" to trigger automatic radius detection" << std::endl;
		std::cout << "                   -cluster:gdtmm                              Cluster by gdtmm instead of rms" << std::endl;
		std::cout << "                   -cluster:input_score_filter  <float>        Ignore structures above certain energy " << std::endl;
		std::cout << "                   -cluster:exclude_res <int> [<int> <int> ..] Exclude residue numbers               " << std::endl;
		std::cout << "                   -cluster:radius        <float>              Cluster radius" << std::endl;
		std::cout << "                   -cluster:limit_cluster_size      <int>      Maximal cluster size" << std::endl;
		std::cout << "                   -cluster:limit_cluster_size_percent <float> Maximal cluster size by percentage" << std::endl;
		std::cout << "                   -cluster:random_limit_cluster_size_percent <float> Maximal cluster size by percentage, cut randomly" << std::endl;
		std::cout << "                   -cluster:limit_clusters          <int>      Maximal number of clusters" << std::endl;
		std::cout << "                   -cluster:limit_total_structures  <int>      Maximal number of structures in total" << std::endl;
		std::cout << "                   -cluster:sort_groups_by_energy              Sort clusters by energy." << std::endl;
		std::cout << "                   -cluster:remove_highest_energy_member       Remove highest energy member of each cluster" << std::endl;
		std::cout << "                   -symmetry:symmetric_rmsd       \t\t\t\t\t\t For symmetric systems find the lowest rms by testing all chain combinations. Works only with silent file input that contain symmetry info and with all CA rmsd" << std::endl;
		std::cout << " Examples: " << std::endl;
		std::cout << "   cluster -database ~/minirosetta_database -in:file:silent silent.out -in::file::binary_silentfile -in::file::fullatom -native 1a19.pdb " << std::endl;
		std::cout << "clustered Poses are given output names in the form of:" << std::endl;
		std::cout << " c.i.j, which denotes the jth member of the ith cluster." << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;

		if ( !option[ out::output ].user() ) {
			option[ out::nooutput ].value( true );
		}

		int time_start = time(nullptr);

		MoverOP mover;
		core::scoring::ScoreFunctionOP sfxn;
		sfxn = core::scoring::get_score_function();
		if ( option[ basic::options::OptionKeys::symmetry::symmetric_rmsd ]() ) {
			core::scoring::ScoreFunctionOP sfxn_sym =
				core::scoring::symmetry::symmetrize_scorefunction( *sfxn );
			sfxn = sfxn_sym;
		}

		ClusterPhilStyleOP clustering;

		if ( option[ basic::options::OptionKeys::cluster::loops ]() ) {
			loops::Loops loops( true );
			clustering = ClusterPhilStyleOP( new ClusterPhilStyle_Loop(loops ) );
		} else {
			clustering = ClusterPhilStyleOP( new ClusterPhilStyle() );
		}

		clustering->set_score_function( sfxn );
		if ( option[ basic::options::OptionKeys::cluster::input_score_filter ].user() ) {
			clustering->set_filter( option[ basic::options::OptionKeys::cluster::input_score_filter ] );
		}
		clustering->set_cluster_radius(
			option[ basic::options::OptionKeys::cluster::radius ]()
		);
		clustering->set_population_weight(
			option[ basic::options::OptionKeys::cluster::population_weight ]()
		);


		// Figure out the thinout factor,
		//mjo commenting out 'thinout_factor' because it is unused and causes a warning
		//core::Real  thinout_factor = 0.0;

		/// std::cerr << "STAARTING to cluster" << std::endl;
		/// if( option[ basic::options::OptionKeys::cluster::thinout_factor ].user() ){
		///  thinout_factor = option[ basic::options::OptionKeys::cluster::thinout_factor ]();
		/// } else {
		///  // autodetermine a good thinout factor
		///  core::Size max_O2_clustering = 400;
		///  core::Size expected_n_decoys = n_input_structures();
		///  // figure out the number of decoys:
		///  thinout_factor = 1.0 - core::Real( max_O2_clustering ) / core::Real( expected_n_decoys );
		///  if( thinout_factor < 0 ) thinout_factor = 0.0; // few decoys, no thinout required
		///  std::cout << "Estimated thinout factor to seed with " << max_O2_clustering << " structures: " << thinout_factor << " ( " << expected_n_decoys << " ) " << std::endl;
		/// }

		core::chemical::ResidueTypeSetCOP rsd_set;
		if ( option[ in::file::fullatom ]() ) {
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		} else {
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "centroid" );
		}


		// Cluster the first up-to-400 structures by calculating a full rms matrix
		core::import_pose::pose_stream::MetaPoseInputStream input = core::import_pose::pose_stream::streams_from_cmd_line();
		while ( input.has_another_pose() && (clustering->nposes() < 400 ) ) {
			core::pose::Pose pose;
			input.fill_pose( pose, *rsd_set );
			clustering->apply( pose );
		}

		mover = clustering;

		int time_readin = time(nullptr);
		clustering->do_clustering( option[ OptionKeys::cluster::max_total_cluster ]() );
		int time_initialc = time(nullptr);
		clustering->do_redistribution();

		// Process any remaining structures by asigning to clusters or forming new clusters
		std::cout << "Assigning extra structures ... " << std::endl;
		AssignToClustersMoverOP mover_add_structures( new AssignToClustersMover( clustering ) );
		mover_add_structures->set_score_function( sfxn );
		//mover_add_structures->set_cluster_radius( clustering.get_cluster_radius() );

		while ( input.has_another_pose() ) {
			core::pose::Pose pose;
			input.fill_pose( pose, *rsd_set );
			mover_add_structures->apply( pose );
		}


		clustering->print_summary();

		// Post processing
		int time_total = time(nullptr);
		// clustering->sort_each_group_by_energy();
		if ( option[ sort_groups_by_energy ] ) {
			clustering->sort_groups_by_energy();
		}
		if ( option[ remove_singletons ] ) {
			clustering->remove_singletons();
		}
		if ( option[ limit_cluster_size ].user() ) {
			clustering->limit_groupsize( option[ limit_cluster_size ] );
		}
		if ( option[ limit_cluster_size_percent ].user() ) {
			clustering->limit_groupsize( option[ limit_cluster_size_percent ] );
		}
		if ( option[ random_limit_cluster_size_percent ].user() ) {
			clustering->random_limit_groupsize( option[ random_limit_cluster_size_percent ] );
		}
		if ( option[ limit_clusters ].user() ) {
			clustering->limit_groups( option[ limit_clusters ] );
		}
		if ( option[ limit_total_structures ].user() ) {
			clustering->limit_total_structures( option[ limit_total_structures] );
		}
		if ( option[ remove_highest_energy_member ] ) {
			clustering->remove_highest_energy_member_of_each_group();
		}

		if ( option[ export_only_low ] ) {
			clustering->sort_each_group_by_energy();
			clustering->sort_groups_by_energy( );
			clustering->export_only_low( option[ export_only_low ]() );
		}
		// --------------------------------------------------------------------
		// Results:
		clustering->print_summary();
		if ( option[ basic::options::OptionKeys::out::file::silent ].user() ) {
			clustering->print_clusters_silentfile( option[ out::prefix ]() );
		} else {
			clustering->print_cluster_PDBs( option[ out::prefix ]() );
		}

		if ( option[ basic::options::OptionKeys::cluster::make_ensemble_cst]() ) {
			EnsembleConstraints_Simple cec( 1.0 );
			clustering->create_constraints( option[ out::prefix ](), cec );
		}

		clustering->print_cluster_assignment();
		std::vector < Cluster >  const & clusterlist=clustering->get_cluster_list();
		std::list < int > sorted_list;
		for (const auto & i : clusterlist) {
			for ( int j=0; j<(int)i.size(); j++ ) {
				sorted_list.push_back(  i[j]  );
			}
		}

		sorted_list.sort();

		std::cout << "Timing: " << std::endl;
		std::cout <<    "  Readin:" << time_readin - time_start
			<< "s\n  Cluster: "    << time_initialc - time_readin
			<< "s\n  Additional Clustering: " << time_total - time_initialc
			<< "s\n  Total: " << time_total - time_start
			<< std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

