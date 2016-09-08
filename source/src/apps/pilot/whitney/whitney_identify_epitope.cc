// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// @author jk

#include <iostream>
#include <iomanip>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>

#include <core/scoring/rms_util.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using namespace core;
using namespace core::scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_KEY( String, identify_epitope_chain )
OPT_KEY( String, constant_score_chain )
OPT_KEY( String, native_pdb )

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.whitney_identify_epitope.main" );

void add_scores (core::pose::Pose & pose,
								 utility::vector1 <core::Real> scores,
								 core::Real kT,
								 char const epitope_chain,
								 char const const_chain ) {

	// find the min and max weights from the protein
	core::Real min = 99999;
	core::Real max = -9999;

	Size nres = pose.size();


	for ( Size i = 1; i <= nres; ++i ) {
		if ( scores.at(i) < min ) min = scores.at(i);
		if ( scores.at(i) > max ) max = scores.at(i);
	}

	TR << "Max score " << max << std::endl;
	TR << "Min score " << min << std::endl;

	//calculate average score
	core::Real average_score =  ( ( max - min ) / 2) * (99.99 / max );

	for ( Size i=1; i <= nres; ++i ) {

		//calculate new score
		core::Real new_score = (( scores.at(i) - min) * ( 99.99 /max ));

		for ( Size j=1; j <= pose.residue( i ).natoms(); j++ ) {

			if ( pose.pdb_info()->chain(i) == const_chain ) {

				pose.pdb_info()->temperature(i, j, average_score );

			}

			if ( pose.pdb_info()->chain(i) == epitope_chain ) {

				pose.pdb_info()->temperature( i, j, new_score );
			}

		}
	}

	//dump pdb file with scores
	std::ostringstream outPDB_name;
	outPDB_name << "identify_epitope" << kT << ".pdb";

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(outPDB_name.str(), std::ios::out);
	pose.dump_pdb( outPDB_stream );
	outPDB_stream.close();
	outPDB_stream.clear();

}


/// General testing code
int
main( int argc, char * argv [] )
{

	try {


  NEW_OPT( identify_epitope_chain, "chain to calculate scores for", "B" );
	NEW_OPT( constant_score_chain, "chain will have the avg score", "A" );
	NEW_OPT( native_pdb, "chain with original pdb", "native.pdb" );

	devel::init(argc, argv);

	std::string tmp_chain = option[ identify_epitope_chain ];
	std::string tmp2_chain = option[ constant_score_chain ];

	std::string const native_pdb_name = option[ native_pdb ] ;

	if ( tmp_chain.length() != 1  || tmp2_chain.length() != 1 ) {
	scoring::ScoreFunctionOP scorefxn( get_score_function() );

	std::cerr << "ERROR!! Chain ID should be one character!" << std::endl;
	exit(1);
	}
	// create poses from native pdb for different kT
	pose::Pose native_pose, first_pose;
	core::import_pose::pose_from_file( native_pose, native_pdb_name , core::import_pose::PDB_file);
	core::import_pose::pose_from_file( first_pose, native_pdb_name , core::import_pose::PDB_file);

	TR << "Starting calculating scores for residues" << std::endl;

	// set options from command line
	char const epitope_chain = tmp_chain[0];
	char const const_chain = tmp2_chain[0];

	//kT values
	core::Real const low_kT = 0.2;
	core::Real const avg_kT = 0.59;
	core::Real const high_kT = 1000;

	core::Size const res_num = native_pose.size();

	// create and open a file for output of scores
	//utility::io::ozstream score_outstream;
	//score_outstream.open( "residue_scores.out", std::ios::out );

	//ws debug
	TR << "Processing Energy minimum" << std::endl;

  core::Real energy_min = 9999999;
	// loop to calculate the minimum score for the list of decoys
 	for (core::Size f=1; f <=  basic::options::start_files().size(); f++) {

		std::string const curr_decoy_fname = basic::options::start_files().at(f);

		// ws debug
		TR << "Current decoy " << curr_decoy_fname << std::endl;

		// scoring function
		scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( "interchain_cen" ) );

		// calculate the score for the current decoy
		pose::Pose curr_pose;
		core::import_pose::centroid_pose_from_pdb( curr_pose, curr_decoy_fname , core::import_pose::PDB_file);
		(*scorefxn)(curr_pose);
		core::Real curr_energy = curr_pose.energies().total_energies()[ total_score ];


		if ( curr_energy < energy_min ) energy_min = curr_energy;

	}

	TR << "The minimum energy for the list of decoys is: " << energy_min << std:: endl;

	// vectors for all scores calculated with different kT values
	utility::vector1 <core::Real> low_kT_scores;
	utility::vector1 <core::Real> avg_kT_scores;
	utility::vector1 <core::Real> high_kT_scores;

	low_kT_scores.resize(res_num, 0 );
	avg_kT_scores.resize(res_num, 0 );
	high_kT_scores.resize(res_num, 0 );

	// calculate the score for each decoy and add that score to each residue
	for (core::Size d=1; d <= basic::options::start_files().size(); d++) {

		std::string const curr_decoy_fname = basic::options::start_files().at(d);

		// scoring function
		scoring::ScoreFunctionOP scorefxn( ScoreFunctionFactory::create_score_function( "interchain_cen" ) );

		// score current decoy
		pose::Pose curr_pose;
		core::import_pose::centroid_pose_from_pdb( curr_pose, curr_decoy_fname , core::import_pose::PDB_file);
		(*scorefxn)(curr_pose);
		core::Real curr_energy = curr_pose.energies().total_energies()[ total_score ];

		// calculate score using Boltzmann's Principle
		core::Real energy_diff = curr_energy - energy_min;
		core::Real low_kTscore =  std::exp( - ( energy_diff / low_kT ) );
		core::Real avg_kTscore =  std::exp( - ( energy_diff / avg_kT ) );
		core::Real high_kTscore =  std::exp( - ( energy_diff / high_kT ) );

		// vector for residues in the interface
		utility::vector1 <bool> interface;
		interface.resize( curr_pose.size(), false );

		// define interface
		core::Real interface_dist = 8.0;
		core::Size rb_jump = 1;

		pack::task::TaskFactory tf;
		tf.push_back( new protocols::toolbox::task_operations::RestrictToInterface( rb_jump, interface_dist ) );
		pack::task::PackerTaskOP task = tf.create_task_and_apply_taskoperations( curr_pose );
		for ( core::Size i=1; i <= curr_pose.size(); ++i ) {
			if ( task->pack_residue(i) && curr_pose.pdb_info()->chain(i) == epitope_chain ) interface.at(i)=true;
		}

		// add scores at interface residues
		for ( Size i=1; i <= curr_pose.size(); ++i ) {

			if ( interface.at(i) ) {

				low_kT_scores.at(i) = low_kT_scores.at(i) + low_kTscore;
				avg_kT_scores.at(i) = avg_kT_scores.at(i) + avg_kTscore;
				high_kT_scores.at(i) = high_kT_scores.at(i) + high_kTscore;

			}
		}

	}

	// calculate scores and dump pdb
	add_scores( native_pose, low_kT_scores, low_kT, epitope_chain, const_chain );

	native_pose = first_pose;
	add_scores( native_pose, avg_kT_scores, avg_kT, epitope_chain, const_chain );

	native_pose = first_pose;
	add_scores( native_pose, high_kT_scores, high_kT, epitope_chain, const_chain );

	TR << "Done computing residue scores" << std::endl;

	return 0;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

