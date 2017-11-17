// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/doug/ncaa/make_rot_lib.cc
/// @brief Given an input file this app will produce a Dunbrack02 rotlib file
/// @author P. Douglas Renfrew

// core headers
#include <core/types.hh>
#include <basic/options/option.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

// devel headers
#include <devel/init.hh>

// unit headers
#include <protocols/make_rot_lib/RotData.hh>
#include <protocols/make_rot_lib/MakeRotLib.hh>
#include <devel/init.hh>

// utility headers
#include <utility/vector1.hh>

// numeric headers
#include <numeric/angle.functions.hh>

// c++ headers
#include <string>
#include <utility/excn/Exceptions.hh>

using namespace core;
using namespace utility;
using namespace protocols::MakeRotLib;
using namespace scoring;
using namespace basic::options;
using namespace basic::options::OptionKeys;

// application specific options
namespace mrlo {
FileOptionKey const rot_lib_options_file( "mrlo::rot_lib_options_file" );
BooleanOptionKey const peptoid( "mrlo::peptoid" );
RealOptionKey const omega_start_val( "mrlo::omega_start_val" );
RealOptionKey const epsilon_start_val( "mrlo::epsilon_start_val" );
BooleanOptionKey const asp_hack( "mrlo::asp_hack" );
BooleanOptionKey const glu_hack( "mrlo::glu_hack" );
BooleanOptionKey const phe_tyr_hack( "mrlo::phe_tyr_hack" );
}


int
main( int argc, char * argv [] )
{
	try {

	// add application specific options to core options system
	option.add( mrlo::rot_lib_options_file, "Input file for make_rot_lib protocol" );
	option.add( mrlo::peptoid, "Patch NCAA with the patches to make the dipeptoid" ).def( 0 );
	option.add( mrlo::omega_start_val, "Starting value for preceding omega angle in patched dipeptide/dipeptoid" ).def ( 180 );
	option.add( mrlo::epsilon_start_val, "Starting value for epsilon angle in patched dipeptide/dipeptoid" ).def ( 180 );
	option.add( mrlo::asp_hack, "" ).def( false );
	option.add( mrlo::glu_hack, "" ).def( false );
	option.add( mrlo::phe_tyr_hack, "" ).def( false );

	//init
	devel::init(argc, argv);

	// get options file name
	std::string options_filename( option[ mrlo::rot_lib_options_file ]() );

	// create main arrays
	RotVec rotamers, centroids, final_rotamers;
	Size ncluster( 0 );
	std::string aa_name;

	// create score function
	std::cout << "Creating scorefunction..." << std::flush << std::endl;
	ScoreFunctionOP scrfxn( get_score_function_legacy( PRE_TALARIS_2013_STANDARD_WTS ) );
	scrfxn->set_weight( fa_dun, 0.0 );
	scrfxn->set_weight( p_aa_pp, 0.0 );
	scrfxn->set_weight( rama, 0.0 );
	scrfxn->set_weight( fa_intra_rep, 0.0 );
	scrfxn->set_weight( fa_intra_atr, 0.0 );
	scrfxn->set_weight( fa_rep, 0.0 );
	scrfxn->set_weight( fa_atr, 0.0 );
  scrfxn->set_weight( mm_twist, 1.0 );
 	scrfxn->set_weight( mm_lj_inter_rep, 1.0 );
 	scrfxn->set_weight( mm_lj_inter_atr, 1.0 );
 	scrfxn->set_weight( mm_lj_intra_rep, 1.0 );
 	scrfxn->set_weight( mm_lj_intra_atr, 1.0 );

	// initialize rotamers and centroids
	init_rotamers_centroids( rotamers, centroids, ncluster, options_filename, aa_name, option[ mrlo::peptoid ](), option[ mrlo::omega_start_val ](), option[ mrlo::epsilon_start_val ]() );

	// minimize rotamers
	min_rotamers( rotamers, scrfxn, aa_name );

	// hack hack hack, lump symetric angles of symetric AAs together (ie. phe 180, 90 == phe 180, -90)
	if ( option[ mrlo::asp_hack ]() == true ) {
		asp_corrections( rotamers );
	}

	if ( option[ mrlo::glu_hack ]() == true ) {
		glu_corrections( rotamers );
	}

	if ( option[ mrlo::phe_tyr_hack ]() == true ) {
		phe_tyr_corrections( rotamers );
	}

	// seed main loop
	calc_all_dist( rotamers, centroids );
	calc_rotamer_clusters( rotamers );
	calc_centroids( rotamers, centroids );

	// main loop
	//bool rot_change=true;  // unused ~Labonte
	//bool cen_change=true;  // unused ~Labonte
	Size num_iter(0);
	//while ( (rot_change && cen_change) || num_iter <= 500 ) {
	while ( num_iter <= 500 ) {
		++num_iter;
		calc_all_dist( rotamers, centroids  );
		/*rot_change =*/ calc_rotamer_clusters( rotamers );  // unused ~Labonte
		/*cen_change =*/ calc_centroids( rotamers, centroids );  // unused ~Labonte
		std::cout << "ITER:" << num_iter << std::endl;
	}

	std::cout << "----- ROTAMERS-----" << std::endl;
	for(Size i = 1; i<=rotamers.size(); ++i) {
		pretty_print_rd( rotamers[i] );
	}

	//check score of cluster
	std::cout << "AVG_CLSTR_CN_DST:" << avg_cluster_cen_dist( rotamers, ncluster ) << std::endl;

	// pull out best rots from rotamers and add them to final_rotamers
	get_final_rots(rotamers, final_rotamers, ncluster);

	// calc probabilities for final rots
	get_final_rot_probs(final_rotamers);

	std::cout << "----- FINAL CENTROIDS-----" << std::endl;
	for(Size i = 1; i<=centroids.size(); ++i) {
		pretty_print_rd( centroids[i] );
	}

	// get std dev from "energy walk"
	calc_std_dev( final_rotamers, scrfxn, aa_name );

	std::cout << "---- FINAL ROTAMERS -----" << std::endl;
	for(Size i = 1; i<=final_rotamers.size(); ++i) {
		pretty_print_rd( final_rotamers[i] );
	}

	dunbrack_print( final_rotamers, centroids, aa_name );

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

 return 0;
}
