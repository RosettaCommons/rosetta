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
/// @authors Rhiju Das

// libRosetta headers
#include <core/types.hh>
#include <core/sequence/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/init/init.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>

//RNA stuff -- move this?
#include <protocols/stepwise/enumerate/general/StepWiseClusterer.hh>

// C++ headers
#include <iostream>
#include <string>


// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>

//Auto Headers
#include <core/chemical/ResidueType.hh>
#include <utility/excn/Exceptions.hh>

using namespace protocols;
using namespace core;
using utility::vector1;


//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
void
cluster_test(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::cluster;


	if ( ! option[ in::file::silent ].user() ) utility_exit_with_message( "The rna_cluster executable requires silent input [with -in:file:silent], and models need to be scored. The clustering algorithm starts with the lowest scoring models and works its way up." );

	utility::vector1< std::string > const silent_files_in( option[ in::file::silent ]() );
	protocols::stepwise::enumerate::general::StepWiseClusterer stepwise_clusterer( silent_files_in );

	Size max_decoys( 400 );
	if ( option[ out::nstruct].user() )	 max_decoys =  option[ out::nstruct ];
	stepwise_clusterer.set_max_decoys( max_decoys );


	Real cluster_radius ( 2.0 );
	if ( option[ OptionKeys::cluster::radius ].user() ) cluster_radius = option[ OptionKeys::cluster::radius ]();
	stepwise_clusterer.set_cluster_radius(	cluster_radius	);

	stepwise_clusterer.set_cluster_by_all_atom_rmsd( true );
	stepwise_clusterer.set_score_diff_cut( option[ score_diff_cut ] );
	stepwise_clusterer.set_rename_tags( true /*option[ rename_tags ]*/ );
	stepwise_clusterer.set_rsd_type_set( "rna" );
	stepwise_clusterer.set_auto_tune( option[ auto_tune ] );

	stepwise_clusterer.cluster();

	std::string const silent_file_out( option[ out::file::silent  ]() );
	stepwise_clusterer.output_silent_file( silent_file_out );

	std::cout << "Maximum number of models: " << max_decoys << ". [use -nstruct to set higher]" << std::endl;
	std::cout << "Used all-heavy-atom clustering radius of " << cluster_radius << " Angstroms. [use -cluster:radius to change.]" << std::endl;

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	cluster_test();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
try {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::cout << std::endl << "Basic usage:  " << argv[0] << " -in:file:silent <input silent file> -out:file:silent <output silent file> -cluster:radius <RMSD threshold in Angstroms>" << std::endl;
	std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

	option.add_relevant(  in::file::silent );
	option.add_relevant(  out::file::silent );
	option.add_relevant(  out::nstruct );
	option.add_relevant(  cluster::radius );
	option.add_relevant(  cluster::score_diff_cut );
	option.add_relevant(  cluster::auto_tune );


	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	core::init::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );
} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
}
}
