// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  analyse sets of structures
/// @details This tool allows to superimpose structures using the wRMSD method [ Damm&Carlson, Biophys J (2006) 90:4558-4573 ]
/// @details Superimposed structures can be written as output pdbs and the converged residues can be determined
/// @author Oliver Lange

#include <devel/init.hh>
#include <core/types.hh>
#include <protocols/toolbox/DecoySetEvaluation.hh>
#include <protocols/toolbox/Cluster.hh>
#include <protocols/toolbox/Cluster.impl.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/LoopsFileIO.hh>

// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/io/ozstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>

// ObjexxFCL includes
#include <ObjexxFCL/FArray1D.hh>

#include <string>

//Auto Headers
#include <utility/excn/EXCN_Base.hh>


static basic::Tracer tr( "main" );

using namespace core;
using namespace protocols;
//using namespace pose;
using namespace protocols::toolbox;
using namespace ObjexxFCL; //Farray
//using namespace ObjexxFCL::format;

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core::io::silent;

OPT_2GRP_KEY( File, out, file, cluster )
OPT_1GRP_KEY( File, rigid, in )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	OPT( in::file::silent );
	OPT( in::file::native );
	OPT( cluster::limit_cluster_size );
	NEW_OPT( out::file::cluster, "write clustered structures to silent file with this name", "cluster.out" );
	NEW_OPT( rigid::in, "residues that are considered for clustering", "rigid.loop");
}

void read_input_weights( FArray1D_double& weights, Size natoms ) {
	if ( !option[ rigid::in ].user() ) return;
	loops::PoseNumberedLoopFileReader loop_file_reader;
	loop_file_reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );
	std::ifstream is( option[ rigid::in ]().name().c_str() );
	if ( !is ) utility_exit_with_message( "[ERROR] Error opening RBSeg file '" + option[ rigid::in ]().name() + "'" );
	loops::SerializedLoopList loops = loop_file_reader.read_pose_numbered_loops_file(is, option[ rigid::in ](), false );
	loops::Loops rigid = loops::Loops( loops );
	for ( Size i=1; i<=natoms; ++i ) {
		if ( rigid.is_loop_residue( i ) ) weights( i )=1.0;
		else weights( i )=0.0;
	}
}

void read_structures( SilentFileData &sfd, DecoySetEvaluation& ensemble ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( !option[ in::file::silent ].user() ) {
		utility_exit_with_message("fast_clustering works only with -in:file:silent as input option");
	}
	utility::vector1< utility::file::FileName > const silent_files( option[ in::file::silent ]() );
	for ( utility::vector1< utility::file::FileName >::const_iterator current_fn_ = silent_files.begin(); current_fn_ != silent_files.end(); ++current_fn_ ) {
		tr.Debug << "reading " << *current_fn_ << std::endl;
		sfd.read_file( *current_fn_ );
	}
	ensemble.push_back_CA_xyz_from_silent_file( sfd, true /*store energies*/ );
}

void run() {
	//bool store_energies ( true );  // unused ~Labonte

	DecoySetEvaluation ensemble;
	SilentFileOptions opts; // initialized from the command line
	SilentFileData sfd( opts );
	read_structures( sfd, ensemble );

	//initialize wRMSD weights with 1.0 unless we have -rigid:in file active
	FArray1D_double weights( ensemble.n_atoms() , 1.0 );
	read_input_weights( weights, ensemble.n_atoms() );
	ensemble.set_weights( weights );

	SilentFileData kept_decoys( opts );
	cluster_silent_structs( ensemble, sfd.begin(), sfd.end(), kept_decoys, ClusterOptions( true ) );
	std::string out_filename=option[ out::file::cluster ]();
	{ utility::io::ozstream out( out_filename ); } //open and close, to empty the file.
	kept_decoys.write_all( out_filename );

	return;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// =============================== MAIN ============================================================
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		register_options();
		devel::init( argc, argv );

		try{
			run();
		} catch ( utility::excn::EXCN_Base& excn ) {
			excn.show( std::cerr );
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}


