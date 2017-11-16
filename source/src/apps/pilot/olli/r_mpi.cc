// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief  check quality of fragments against input structure
/// @author Oliver Lange

// #include <protocols/abinitio/Templates.hh>
// #include <protocols/abinitio/TemplateJumpSetup.hh>
// #include <protocols/abinitio/PairingStatistics.hh>
// #include <protocols/abinitio/StrandConstraints.hh>
// #include <protocols/moves/Mover.hh>

#include <core/pose/Pose.fwd.hh>
//#include <core/pose/util.hh>
#include <devel/init.hh>
//#include <core/io/pdb/pdb_writer.hh>

#include <core/scoring/ResidualDipolarCoupling.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
//#include <core/scoring/constraints/BoundConstraint.hh>
#include <protocols/jd2/NoOutputJobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/archive/ArchiveManager.hh>

//for derivative check
#include <protocols/jd2/MpiFileBuffer.hh>

// Utility headers
#include <basic/options/option_macros.hh>
#include <utility/io/ozstream.hh>
//#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/exit.hh>
#include <utility/excn/Exceptions.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/io/silent/SilentFileData.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

static basic::Tracer tr( "main" );

using namespace core;
using namespace protocols;
using namespace protocols::jd2;
//using namespace abinitio;
//using namespace jumping;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace scoring;
//using namespace scoring::constraints;
using namespace utility::io;

OPT_1GRP_KEY( Boolean, score_app, linmin )
OPT_KEY( Boolean, dump_all )
OPT_KEY( File, average )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//  Templates::register_options();
	OPT( in::file::s );
	NEW_OPT( score_app::linmin, "Do a quick round of linmin before reporting the score", false );
	NEW_OPT( dump_all, "show individual RDCs for each structure on screen", false );
	NEW_OPT( average, "compute average RDC and mean_deviation for each residue", "average.dat" );
}

void run() {
	protocols::jd2::archive::ArchiveManager archive( 0,0,0 );
	archive.go();
	//   WriteOut_MpiFileBuffer buffer( 0 );
	//   std::cout << "hello\n" << std::endl;
	//  buffer.run();
	//  ozstream out;
	//  out.open("bla");
	//  for ( Size i = 1; i<10; i++ ) {
	//   out << A( 25, "shello") << i << F( 20, 5, 20.42 ) << std::endl;
	//   for ( Size j = 1; j<20000000; j++ ) {
	//    Real k = exp( j );
	//   }
	//  }
	//  std::cerr<<"loop-finished" << std::endl;
	//  using namespace basic::options;
	//  using namespace basic::options::OptionKeys;
	//  if ( option[ in::file::s ].user() ) {
	//   std::cerr << "here" << std::endl;
	//   core::pose::Pose pose;
	//   core::import_pose::pose_from_file( pose, "S_2_0091.pdb" , core::import_pose::PDB_file);
	//   core::io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_silent_struct_out();
	//   pss->fill_struct( pose, "mala" );
	//    std::cerr << "start output..." << std::endl;
	//    pss->print_header( out );
	//    std::cerr << "header completed" << std::endl;
	//    pss->print_scores( out );
	//    std::cerr << "scores completed" << std::endl;
	//    pss->print_conformation( out );
	//    std::cerr << "conformation completed" << std::endl;
	//   out << pss << std::endl;
	//   out.close();

	//   core::io::silent::SilentFileData sfd;
	//   sfd.add_structure( pss );
	//   sfd.write_all( "new_file" );
	//  }
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/// =============================== MAIN ============================================================
/////////////////////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
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


// Forward
