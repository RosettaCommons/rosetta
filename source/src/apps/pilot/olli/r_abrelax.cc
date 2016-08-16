// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#include <protocols/moves/Mover.hh>
#include <basic/options/option.hh>
#include <core/types.hh>
#include <devel/init.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/abinitio/IterativeAbrelax.hh>
#include <protocols/jd2/archive/ArchiveManager.hh>
#include <protocols/abinitio/BrokerMain.hh>

#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/jumps.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/broker.OptionKeys.gen.hh>
#include <basic/MemTracer.hh>
#include <utility/excn/Exceptions.hh>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>

void register_options_broker() {
  using namespace protocols;
  using namespace abinitio;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace core;
  using namespace fragment;

  option.add_relevant( broker::setup );

  option.add_relevant(in::file::native);
  option.add_relevant(in::file::silent);
  option.add_relevant(in::file::frag3);
  option.add_relevant(in::file::frag9);
  option.add_relevant(in::file::fasta);

  option.add_relevant(out::file::silent);
  option.add_relevant(out::nstruct);

  option.add_relevant(run::test_cycles);
  // constraints
  option.add_relevant(constraints::cst_file);

  option.add_relevant(frags::nr_large_copies);
  option.add_relevant(frags::annotate);

  protocols::abinitio::AbrelaxApplication::register_options();
  protocols::abinitio::IterativeAbrelax::register_options();
  protocols::jd2::archive::ArchiveManager::register_options();

}


int
main( int argc, char * argv [] )
{
	try{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using std::string;
  using utility::vector1;

	// basic::get_usage_from_procfilesystem( std::cerr );
// 	std::cerr << " @ program start" << std::endl;

	int rank( 0 );
#ifdef USEMPI
	{ // scope
	int already_initialized( 0 );
	MPI_Initialized( & already_initialized );
	if ( already_initialized == 0 ) MPI_Init(&argc, &argv);
	MPI_Comm_rank( MPI_COMM_WORLD, ( int* )( &rank ) );
	}
#endif
	if ( rank == 0 ) {
		basic::get_usage_from_procfilesystem( std::cerr );
		std::cerr << " @ program start" << std::endl;
	}

	register_options_broker();

	if ( rank == 0 ) {
		basic::get_usage_from_procfilesystem( std::cerr );
		std::cerr << " @ register_options" << std::endl;
	}

  try {
    devel::init( argc, argv );//
		basic::mem_tr << "devel::init" << std::endl;
		//    protocols::abinitio::Broker_main();
		protocols::jd2::JobDistributor::get_instance()->go( NULL );
  } catch ( utility::excn::EXCN_Base& excn ) {
    std::cerr << "Exception : " << std::endl;
    excn.show( std::cerr );
  }
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
  return 0;
}
