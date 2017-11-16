// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/public/cst_info.cc
/// @brief Output details about constraint satisfaction of a Pose.
/// @author Rocco Moretti (rmorettiase@gmail.com)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/simple_moves/CstInfoMover.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/moves/MoverContainer.hh>

// core headers
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>

// utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>

static basic::Tracer TR("cst_info");

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( in::file::s );
	option.add_relevant( in::file::l );
	option.add_relevant( in::file::silent );
	option.add_relevant( constraints::cst_file );
	option.add_relevant( constraints::cst_fa_file );

}

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols;
		using namespace protocols::jd2;
		using namespace protocols::simple_moves;

		devel::init( argc, argv );

		if ( ! option[ out::no_nstruct_label ] && ! option[ out::nstruct ].user() ) {
			option[ out::no_nstruct_label ].value( true ); // Don't add the nstruct label unless the user requests it.
		}
		register_options();

		utility::vector1< std::string > cstfiles;
		if ( option[ constraints::cst_file ].user() ) {
			cstfiles.append( option[ constraints::cst_file ]() );
		}
		if ( option[ constraints::cst_fa_file ].user() ) {
			cstfiles.append( option[ constraints::cst_fa_file ]() );
		}

		// Setup mover here

		// (I would have liked to set things up such that the constraint definitions are output only once,
		// in the initial setup phase, but as constraints require a pose to be read in, that's a non-starter.)

		moves::SequenceMoverOP protocol( new protocols::moves::SequenceMover );
		if ( cstfiles.size() < 1  ) {
			utility_exit_with_message("Must specify either -constraints:cst_file or -constraints:cst_fa_file for the cst_info application!");
		} else {
			for ( core::Size ii(1); ii <= cstfiles.size(); ++ii ) {
				CstInfoMoverOP cstinfomover( new protocols::simple_moves::CstInfoMover() );
				cstinfomover->cst_file( cstfiles[ii] );
				cstinfomover->recursive( true );
				std::string prefix( "CstFile"+ utility::to_string( ii ) );
				cstinfomover->prefix( prefix );
				TR << "Loading constraint file " << cstfiles[ii] << " as " << prefix << std::endl;
				protocol->add_mover( cstinfomover );
			}
			simple_moves::ScoreMoverOP rescore( new simple_moves::ScoreMover("empty") );
			rescore->set_verbose( false );
			protocol->add_mover( rescore );
		}

		// Set up a job outputter that writes a scorefile and no PDBs and no Silent Files.
		SilentFileJobOutputterOP jobout( new SilentFileJobOutputter );
		jobout->set_write_no_structures();
		jobout->set_write_separate_scorefile(true);

		// If the user chooses something else, then so be it, but by default only create a score file and nothing else.
		JobDistributor::get_instance()->set_job_outputter( JobDistributorFactory::create_job_outputter( jobout ) );

		JobDistributor::get_instance()->go( protocol );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
