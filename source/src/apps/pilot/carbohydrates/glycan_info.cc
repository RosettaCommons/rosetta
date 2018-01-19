// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/carbohydrates/glycan_relax.cc
/// @brief Application for Glycan Relax.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) and Jason W. Labonte (JWLabonte@gmail.com)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>
#include <protocols/analysis/GlycanInfoMover.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <utility/string_util.hh>

static basic::Tracer TR("glycan_info");


void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	option.add_relevant( in::file::s );
	option.add_relevant( in::file::l );

}


int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace protocols::jd2;

		devel::init( argc, argv );
		//register_options();

		if ( ( ! option [ in::file::l ].user() ) && ( ! option [ in::file::s ].user() ) ) {
			utility_exit_with_message("Please specify either -s or -l to specify the input PDB.");
		}



		// Make sure the default JobOutputter is SilentJobOutputter to ensure that when this
		// is called with default arguments is prints a proper scorefile and not the hacky thing that
		// the  JobOutputter scorefile() function produces (which for example skips Evaluators!!)

		// Set up a job outputter that writes a scorefile and no PDBs and no Silent Files.
		// Copied from score_jd2.
		SilentFileJobOutputterOP jobout( new SilentFileJobOutputter );
		jobout->set_write_no_structures();
		jobout->set_write_separate_scorefile(true);


		protocols::jd2::JobDistributor::get_instance()->set_job_outputter( JobDistributorFactory::create_job_outputter( jobout ));


		protocols::analysis::GlycanInfoMoverOP mover_protocol( new protocols::analysis::GlycanInfoMover() );

		protocols::jd2::JobDistributor::get_instance()->go( mover_protocol );


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
