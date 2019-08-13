// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    calc_nmr_data.cc
/// @brief   Uses ParaNMRScoreMover to score input pose(s) with PCSs, RDCs and PREs
///          and write scores and Q-factors to scorefile. Optionally, the predicted
///          NMR values and infos about the NMR tensor and spinlabel can be written
///          to extra files.
/// @details last Modified: 04/10/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Project headers
#include <protocols/nmr/ParaNMRScoreMover.hh>

// Core headers
#include <devel/init.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobDistributorFactory.hh>
#include <protocols/jd2/internal_util.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/nmr.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

static basic::Tracer TR( "apps.public.nmr.score_with_para_nmr_data" );

int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols;
		using namespace protocols::moves;
		using namespace protocols::jd2;
		using namespace protocols::nmr;
		using namespace basic::options;

		OPT( in::file::silent );
		OPT( in::file::s );
		OPT( in::file::l );
		OPT( in::file::silent_read_through_errors );
		OPT( out::file::silent );
		OPT( out::file::scorefile );
		OPT( nmr::score::verbose );
		OPT( nmr::score::output_exp_calc );
		OPT( nmr::score::output_tensor_info );

		// initialize core
		devel::init(argc, argv);

		TR << " * * * Evaluating pose with Paramagnetic NMR Data * * * " << std::endl;

		// Create scoremover
		ParaNMRScoreMoverOP pnem( new ParaNMRScoreMover );
		bool verbosity(option[ OptionKeys::nmr::score::verbose ]);
		bool calc_vals(option[ OptionKeys::nmr::score::output_exp_calc ]);
		bool tensor_info(option[ OptionKeys::nmr::score::output_tensor_info ]);
		pnem->set_verbose_scorefile(verbosity);
		pnem->write_calc_values(calc_vals);
		pnem->write_tensor_info(tensor_info);

		// Set up a job outputter that writes a scorefile and no PDBs and no Silent Files.
		SilentFileJobOutputterOP jobout( new SilentFileJobOutputter );
		jobout->set_write_no_structures();
		jobout->set_write_separate_scorefile(true);

		JobDistributor::get_instance()->set_job_outputter( JobDistributorFactory::create_job_outputter( jobout ));
		JobDistributor::get_instance()->go( pnem );

		TR << " * * * DONE * * * " << std::endl;

	} catch ( utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;
}
