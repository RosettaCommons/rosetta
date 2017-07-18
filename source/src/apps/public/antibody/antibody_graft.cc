// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/pilot/jianqing/antibody_graft.cc
/// @brief Graft antibody framework and CDR loops.  Use with Rosetta/tools/antibody/antibody.py to select templates.
/// @author Jianqing Xu (xubest@gmail.com), Daisuke Kuroda, Jeff Gray
/// 09/09/2011 - 2013


#include <protocols/antibody/GraftCDRLoopsProtocol.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

#include <devel/init.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
//#include <utility/tools/make_vector1.hh>

// option key includes
#include <basic/options/option.hh>
#include <string>
#include <basic/Tracer.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.antibody" );


int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace protocols::antibody;
		using namespace protocols::jd2;

		GraftCDRLoopsProtocol::register_options();
		protocols::jd2::register_options();

		devel::init(argc, argv);

		GraftCDRLoopsProtocolOP abm( new GraftCDRLoopsProtocol() );

		TR<< *abm << std::endl;  // exit(-1);

		JobDistributor::get_instance()->go(abm);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
