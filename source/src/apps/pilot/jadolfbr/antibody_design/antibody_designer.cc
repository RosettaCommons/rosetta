// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/antibody_designer.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/antibody/design/AntibodyDesignProtocol.hh>
#include <protocols/viewer/viewers.hh>

#include <utility/excn/Exceptions.hh>
#include <devel/init.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>

using namespace protocols::antibody::design;
using namespace basic::options;

// FUNCTION DECLARATIONS
void*
my_main( void *)
{
	protocols::jd2::JobDistributor::get_instance()->go(protocols::moves::MoverOP( new AntibodyDesignProtocol ));
	return 0 ;
}

//Description:  This application will become the Rosetta Antibody Designer.  Main code is handled through AntibodyDesignProtocol
int main(int argc, char* argv[])
{
	try{
		devel::init(argc, argv);

		// Check to make sure we pass a PDB.  This may need to change for silent files, etc. etc., but there are no good ways to exit if this is not set.
		if ( ( ! option [ OptionKeys::in::file::l ].user() ) && ( ! option [ OptionKeys::in::file::s ].user() ) ) {
			utility_exit_with_message("Please specify either -s or -l to specify the input PDB.");
		}


		if ( option [ OptionKeys::antibody::design::view ]() ) {
			protocols::viewer::viewer_main( my_main );
		}
		else {
			protocols::jd2::JobDistributor::get_instance()->go(protocols::moves::MoverOP( new AntibodyDesignProtocol ));
		}

	} catch( utility::excn::EXCN_Base& excn ) {
		basic::Error()
			<< "ERROR: Exception caught by antibody_design application:"
			<< excn << std::endl;
		std::exit( 1 );
	}

}
