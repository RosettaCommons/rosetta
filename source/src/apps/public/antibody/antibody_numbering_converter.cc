// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/public/antibody/antibody_numbering_converter.cc
/// @brief App to convert antibodies to a different numbering scheme.  Antibodies must have L and/or H chains.  Only works on a FAB, currently, just like the rest of RosettaAntibody.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// devel headers
#include <devel/init.hh>

// protocol headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/antibody/AntibodyNumberingConverterMover.hh>


// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>

static basic::Tracer TR("antibody_numbering_converter");


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

		devel::init( argc, argv );
		register_options();




		if ( ( ! option [ in::file::l ].user() ) && ( ! option [ in::file::s ].user() ) ) {
			utility_exit_with_message("Please specify either -s or -l to specify the input PDB.");
		}


		protocols::antibody::AntibodyNumberingConverterMoverOP mover_protocol( new protocols::antibody::AntibodyNumberingConverterMover() );

		protocols::jd2::JobDistributor::get_instance()->go( mover_protocol );


	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
