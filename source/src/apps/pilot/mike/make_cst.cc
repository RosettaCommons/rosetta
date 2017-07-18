// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file Mike Tyka
/// @brief


// libRosetta headers
#include <protocols/simple_moves/ScoreMover.hh>
#include <core/chemical/ChemicalManager.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/NullMover.hh>
#include <protocols/simple_moves/MinMover.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility/file/FileName.hh>
#include <utility/excn/Exceptions.hh>
#include <protocols/cluster/cluster.hh>
#include <devel/init.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>


using namespace core;
using namespace basic::options;


int
main( int argc, char * argv [] )
{
	try {
		using namespace protocols;
		using namespace protocols::jobdist;
		using namespace protocols::moves;
		using namespace core::scoring;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace utility::file;

		// scoring should by default not produce output files - that's so annoying
		// unless of coures the user insists.

		// initialize core
		devel::init(argc, argv);
		jd2::register_options();
		if ( !option[ out::output ].user() ) {
			option[ out::nooutput ].value( true );
		}

		protocols::cluster::EnsembleConstraints_Simple *cec = new protocols::cluster::EnsembleConstraints_Simple( 1.0 );
		MoverOP mover = cec;
		protocols::jd2::JobDistributor::get_instance()->go( mover );
		cec->createConstraints( std::cout );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

