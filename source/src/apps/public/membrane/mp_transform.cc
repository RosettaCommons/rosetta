// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/membrane/mp_transform.cc
/// @brief  RosettaMP Transform protein into a membrane
/// @author  Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @note   Last Updated: 10/14/15

// App headers
#include <devel/init.hh>

// Project Headers
#include <core/conformation/membrane/MembraneInfo.hh>


#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/OptimizeProteinEmbeddingMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/visualize/VisualizeEmbeddingMover.hh>

// Package Headers
#include <apps/benchmark/performance/init_util.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <basic/Tracer.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/string_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>

// C++ Headers
#include <cstdlib>
#include <string>
#include <cmath>

using basic::Error;
using basic::Warning;

static basic::Tracer TR( "apps.public.membrane.mp_transform" );

///////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {
	try {

		using namespace protocols::jd2;
		using namespace basic::options;
		using namespace protocols::membrane;
		using namespace core;
		using namespace core::pose;
		using namespace protocols::moves;
		using namespace protocols::membrane::visualize;
		using namespace protocols::membrane::geometry;

		devel::init(argc, argv);

		// read option to optimize membrane position
		bool optimize( false );
		if ( option[ OptionKeys::mp::transform::optimize_embedding ].user() ) {
			optimize = option[ OptionKeys::mp::transform::optimize_embedding ]();
		}

		// warning for using PDB lists
		if ( option[ OptionKeys::in::file::l ].user() ) {
			TR << "WARNING: When using lists, only provide PDBs of the same protein (for instance from multiple models) which have an identical span file! Only a single spanfile is allowed!" << std::endl;
		}

		// if optimizing membrane position
		if ( optimize == true ) {

			// create movers and add them in a SequenceMover
			AddMembraneMoverOP addmem( new AddMembraneMover() );
			OptimizeProteinEmbeddingMoverOP opt( new OptimizeProteinEmbeddingMover() );
			VisualizeEmbeddingMoverOP vis( new VisualizeEmbeddingMover() );
			SequenceMoverOP seq( new SequenceMover( addmem, opt, vis ) );

			// call SequenceMover in JD2
			JobDistributor::get_instance()->go( seq );

		} else {

			// create movers and add them in a SequenceMover
			AddMembraneMoverOP addmem( new AddMembraneMover() );
			TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover() );
			VisualizeEmbeddingMoverOP vis( new VisualizeEmbeddingMover() );
			SequenceMoverOP seq( new SequenceMover( addmem, transform, vis ) );

			// call SequenceMover in JD2
			JobDistributor::get_instance()->go( seq );
		}
	}
catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}
}
