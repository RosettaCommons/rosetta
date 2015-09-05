// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  apps/pilot/membrane/mp_transform.cc
///
/// @brief  RosettaMP Transform protein into a membrane
///
/// @author  Julia Koehler Leman (julia.koehler1982@gmail.com)
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @note   Last Updated: 5/18/15

// App headers
#include <devel/init.hh>

// Project Headers
#include <core/conformation/membrane/MembraneInfo.hh>


#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/OptimizeMembranePositionMover.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/membrane/visualize/VisualizeEmbeddingMover.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>

// Package Headers
#include <apps/benchmark/performance/init_util.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
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

/// @brief Create output filename
std::string output_filename( std::string infilename ) {

	using namespace basic::options;
	using namespace utility;

	// create output filename to /my/path/myprotein.pdb to myprotein
	const std::string tmp( file_basename( infilename ) );
	std::string output("");

	// get pathname
	if ( option[ OptionKeys::out::path::pdb ].user() ) {
		output = option[ OptionKeys::out::path::pdb ]().path() + utility::replace_in( infilename, ".pdb", "" );
	} else {
		output = utility::replace_in( infilename, ".pdb", "" );
	}

	// add counter to output filename
	output += "_0001.pdb";

	return output;

} // output filename

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

		// cry if PDB not given
		if ( ! option[OptionKeys::in::file::s].user() ) {
			throw new utility::excn::EXCN_Msg_Exception("Please provide PDB file!");
		}

		// read in pose
		Pose pose;
		std::string infilename( option[OptionKeys::in::file::s](1) );
		core::import_pose::pose_from_pdb( pose, infilename );

		// read option to optimize membrane position
		bool optimize( false );
		if ( option[ OptionKeys::mp::transform::optimize_embedding ].user() ) {
			optimize = option[ OptionKeys::mp::transform::optimize_embedding ]();
		}

		// if optimizing membrane position
		if ( optimize == true ) {

			// add membrane
			AddMembraneMoverOP addmem( new AddMembraneMover() );
			addmem->apply( pose );

			// transforms into fixed membrane and optimizes MEM (=dangling)
			OptimizeMembranePositionMoverOP opt( new OptimizeMembranePositionMover() );
			opt->apply( pose );

			// get the optimized embedding from flexible MEM
			core::Vector center( pose.conformation().membrane_info()->membrane_center() );
			core::Vector normal( pose.conformation().membrane_info()->membrane_normal() );
			EmbeddingDefOP current_emb( new EmbeddingDef( center, normal ) );

			// use current embedding and transform into default membrane (fixed)
			TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover( current_emb ) );
			transform->use_default_membrane( true );
			transform->apply( pose );

			// reset the membrane residue because it was changed during optimization
			SetMembranePositionMoverOP set( new SetMembranePositionMover() );
			set->apply( pose );

			// Visualize embedding
			VisualizeEmbeddingMoverOP vis( new VisualizeEmbeddingMover() );
			vis->apply( pose );

			// dump PDB
			pose.dump_pdb( output_filename( infilename ) );

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
