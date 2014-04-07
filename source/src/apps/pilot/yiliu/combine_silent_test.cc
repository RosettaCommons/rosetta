// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author James Thompson

// libRosetta headers


#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <protocols/moves/NullMover.hh>
#include <protocols/jobdist/standard_mains.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>

#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

using basic::T;
using basic::Warning;
using basic::Error;

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>



using namespace core;
using namespace core::chemical;
using utility::file::FileName;
using utility::vector1;
static basic::Tracer TR("apps.pilot.yiliu.silent");

int
main( int argc, char* argv [] ) {
	try {
	// options, random initialization
	devel::init( argc, argv );

	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;

  //YL, declare silent file data
	SilentFileData sfd, sfd_out;

  //YL, declare residue type
	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
		option[ in::file::residue_type_set ]()
	);

	core::pose::Pose native_pose, current_pose;
	if ( option[ in::file::native ].user() ) {
		TR << option[ in::file::native].user() << std::endl;
		core::import_pose::pose_from_pdb( native_pose, *rsd_set, option[ in::file::native ]() );
	}

	//YL, declare score function
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
	//YL, declare NullMover
	protocols::moves::NullMover mover;
	protocols::jobdist::not_universal_main( mover );

//	MetaPoseInputStream input;
//	if ( option[ in::file::s ].user() ) {
//		PoseInputStreamOP pdb_input( new PDBPoseInputStream( option[ in::file::s ]() ) );
//		input.add_pose_input_stream( pdb_input );
//	}
//
//	if ( option[ in::file::silent ].user() ) {
//		utility::vector1< std::string > tags;
//		PoseInputStreamOP silent_input;
//		if ( option[ in::file::tags ].user() ) {
//			tags = option[ in::file::tags ]();
//			silent_input = new SilentFilePoseInputStream( option[ in::file::silent ](), option[ in::file::tags ]() );
//		} else {
//			silent_input = new SilentFilePoseInputStream( option[ in::file::silent ]() );
//		}
//		input.add_pose_input_stream( silent_input );
//	}
//	// if ( option[ in::file::l ].user() )  pdb_input->add_file_list( option[ in::file::l ]() );
//
//
//	while( input.has_another_pose() ) {
//		input.fill_pose( current_pose, *rsd_set );
//		(*scorefxn)( current_pose );
//
//		SilentStructOP ss = SilentStructFactory::get_instance()->get_silent_struct();
//		ss->fill_struct( current_pose );
//
//		// add rmsd information if a native was given.
//		if ( option[ in::file::native ].user() ) {
//			// calculate RMS
//			Real rmsd = core::scoring::CA_rmsd( native_pose, current_pose );
//			ss->add_energy( "CA_rmsd", rmsd );
//		}
//		std::cerr << "writing " << ss->decoy_tag() << std::endl;
//		sfd_out.write_silent_struct( *ss, option[ out::file::silent ]() );
//	}

	return 0;
	 } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
} // main
