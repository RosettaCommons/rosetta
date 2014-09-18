// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:f;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file   FlexPepDock.cc
//
/// @brief Main file for applying flexible peptide docking protocol
/// @author Barak Raveh
/// @date August 05, 2008

//#define GL_GRAPHICS

// Mini-Rosetta headers

#include <numeric/random/random.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
//#include <core/chemical/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/threadsc.OptionKeys.gen.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
//#include <core/util/basic.hh>
#include <basic/Tracer.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using basic::T;
using basic::Error;
using basic::Warning;
using core::pose::Pose;



static thread_local basic::Tracer TR( "pilot_app.phoshporylate_position" );




///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

  using namespace core;
  using namespace basic::options;
  using namespace std;

  devel::init(argc, argv);

  pose::Pose pose;

	// verify inputs
  if ( !option[ OptionKeys::out::file::o ].user() ||
    !option[ OptionKeys::run::chain ].user() ||
	 !option[ OptionKeys::threadsc::nres ].user()  // # TODO: I "ride" here over the old nres, better find a param of our own...
  ) {
    TR  << "Usage: " << endl <<
      argv[0] << endl <<
      "  -s <fname> -chain <chain> -threadsc::nres <resid>" << endl <<
			"  -o <outfile> -database <minirosetta_db>"  << endl;
    exit(-1);
  }

  // read params and poses
	std::string start_file = option[ OptionKeys::in::file::s ][0];
	core::import_pose::pose_from_pdb( pose, start_file );
  string chain = option[ OptionKeys::run::chain ];
  Size pdb_res = option[ OptionKeys::threadsc::nres ];
  string output_fname = option[ OptionKeys::out::file::o ];

	// phosphorylate
  core::pose::PDBInfoCOP pdbinfo = pose.pdb_info();
  Size pose_res = pdbinfo->pdb2pose(chain[0], pdb_res);
	core::pose::add_variant_type_to_pose_residue( pose , chemical::PHOSPHORYLATION, pose_res );
  core::io::pdb::dump_pdb(pose, output_fname);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
