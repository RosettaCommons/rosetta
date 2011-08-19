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
#include <core/chemical/util.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/options/util.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/run.OptionKeys.gen.hh>
#include <core/options/keys/threadsc.OptionKeys.gen.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/util/basic.hh>
#include <core/util/Tracer.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

using core::util::T;
using core::util::Error;
using core::util::Warning;
using core::pose::Pose;


static numeric::random::RandomGenerator RG(12321); // <- Magic number, do not change it!!!

static core::util::Tracer TR("pilot_app.phoshporylate_position");




///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
  using namespace core;
  using namespace options;
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
  io::pdb::pose_from_pdb( pose, options::start_file() );
  string chain = option[ OptionKeys::run::chain ];
  Size pdb_res = option[ OptionKeys::threadsc::nres ];
  string output_fname = option[ OptionKeys::out::file::o ];

	// phosphorylate
  core::pose::PDBInfoCOP pdbinfo = pose.pdb_info();
  Size pose_res = pdbinfo->pdb2pose(chain[0], pdb_res);
	chemical::add_variant_type_to_pose_residue( pose , chemical::PHOSPHORYLATION, pose_res );
  core::io::pdb::dump_pdb(pose, output_fname);

  exit(0);

}
