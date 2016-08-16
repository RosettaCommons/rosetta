// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   FlexPepDock.cc
//
/// @brief compare structures by various means, using FlexPepDock params file
//  @author Barak Raveh
/// @date August 05, 2008

//#define GL_GRAPHICS

// Mini-Rosetta headers
#include <devel/FlexPepDocking/FlexPepDockingProtocol.hh>

#include <numeric/random/random.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/util.hh>//option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/threadsc.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>

//#include <protocols/jobdist/standard_mains.hh>
//#include <protocols/moves/Mover.hh>

#include <ObjexxFCL/format.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>

using basic::T;
using basic::Error;
using basic::Warning;
using core::pose::Pose;


static THREAD_LOCAL basic::Tracer TR( "thread_sidechains" );


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

  using namespace core;
  using namespace basic::options;
  using namespace std;

  devel::init(argc, argv);

  pose::Pose src_pose,trg_pose;
  Size sfirst_res, tfirst_res, nres;
  // verify input params
  if ( !option[ OptionKeys::in::file::native ].user() ||
    !option[ OptionKeys::out::file::o ].user() ||
    !option[ OptionKeys::threadsc::src_chain ].user() ||
    !option[ OptionKeys::threadsc::trg_chain ].user() ||
    !option[ OptionKeys::threadsc::src_first_resid ].user() ||
    !option[ OptionKeys::threadsc::trg_first_resid ].user() ||
    !option[ OptionKeys::threadsc::nres ].user()
  ) {
    TR  << "Usage: " << endl <<
      argv[0] << endl <<
      "  -in:file:native:<fname> -s <fname> -out:file:<output_file> " <<
      "  -threadsc:src_chain:<chain> -threadsc:src_first_resid:<resid> " << endl <<
      "  -threadsc:trg_chain:<chain> -threadsc:trg_first_resid:<resid>" << endl <<
      "  -database <minirosetta_db>"  << endl;
    exit(-1);
  }
  // read poses
  core::import_pose::pose_from_file( trg_pose, basic::options::start_file() , core::import_pose::PDB_file);
  string native_fname = option[ OptionKeys::in::file::native ];
  core::import_pose::pose_from_file( src_pose, native_fname, core::import_pose::PDB_file);
  string output_fname = option[ OptionKeys::out::file::o ];
  // compute threading range
  string schain_pdb = option[ OptionKeys::threadsc::src_chain ];
  string tchain_pdb = option[ OptionKeys::threadsc::trg_chain ];
  Size sresi_pdb =  option[ OptionKeys::threadsc::src_first_resid ];
  Size tresi_pdb =  option[ OptionKeys::threadsc::trg_first_resid ];
  sfirst_res = src_pose.pdb_info()->pdb2pose(schain_pdb[0],sresi_pdb);
  tfirst_res = trg_pose.pdb_info()->pdb2pose(tchain_pdb[0],tresi_pdb);
  nres = option[ OptionKeys::threadsc::nres ];
  // TODO: also verify no overflow within the same chain
  if(sfirst_res + nres - 1 > src_pose.total_residue() ||
    tfirst_res + nres - 1 > src_pose.total_residue()) {
    TR << "ERROR: specified resiude IDs beyond range" << std::endl;
    exit(-1);
  }

  // thread residues
  TR << "Threading residues from [" << native_fname << "] to [" << basic::options::start_file() << "]" << endl;
  for(Size i=0; i < nres; i++) {
    Size sresi = sfirst_res + i;
    Size tresi = tfirst_res + i;
    if(!src_pose.residue(sresi).is_protein() ||
      !trg_pose.residue(tresi).is_protein()) {
      TR << "ERROR: non-protein residue encountered" << std::endl;
      exit(-1);
    }
    conformation::Residue src_res = src_pose.residue(sresi);
    bool orient_bb = true; // orient new residue to existing backbone
    trg_pose.replace_residue(tresi,src_res,orient_bb);
  }

  // output overlayed pose
  TR << "Output to [" << output_fname << "]" << endl;
  core::io::pdb::dump_pdb(trg_pose, output_fname);

  	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
