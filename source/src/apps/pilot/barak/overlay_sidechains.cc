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
/// @brief Main file for applying flexible peptide docking protocol
/// @author Barak Raveh
/// @date August 05, 2008

//#define GL_GRAPHICS

// Mini-Rosetta headers
#include <devel/FlexPepDocking/FlexPepDockingProtocol.hh>

#include <numeric/random/random.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/util.hh>//option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/pose/Pose.hh>
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


static THREAD_LOCAL basic::Tracer TR( "pilot_app.barak.overlay_sidechains" );

bool verify_identical(Pose& pose1, Pose& pose2)
{
  if(pose1.total_residue() != pose2.total_residue()) {
    return false;
  }
  for(core::Size resi=1; resi < pose1.total_residue(); resi++) {
    bool p1_protein = pose1.residue(resi).is_protein();
    bool p2_protein = pose2.residue(resi).is_protein();
    if(p1_protein != p2_protein) {
      return false;
    }
    if(!p1_protein){ // chis are relevant only for protein residues
      continue;
    }
    if(pose1.residue_type(resi).nchi() != pose2.residue_type(resi).nchi()) {
      return false;
    }
  }
  return true;
}


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

  // read params and poses
  core::import_pose::pose_from_file( pose, basic::options::start_file() , core::import_pose::PDB_file);
  if ( !option[ OptionKeys::in::file::native ].user() ) {
    TR << "specify reference native file (-native option)" << std::endl;
    exit(-1);
  }
  string native_fname = option[ OptionKeys::in::file::native ];
  pose::Pose ref_pose;
  core::import_pose::pose_from_file( ref_pose, native_fname, core::import_pose::PDB_file);
  if ( !option[ OptionKeys::out::file::o ].user() ) {
    TR << "specify output file (-o option)" << std::endl;
    exit(-1);
  }
  string output_fname = option[ OptionKeys::out::file::o ];

  // verify same number of chi angles in pose and ref_pose
  if(!verify_identical(pose, ref_pose)){
    TR << "mismatching # of residues" << std::endl;
    exit(-1);
  }

  // overlay chi values
  TR << "Impose native chi angles of [" << native_fname << "] to start structure:" << endl;
  for(core::Size resi=1; resi < pose.total_residue(); resi++)
    {
      if(!pose.residue(resi).is_protein()){
	continue;
      }
      for(core::Size j_chi=1 ; j_chi <= pose.residue_type(resi).nchi();
	  j_chi++) {
	Real ref_chi = ref_pose.chi(j_chi,resi);
	pose.set_chi(j_chi, resi, ref_chi);
      }
    }

  // output overlayed pose
  TR << "Output to [" << output_fname << "]" << endl;
  core::io::pdb::traced_dump_pdb(TR, pose, output_fname);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
