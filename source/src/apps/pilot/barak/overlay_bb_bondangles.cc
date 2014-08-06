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
#include <devel/FlexPepDocking/FlexPepDockingProtocol.hh>

#include <numeric/random/random.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/util.hh>//option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/pose/PDBInfo.hh>
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


static numeric::random::RandomGenerator RG(12321); // <- Magic number, do not change it!!!

static basic::Tracer TR("pilot_app.barak.overlay_sidechains");


void overlay_bb_bondangles(
  Pose& ref_pose,
  Pose& pose,
  core::Size ref_fromres,
  core::Size fromres,
  core::Size nres)
{
  using namespace core;
  using namespace core::conformation;
  using namespace id;

  core::conformation::Conformation& conformation
    = pose.conformation();
  core::conformation::Conformation const& ref_conformation
    = ref_pose.conformation();
  for(Size k=0; k < nres ; k++)
    {
      Size resid = fromres + k;
      Size ref_resid = ref_fromres + k;
      Residue const & rsd( conformation.residue( resid ) );
      Residue const & ref_rsd( ref_conformation.residue( ref_resid ) );

      // bond-angle just before residue (C=N-CA if protein)
      if ( ref_resid>1 && ref_rsd.is_polymer() && !ref_rsd.is_lower_terminus() )
	{
	  Residue const& prev_ref_rsd( ref_conformation.residue( ref_resid-1 ) );
	  Size const nbb_prev( prev_ref_rsd.n_mainchain_atoms() );
	  AtomID
	    bb1_ref    ( prev_ref_rsd.mainchain_atom( nbb_prev ),  ref_resid-1 ),
	    bb2_ref    ( ref_rsd.mainchain_atom( 1 ),  ref_resid ),
	    bb3_ref    ( ref_rsd.mainchain_atom( 2 ),  ref_resid );
	  Residue const& prev_rsd( conformation.residue( resid-1 ) );
	  AtomID
	    bb1    ( prev_rsd.mainchain_atom( nbb_prev ),  resid-1 ),
	    bb2    ( rsd.mainchain_atom( 1 ),  resid ),
	    bb3    ( rsd.mainchain_atom( 2 ),  resid );
	  double ref_angle = ref_conformation.bond_angle( bb1_ref, bb2_ref, bb3_ref) ;
	  conformation.set_bond_angle(bb1, bb2, bb3, ref_angle);
	}

      // intra-residue mainchain bonds and angles
      Size const nbb( ref_rsd.n_mainchain_atoms() );
      assert( nbb >= 2 ); // or logic gets a bit trickier
      for ( Size i=2; i < nbb; ++i ) {
	AtomID
	  bb1_ref    ( ref_rsd.mainchain_atom(i-1),  ref_resid ),
	  bb2_ref    ( ref_rsd.mainchain_atom(  i),  ref_resid ),
	  bb3_ref    ( ref_rsd.mainchain_atom(i+1),  ref_resid );
	AtomID
	  bb1    ( rsd.mainchain_atom(i-1),  resid ),
	  bb2    ( rsd.mainchain_atom(  i),  resid ),
	  bb3    ( rsd.mainchain_atom(i+1),  resid );
	double ref_angle = ref_conformation.bond_angle( bb1_ref, bb2_ref, bb3_ref) ;
	conformation.set_bond_angle(bb1, bb2, bb3, ref_angle);
      }

      // bond angle between residues (CA-C=N if protein)
      if ( ref_resid < ref_conformation.size() && ref_rsd.is_polymer() && !ref_rsd.is_upper_terminus() )
	{
	  Residue const & next_ref_rsd(  ref_conformation.residue( ref_resid+1 ) );
	  AtomID
	    bb1_ref    ( ref_rsd.mainchain_atom( nbb-1 ),  ref_resid ),
	    bb2_ref    ( ref_rsd.mainchain_atom( nbb   ),  ref_resid ),
	    bb3_ref    ( next_ref_rsd.mainchain_atom( 1 ),  ref_resid+1 );
	  Residue const& next_rsd( conformation.residue( resid+1 ) );
	  AtomID
	    bb1    ( rsd.mainchain_atom( nbb-1 ),  resid ),
	    bb2    ( rsd.mainchain_atom( nbb   ),  resid ),
	    bb3    ( next_rsd.mainchain_atom( 1 ),  resid+1 );
	  double ref_angle = ref_conformation.bond_angle( bb1_ref, bb2_ref, bb3_ref) ;
	  conformation.set_bond_angle(bb1, bb2, bb3, ref_angle);
	}
    }
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
  core::import_pose::pose_from_pdb( pose, basic::options::start_file() );
  if ( !option[ OptionKeys::in::file::native ].user() ) {
    TR << "specify reference native file (-native option)" << std::endl;
    exit(-1);
  }
  string native_fname = option[ OptionKeys::in::file::native ];
  pose::Pose ref_pose;
  core::import_pose::pose_from_pdb( ref_pose, native_fname);
  if ( !option[ OptionKeys::out::file::o ].user() ) {
    TR << "specify output file (-o option)" << std::endl;
    exit(-1);
  }
  string output_fname = option[ OptionKeys::out::file::o ];

  // overlay bond angles
  core::pose::PDBInfoCOP ref_pdbinfo = ref_pose.pdb_info();
  core::pose::PDBInfoCOP pdbinfo = pose.pdb_info();
  Size ref_fromres = ref_pdbinfo->pdb2pose('C',0);
  Size trg_fromres = pdbinfo->pdb2pose('C',1);
  Size nres = 10;
  overlay_bb_bondangles(ref_pose, pose, ref_fromres, trg_fromres, nres);

  // output overlayed pose
  TR << "Output to [" << output_fname << "]" << endl;
  core::io::pdb::traced_dump_pdb(TR, pose, output_fname);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
