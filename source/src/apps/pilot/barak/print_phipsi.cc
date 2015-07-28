// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
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

#include <numeric/constants.hh>
#include <numeric/random/random.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/util.hh>//option.hh>
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


// convert radians to degrees
double rad2deg(double rad)
{
  return rad * 360 / numeric::constants::d::pi_2;
}


// enumerate and print residue inter and intra bond angles
void
printBondAngles(core::Size const seqpos, core::pose::Pose & pose)
{
  using namespace core;
  using namespace core::conformation;
  using namespace id;

  core::conformation::Conformation const& conformation = pose.conformation();
  char ch =  pose.pdb_info()->chain(seqpos);
  Size resid = pose.pdb_info()->number(seqpos);
  Residue const & rsd( conformation.residue( seqpos ) );

  // bond-angle just before residue (C=N-CA if protein)
  if ( seqpos>1 && rsd.is_polymer() && !rsd.is_lower_terminus() ) {
    Residue const& prev_rsd( conformation.residue( seqpos-1 ) );
    Size const nbb_prev( prev_rsd.n_mainchain_atoms() );
    AtomID
      bb1    ( prev_rsd.mainchain_atom( nbb_prev ),  seqpos-1 ),
      bb2    ( rsd.mainchain_atom( 1 ),  seqpos ),
      bb3    ( rsd.mainchain_atom( 2 ),  seqpos );
    std::string
      bb1_name ( prev_rsd.atom_name(prev_rsd.mainchain_atom(nbb_prev)) ),
      bb2_name ( rsd.atom_name(rsd.mainchain_atom( 1 )) ),
      bb3_name ( rsd.atom_name(rsd.mainchain_atom( 2 )) );
    Size prev_resid = pose.pdb_info()->number(seqpos-1);
    double angle = rad2deg( conformation.bond_angle( bb1, bb2, bb3) );
    std::cout << "Chain " << ch << " bond-angle ["
	      << bb1_name << prev_resid << ", "
	      << bb2_name << resid << ", "
	      << bb3_name << resid
	      << "]  =  " << angle << " Deg" << std::endl;

  }

  // intra-residue mainchain bonds and angles
  Size const nbb( rsd.n_mainchain_atoms() );
  assert( nbb >= 2 ); // or logic gets a bit trickier
  for ( Size i=2; i<= nbb; ++i ) {
    AtomID
      bb1    ( rsd.mainchain_atom(i-1),  seqpos ),
      bb2    ( rsd.mainchain_atom(  i),  seqpos );
    std::string
      bb1_name ( rsd.atom_name(rsd.mainchain_atom(i-1)) ),
      bb2_name ( rsd.atom_name(rsd.mainchain_atom(i  )) );
    double length = conformation.bond_length( bb1, bb2);
    std::cout << "Chain " << ch << " bond-length ["
	      << bb1_name << resid << ", "
	      << bb2_name << resid
	      << "]  =  " << length << " A" << std::endl;
    if ( i<nbb )
      {
	AtomID  bb3    ( rsd.mainchain_atom(i+1),  seqpos );
	std::string bb3_name ( rsd.atom_name(rsd.mainchain_atom(i+1)) );
	double angle = rad2deg( conformation.bond_angle( bb1, bb2, bb3 ) );
    std::cout << "Chain " << ch << " bond-angle ["
	      << bb1_name << resid << ", "
	      << bb2_name << resid << ", "
	      << bb3_name << resid
	      << "]  =  " << angle << " Deg" << std::endl;
      }
  }

  // bond angle between residues (CA-C=N if protein)
  if ( seqpos < conformation.size() && rsd.is_polymer() && !rsd.is_upper_terminus() ) {
    Residue const & next_rsd(  conformation.residue( seqpos+1 ) );
    AtomID
      bb1    ( rsd.mainchain_atom( nbb-1 ),  seqpos ),
      bb2    ( rsd.mainchain_atom( nbb   ),  seqpos ),
      bb3    ( next_rsd.mainchain_atom( 1 ),  seqpos+1 );
    std::string
      bb1_name ( rsd.atom_name(rsd.mainchain_atom(nbb-1)) ),
      bb2_name ( rsd.atom_name(rsd.mainchain_atom(nbb  )) ),
      bb3_name ( next_rsd.atom_name(next_rsd.mainchain_atom(1)) );
    Size next_resid = pose.pdb_info()->number(seqpos+1);
    double length = conformation.bond_length(bb2, bb3);
    std::cout << "Chain " << ch << " bond-length ["
	      << bb2_name << resid      << ", "
	      << bb3_name << next_resid
	      << "]  =  " << length << " A" << std::endl;
    double angle = rad2deg( conformation.bond_angle( bb1, bb2, bb3) );
    std::cout << "Chain " << ch << " bond-angle ["
	      << bb1_name << resid      << ", "
	      << bb2_name << resid      << ", "
	      << bb3_name << next_resid
	      << "]  =  " << angle << " Deg" << std::endl;
  }
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
  using namespace core;
  using namespace std;

  devel::init(argc, argv);

  pose::Pose pose;
  core::import_pose::pose_from_pdb( pose, basic::options::start_file() );
  Size const nres( pose.total_residue() );

  std::cout << "chain tor_type res_id tor_angle" << std::endl;
  for(Size i = 1; i <= nres; ++i)
    {
      std::cout << setw(5) << pose.pdb_info()->chain(i)
		<<  " phi      " << setw(6) << pose.pdb_info()->number(i)
		<<  " " << pose.phi(i) << std::endl;
      std::cout << setw(5) << pose.pdb_info()->chain(i)
		<<  " psi      " << setw(6) << pose.pdb_info()->number(i)
		<<  " " << pose.psi(i) << std::endl;
      std::cout << setw(5) << pose.pdb_info()->chain(i)
		<<  " omega      " << setw(6) << pose.pdb_info()->number(i)
		<<  " " << pose.omega(i) << std::endl;
      printBondAngles(i, pose);
    }

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
