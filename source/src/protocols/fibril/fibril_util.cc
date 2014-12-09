// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/fibril/fibril_util.hh
/// @brief utility functions for handling of fibril symmetry modeling
/// @author Lin Jiang

// Unit headers
#include <protocols/fibril/fibril_util.hh>
#include <core/pose/symmetry/util.hh>
// AUTO-REMOVED #include <core/conformation/symmetry/util.hh>

// AUTO-REMOVED #include <core/conformation/symmetry/SymDof.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
// AUTO-REMOVED #include <core/scoring/symmetry/SymmetricEnergies.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/kinematics/FoldTree.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// AUTO-REMOVED #include <core/init/init.hh>

// Utility functions
#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <numeric/random/random.hh>

// Package Headers
#include <core/kinematics/Edge.hh>

#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>

#include <core/conformation/symmetry/SymmData.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>
//namespaces
using namespace core;
using namespace core::conformation;
using namespace core::conformation::symmetry;
using namespace utility;


namespace protocols {
namespace fibril {

static thread_local basic::Tracer TR( "protocols.moves.symmetry.fibril_util" );

//lin functions for xyz coordinate
bool
set_xyz_coord(
   Vector const & origin,
   Vector const & x,
   Vector const & y,
   Vector & z,
   bool check_angle
)
{
  if( check_angle ) {
    Real const angle_xy ( numeric::angle_radians(x,origin,y) );
    if( angle_xy - 0.5*numeric::NumericTraits<Real>::pi() > 1e-3 ) {
      return false;
    }
  }
  z =  origin + ( x - origin ).cross( y - origin ) ;
  return true;
}


/////////////////////////////////////////////////////////////////////////////
void
reorient_extended_fibril(
  Conformation & src_conformation,
  SymmData & symmdata
)
{
  using namespace core::kinematics;
  using namespace id;

  kinematics::FoldTree f( src_conformation.fold_tree() );
  Size anchor ( utility::string2int(symmdata.get_anchor_residue()) );
  Size next_anchor ( anchor <= src_conformation.size() - 2 ? anchor + 2 : anchor - 2 );
  chemical::ResidueType const& rt1 ( src_conformation.residue_type ( anchor ) );
  chemical::ResidueType const& rt2 ( src_conformation.residue_type ( next_anchor ) );
  AtomID a1( rt1.atom_index ("C") , anchor );
  AtomID a2( rt2.atom_index ("C") , next_anchor );
  AtomID a3( rt1.atom_index ("O") , anchor );
  Vector origin ( src_conformation.xyz(a1) );
  Vector y ( src_conformation.xyz(a2) );
  Vector z ( src_conformation.xyz(a3) );
  Vector x;
  set_xyz_coord( origin, y, z, x, false);
  x = ( x - origin ).normalized() + origin;
  z = ( z - origin ).normalized() + origin;
  Stub const rot_stub ( origin, x, z );
  Stub const src_stub ( Vector(0,0,0), Vector(1,0,0), Vector(0,0,1) );

  for ( Size i = 1; i <= src_conformation.size(); ++i ) {
    for ( Size j = 1; j <= src_conformation.residue_type(i).natoms(); ++j ) {
      AtomID id( j, i );
      Vector const old_xyz( src_conformation.xyz(id) );
      Vector const new_xyz( src_stub.local2global( rot_stub.global2local( old_xyz ) ) );
      src_conformation.set_xyz( id, new_xyz );
    }
  }
  //pose.apply_transform_Rx_plus_v( rt.get_rotation(), rt.get_translation() );
}

/////////////////////////////////////////////////////////////////////////////
void
make_symmetric_fibril(
  pose::Pose & pose
)
{
  using namespace basic::options;
  using namespace basic::options::OptionKeys;


  SymmData symmdata( pose.n_residue(), pose.num_jump() );
  std::string symm_def = option[ OptionKeys::symmetry::symmetry_definition ];
  symmdata.read_symmetry_data_from_file(symm_def);
  if( option[ in::file::native ].user() ) {
    pose::Pose monomer_pose;
		core::import_pose::pose_from_pdb( monomer_pose, option[ in::file::native ]().name() );
    protocols::loops::Loops loop1, loop2;
    bool default_loop_file_is_present = option[ OptionKeys::loops::loop_file ].user();
    loop1 = protocols::loops::Loops( default_loop_file_is_present );
		if(  !default_loop_file_is_present ) {
			//make for monomer_pose
			loop1.push_back( protocols::loops::Loop( 1, monomer_pose.total_residue(), 0,  0.0, false) );
		}
		if( option[ OptionKeys::loops::extended_loop_file ].user() ) {
			loop2 = protocols::loops::Loops( option[ OptionKeys::loops::extended_loop_file ]() );
		} else {
			loop2 = loop1;
		}
		//monomer_pose.dump_pdb("before.pdb");
    protocols::fibril::superimpose_pose_on_subset_bb( pose, monomer_pose, loop1, loop2 );
		//monomer_pose.dump_pdb("after.pdb");
  } else if( option[ in::file::alignment ].user() ) {
    protocols::fibril::reorient_extended_fibril( pose.conformation(), symmdata );
  }
	core::pose::symmetry::make_symmetric_pose( pose, symmdata );
  std::cout<<"fold tree: "<<pose.fold_tree()<<std::endl;
}

/////////////////////////////////////////////////////////////////////////////
void
superimpose_pose_on_subset_bb(
  pose::Pose& pose,
  pose::Pose& ref_pose,
  protocols::loops::Loops core,
  protocols::loops::Loops ref_core
)
{

  using namespace protocols;
  using namespace core::scoring;

  //set atom map for superimpose_pose
  core::id::AtomID_Map< id::AtomID > atom_map;
  core::pose::initialize_atomid_map( atom_map, pose, core::id::BOGUS_ATOM_ID ); // maps every atomid to bogus atom

  utility::vector1< core::id::AtomID > ref_ids;
  utility::vector1< core::id::AtomID > ids;

  for ( loops::Loops::const_iterator region = core.begin(); region != core.end(); ++region ) {
    for ( core::Size i=region->start(); i<= region->stop(); ++i ){
      core::id::AtomID dummy_atomid1( pose.residue(i).atom_index("CA"), i);
      ids.push_back(dummy_atomid1);
      core::id::AtomID dummy_atomid2( pose.residue(i).atom_index("N"), i);
      ids.push_back(dummy_atomid2);
      core::id::AtomID dummy_atomid3( pose.residue(i).atom_index("C"), i);
      ids.push_back(dummy_atomid3);
      core::id::AtomID dummy_atomid4( pose.residue(i).atom_index("O"), i);
      ids.push_back(dummy_atomid4);
    }
  }

  for ( loops::Loops::const_iterator region = ref_core.begin(); region != ref_core.end(); ++region ) {
    for ( core::Size i=region->start(); i<= region->stop(); ++i ){
      core::id::AtomID dummy_atomid1( ref_pose.residue(i).atom_index("CA"), i);
      ref_ids.push_back(dummy_atomid1);
      core::id::AtomID dummy_atomid2( ref_pose.residue(i).atom_index("N"), i);
      ref_ids.push_back(dummy_atomid2);
      core::id::AtomID dummy_atomid3( ref_pose.residue(i).atom_index("C"), i);
      ref_ids.push_back(dummy_atomid3);
      core::id::AtomID dummy_atomid4( ref_pose.residue(i).atom_index("O"), i);
      ref_ids.push_back(dummy_atomid4);
    }
  }
  assert( ids.size()== ref_ids.size());

  for ( core::Size i=1; i<=ids.size(); ++i ) {
    atom_map.set( ids[i], ref_ids[i] );
  }
  superimpose_pose( pose, ref_pose, atom_map );
}

} // fibril
} // protocols
