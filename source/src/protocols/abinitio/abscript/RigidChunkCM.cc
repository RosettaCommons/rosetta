// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/abscript/RigidChunkCMCM.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/RigidChunkCM.hh>
#include <protocols/abinitio/abscript/RigidChunkCMCreator.hh>

// Package headers
#include <core/environment/DofPassport.hh>

#include <protocols/environment/DofUnlock.hh>
#include <protocols/environment/claims/CutBiasClaim.hh>
#include <protocols/environment/claims/JumpClaim.hh>
#include <protocols/environment/claims/XYZClaim.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/AtomTree.hh>

// Project headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.tmpl.hh>
#include <protocols/loops/LoopsFileIO.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

//Utility Headers
#include <utility/tag/Tag.hh>

#include <numeric/random/random.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/xyz.functions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

//Req'd on WN32
#include <basic/datacache/WriteableCacheableMap.hh>

static numeric::random::RandomGenerator RG(1892476846);
static basic::Tracer tr("protocols.abinitio.abscript.RigidChunkCM", basic::t_info);

namespace protocols {
namespace abinitio {
namespace abscript {

using namespace core::environment;
using namespace protocols::environment;

// creator
std::string
RigidChunkCMCreator::keyname() const {
  return RigidChunkCMCreator::mover_name();
}

protocols::moves::MoverOP
RigidChunkCMCreator::create_mover() const {
  return new RigidChunkCM;
}

std::string
RigidChunkCMCreator::mover_name() {
  return "RigidChunkCM";
}

RigidChunkCM::RigidChunkCM():
  Parent()
{}

RigidChunkCM::RigidChunkCM( std::string const& label,
                            loops::Loops const& rigid_core,
                            core::pose::Pose const& template_pose ):
  Parent(),
  label_( label ),
  rigid_core_( rigid_core ),
  template_( new core::pose::Pose(template_pose) ) {}


loops::Loops read_rigid_core( std::string const& file){

  loops::PoseNumberedLoopFileReader reader;
  reader.hijack_loop_reading_code_set_loop_line_begin_token( "RIGID" );

  std::ifstream infile( file.c_str() );

  if (!infile.good()) {
    tr.Error << "[ERROR] Error opening RBSeg file '" << file << "'" << std::endl;
    throw utility::excn::EXCN_RosettaScriptsOption( "[ERROR] Error opening RBSeg file '" + file + "'" );
  }

  return loops::Loops( reader.read_pose_numbered_loops_file( infile, file, false ) );
}

void RigidChunkCM::parse_my_tag( utility::tag::TagCOP tag,
                                 basic::datacache::DataMap&,
                                 protocols::filters::Filters_map const&,
                                 protocols::moves::Movers_map const&,
                                 core::pose::Pose const& ){

  using namespace basic::options;
  label( tag->getOption< std::string >( "label", "BASE" ) );

  // Size min_loop_size_ = tag->getOption("min_loop_size", 0 );

  // bool bUseThreadingJobLoops_ = tag->getOption( "threading_loops", false );

  loops::Loops rigid_core;
  if( tag->hasOption("region_file") ){
    rigid_core = read_rigid_core( tag->getOption< std::string >( "region_file", "[NOT_SET]" ) );
  } else if ( tag->hasOption( "region" ) ){
    runtime_assert( false ); // unimplemented
  } else if ( tag->hasOption("loop_file") ){
    runtime_assert( false ); // unimplemented
  } else {
    std::string const err = "RigidChunkCM did not recieve any information about which chunks to fix. Use options 'region_file' to specify.";
    throw utility::excn::EXCN_BadInput( err );
  }

  if( tag->hasOption("template") ){
    std::string file = tag->getOption< std::string >( "template" );
    if( file == "INPUT" ){
      template_ = NULL;
    } else {
      core::pose::PoseOP p = new core::pose::Pose();
      core::import_pose::pose_from_pdb( *p,
                                        *core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ),
                                        file );
      template_ = p;
    }
  } else {
    throw utility::excn::EXCN_BadInput( "RigidChunkCM requires a template pdb with coordinates for rigid regions.");
  }

  rigid_core_ = select_parts( rigid_core, tag->getOption( "random_grow_by",
                                                          option[ OptionKeys::loops::random_grow_loops_by ]() ) );

}

claims::EnvClaims RigidChunkCM::yield_claims( core::pose::Pose const& in_p,
                                              basic::datacache::WriteableCacheableMapOP ){
  using namespace claims;
  EnvClaims claims;

  if( !template_ ){
    tr.Debug << "Building template from broker-time pose." << std::endl;
    template_ = new core::pose::Pose( in_p );
  }

  if( in_p.is_fullatom() && !template_->is_fullatom() ){
    throw utility::excn::EXCN_BadInput( "If the modelled pose is full atom, the template must also be full atom." );
  }

  loops::Loop loop_prev;

  for ( loops::Loops::const_iterator loop_it = rigid_core_.begin();
        loop_it != rigid_core_.end(); ++loop_it ) {
    XYZClaimOP xyz_claim = new XYZClaim( this,
                                         label(),
                                         std::make_pair( loop_it->start(),
                                                         loop_it->stop() ) );
    xyz_claim->strength( EXCLUSIVE, EXCLUSIVE );
    claims.push_back( xyz_claim );


    XYZClaimOP support_claim = new XYZClaim( this );
    if( loop_it->start() > 1 ){
      support_claim->add_position( LocalPosition( label(), loop_it->start()-1 ) );
    } if( loop_it->stop() < template_->total_residue() ){
      support_claim->add_position( LocalPosition( label(), loop_it->stop()+1 ) );
    }
    support_claim->strength( DOES_NOT_CONTROL, MUST_CONTROL );
    claims.push_back( support_claim );

    claims.push_back( new CutBiasClaim( this, label(), std::make_pair( loop_it->start(), loop_it->stop() ), 0.0 ) );

    if( loop_prev.start() != 0 && loop_prev.stop() != 0 ){
      core::Size jump_start = loop_prev.start() + ( ( loop_prev.stop() - loop_prev.start() ) / 2 );
      core::Size jump_end   = loop_it->start() +  ( ( loop_it->stop()  - loop_it->start()  ) / 2 );


      JumpClaimOP j_claim = new JumpClaim( this,
                                           "RigidChunkJump"+utility::to_string( claims.size()/4 ),
                                           LocalPosition( label(), jump_start ),
                                           LocalPosition( label(), jump_end ) );
      j_claim->strength( EXCLUSIVE, EXCLUSIVE );
      j_claim->physical( false );

      claims.push_back( j_claim );
    }
    loop_prev = *loop_it;
  }

  return claims;
}

bool missing_density( core::pose::Pose const& pose, loops::Loops const& core ){
  bool missing_density = false;

  //sanity check: no missing density in backbon in any of the rigid_core residues?
  for ( loops::Loops::const_iterator it = core.begin(); it!=core.end(); ++it ) {
    for ( core::Size pos = it->start(); pos <=it->stop(); ++pos ) {
      // Do we really have Sidechains ?
      // check this my making sure that no SC atom is more than 20A (?) away from CA
      numeric::xyzVector< core::Real> ca_pos = pose.residue( pos ).atom("CA").xyz();
      numeric::xyzVector< core::Real> n_pos = pose.residue( pos ).atom("N").xyz();
      numeric::xyzVector< core::Real> o_pos = pose.residue( pos ).atom("O").xyz();
      if ( ( n_pos - ca_pos).length() > 20 || ( ( n_pos - o_pos ).length() > 20 ) ) {
        tr.Error << "missing backbone in rigid-chunk at " << pos << std::endl;
        missing_density = true;
      }
    }
  }

  return missing_density;
}

/////////////// MAGIC FPD CODE FROM RIGID CHUNK CLAIMER /////////////////////////////
void fix_internal_coords_of_siblings( core::pose::Pose& pose,
                                        core::pose::Pose const& ref_pose,
                                        core::id::AtomID atom,
                                        core::id::AtomID ref_atom ) {
  using namespace core::id;

  runtime_assert( atom.rsd() >= 1 && atom.rsd() <= pose.total_residue() );
  runtime_assert( pose.conformation().atom_tree().has( atom ) );
  runtime_assert( ref_pose.conformation().atom_tree().has( ref_atom ) );

  bool has_par1( pose.conformation().atom_tree().atom( atom ).parent() );
  bool ref_has_par1( ref_pose.conformation().atom_tree().atom( ref_atom ).parent() );

  //folding direction matters for the angle we have to set...hence find the parent atoms and get the angle
  AtomID par1O;
  AtomID ref_par1O;
  if ( has_par1 && ref_has_par1 ) {
    par1O=pose.conformation().atom_tree().atom( atom ).parent()->id();
    std::string const & aname( pose.residue( par1O.rsd() ).atom_name( par1O.atomno() ));
    ref_par1O=core::id::AtomID( ref_pose.residue( par1O.rsd() ).atom_index( aname ), par1O.rsd() );
  }	else {
    tr.Warning << "cannot fix internal coords of " << atom << " in RigidChunk because 1st parent is missing " << std::endl;
    return;
  }
  bool has_par2( pose.conformation().atom_tree().atom( par1O ).parent() );
  bool ref_has_par2( ref_pose.conformation().atom_tree().atom( ref_par1O ).parent() );
  core::id::AtomID par2O;
  core::id::AtomID ref_par2O;
  if ( has_par2 && ref_has_par2 ) {
    par2O=pose.conformation().atom_tree().atom( par1O ).parent()->id();
    std::string const & aname( pose.residue(par2O.rsd()).atom_name( par2O.atomno() ) );
    ref_par2O=core::id::AtomID( ref_pose.residue( par2O.rsd() ).atom_index( aname ), par2O.rsd() );
  } else {
    tr.Warning << "cannot fix internal coords of " << atom << " in RigidChunk because 2nd parent is missing " << std::endl;
    return;
  }
  runtime_assert( ref_pose.conformation().atom_tree().has( ref_par1O ) );
  runtime_assert( ref_pose.conformation().atom_tree().has( ref_par2O ) );
  runtime_assert( pose.conformation().atom_tree().has( par1O ) );
  runtime_assert( pose.conformation().atom_tree().has( par2O ) );

  core::Real angle( numeric::angle_radians( ref_pose.xyz( ref_atom ), ref_pose.xyz( ref_par1O ), ref_pose.xyz( ref_par2O ) ) );
  tr.Trace << "ref angle direct: " << angle << std::endl;
  pose.conformation().set_bond_angle(  par2O, par1O, atom, angle );

  DOF_ID torsion_offset_dof( atom, PHI );
  DOF_ID ref_torsion_offset_dof( ref_atom, PHI );
  core::Real value( ref_pose.conformation().atom_tree().dof( ref_torsion_offset_dof ) );
  pose.conformation().set_dof( torsion_offset_dof, value );
}

void fix_mainchain_connect( core::pose::Pose& pose,
                            core::pose::Pose const& ref_pose,
                            core::Size upper_residue ) {
  core::conformation::Residue const & prev_rsd( ref_pose.residue( upper_residue-1 ) );
  core::conformation::Residue const &      rsd( ref_pose.residue( upper_residue ) );
  core::Size const nbb_prev( prev_rsd.n_mainchain_atoms() );
  core::id::AtomID bbM1   ( prev_rsd.mainchain_atom( nbb_prev-2 ),  upper_residue-1 );
  core::id::AtomID bb0    ( prev_rsd.mainchain_atom( nbb_prev-1 ),  upper_residue-1 );
  core::id::AtomID bb1    ( prev_rsd.mainchain_atom( nbb_prev   ),  upper_residue-1 );
  core::id::AtomID bb2    (      rsd.mainchain_atom(        1   ),  upper_residue   );
  core::id::AtomID bb3    (      rsd.mainchain_atom(        2   ),  upper_residue   );
  core::id::AtomID bb4    (      rsd.mainchain_atom(        3   ),  upper_residue   );

  core::conformation::Residue const & ref_resi = ref_pose.residue( upper_residue );
  tr.Trace << "mainchain torsion: ref: " << ref_resi.mainchain_torsion( 1 ) << " atom-tree: "
  << ref_pose.conformation().torsion_angle( bb1, bb2, bb3, bb4 ) << std::endl;

  core::conformation::Residue const & resi = pose.residue( upper_residue );
  tr.Trace << "mainchain torsion (before): conf: " << resi.mainchain_torsion( 1 ) << " atom-tree: "
  << pose.conformation().torsion_angle( bb1, bb2, bb3, bb4 ) << std::endl;

  pose.conformation().set_bond_length( bb1, bb2, ref_pose.conformation().bond_length( bb1, bb2 ) );
  pose.conformation().set_bond_angle ( bb0, bb1, bb2, ref_pose.conformation().bond_angle( bb0, bb1, bb2 ) );
  pose.conformation().set_bond_angle ( bb1, bb2, bb3, ref_pose.conformation().bond_angle( bb1, bb2, bb3 ) );
  pose.conformation().set_torsion_angle( bbM1, bb0, bb1, bb2, ref_pose.conformation().torsion_angle( bbM1, bb0, bb1, bb2 ) );
  pose.conformation().set_torsion_angle( bb0, bb1, bb2, bb3, ref_pose.conformation().torsion_angle( bb0, bb1, bb2, bb3 ) );
  pose.conformation().set_torsion_angle( bb1, bb2, bb3, bb4, ref_pose.conformation().torsion_angle( bb1, bb2, bb3, bb4 ) );

  core::conformation::Residue const & new_resi = pose.residue( upper_residue ); //this should trigger update of coords and torsions
  tr.Trace << "mainchain torsion (after): conf: " << new_resi.mainchain_torsion( 1 ) << " atom-tree: "
  << pose.conformation().torsion_angle( bb1, bb2, bb3, bb4 ) << std::endl;

  if ( prev_rsd.has( "O" ) ) {
    core::id::AtomID ref_atomO( prev_rsd.atom_index( "O" ), upper_residue-1 );
    core::id::AtomID atomO( pose.residue_type( upper_residue-1 ).atom_index( "O" ), upper_residue-1 );
    fix_internal_coords_of_siblings( pose, ref_pose, atomO, ref_atomO );
  }
  if ( rsd.has( "H" ) ) {
    core::id::AtomID ref_atomH( rsd.atom_index( "H" ), upper_residue );
    core::id::AtomID atomH( new_resi.atom_index( "H" ), upper_residue );
    runtime_assert( new_resi.has( "H" ) );
    fix_internal_coords_of_siblings( pose, ref_pose, atomH, ref_atomH );
  }

  if ( tr.Trace.visible() ) {
    bool ideal1( core::pose::is_ideal_position( upper_residue, ref_pose ) );
    if ( ideal1 && !core::pose::is_ideal_position( upper_residue, pose ) ) {
      tr.Warning << " pose in RigidChunkClaimer is not ideal at position " << upper_residue << " although template pose was ideal there " << std::endl;
    }

    bool ideal2( core::pose::is_ideal_position( upper_residue-1, ref_pose ) );
    if ( ideal2 && !core::pose::is_ideal_position( upper_residue-1, pose ) ) {
      tr.Warning << " pose in RigidChunkClaimer is not ideal at position " << upper_residue-1 << " although template pose was ideal there " << std::endl;
    }
  }
}

void copy_internal_coords( core::pose::Pose& pose, core::pose::Pose const& ref_pose, loops::Loops core ) {
  ///fpd if there are post modifications to pose (not in ref_pose), we can't just copy ref_pose->pose
  ///fpd    instead ... make xyz copy in rigid regions
  for ( loops::Loops::const_iterator region = core.begin(); region != core.end(); ++region ) {
    for ( core::Size i=region->start(); i<=region->stop(); ++i) {
      core::conformation::Residue const &rsd_i = ref_pose.residue(i);
      pose.replace_residue ( i , rsd_i , false );
    }
  }

  if ( tr.Trace.visible() ) {
    tr.Trace << pose.fold_tree() << std::endl;
    tr.Trace << ref_pose.fold_tree() << std::endl;
  }

  ///fpd fix connections
  ///fpd this requires that the input pose have one flanking residue on each side of each region
  for ( loops::Loops::const_iterator region = core.begin(); region != core.end(); ++region ) {
    core::Size loop_start = region->start();
    core::Size loop_stop  = region->stop();

    bool lower_connect = ( loop_start > 1
                          && !pose.residue(loop_start).is_lower_terminus()
                          && !pose.fold_tree().is_cutpoint( loop_start - 1 ) );
    bool upper_connect = ( loop_stop < pose.total_residue()
                          && !pose.residue(loop_stop).is_upper_terminus()
                          && !pose.fold_tree().is_cutpoint( loop_stop ) );

    if ( lower_connect ) {
      tr.Trace << "fixing lower connection for " << loop_start << std::endl;
      fix_mainchain_connect( pose, ref_pose, loop_start );
    } else {
      tr.Trace << "NOT fixing lower connection for " << loop_start << std::endl;
    }

    if ( upper_connect ) {
      tr.Trace << "fixing upper connection for " << loop_stop << std::endl;
      fix_mainchain_connect( pose, ref_pose, loop_stop+1 );
    } else {
      tr.Trace << "NOT fixing upper connection for " << loop_stop << std::endl;
    }
  }
}
//////////////////////////////// END MAGIC FPD CODE ////////////////////////////////////////////

void RigidChunkCM::initialize( Pose& pose ){

  if ( missing_density( pose, rigid_core_ ) ) {
    throw utility::excn::EXCN_BadInput( " missing density in backbone of rigid-chunk. Check your LOOP definitions.");
  }

  DofUnlock activation( pose.conformation(), passport() );

  for ( loops::Loops::const_iterator region = rigid_core_.begin();
        region != rigid_core_.end(); ++region ) {
    for (Size i=region->start(); i<=region->stop(); ++i) {
      core::conformation::Residue const &rsd_i = template_->residue(i);
      pose.replace_residue ( i , rsd_i , false );
    }
  }

  if ( tr.Trace.visible() ) {
    tr.Trace << pose.fold_tree() << std::endl;
    tr.Trace << template_->fold_tree() << std::endl;
  }

  ///fpd fix connections
  ///fpd this requires that the input pose have one flanking residue on each side of each region
  for ( loops::Loops::const_iterator region = rigid_core_.begin();
        region != rigid_core_.end(); ++region ) {
    Size loop_start = region->start();
    Size loop_stop  = region->stop();

    bool lower_connect = ( loop_start > 1
                          && !pose.residue(loop_start).is_lower_terminus()
                          && !pose.fold_tree().is_cutpoint( loop_start - 1 ) );
    bool upper_connect = ( loop_stop < pose.total_residue()
                          && !pose.residue(loop_stop).is_upper_terminus()
                          && !pose.fold_tree().is_cutpoint( loop_stop ) );

    if ( lower_connect ) {
      tr.Trace << "fixing lower connection for " << loop_start << std::endl;
      fix_mainchain_connect( pose, *template_, loop_start );
    } else {
      tr.Trace << "NOT fixing lower connection for " << loop_start << std::endl;
    }

    if ( upper_connect ) {
      tr.Trace << "fixing upper connection for " << loop_stop << std::endl;
      fix_mainchain_connect( pose, *template_, loop_stop+1 );
    } else {
      tr.Trace << "NOT fixing upper connection for " << loop_stop << std::endl;
    }
  }
}


void RigidChunkCM::apply( core::pose::Pose& ){
}

loops::Loops RigidChunkCM::select_parts( loops::Loops const& rigid_core, core::Size random_grow_loops_by ) {
  loops::Loops current_rigid_core;

  if( rigid_core.size() < 1){
    tr.Error << "Given rigid core had size < 1. Check your loop definitions." << std::endl;
    throw utility::excn::EXCN_BadInput( "Given rigid core had size < 1. Check your loop definitions." );
  }

  for( Size attempts = 1; attempts <= 50 && current_rigid_core.size() != 0; ++attempts ) {
    for ( loops::Loops::const_iterator it = rigid_core.begin(); it != rigid_core.end(); ++it ) {
      if ( RG.uniform() >= it->skip_rate() )  {
        current_rigid_core.push_back( *it );
      }
    }
  }

  if ( current_rigid_core.size() == 0 ) {
    current_rigid_core = rigid_core;
  }

  if ( random_grow_loops_by > 0 ) {
    core::Size nres( current_rigid_core[ current_rigid_core.size() ].stop() + 200 ); //it doesn't matter for this where exactly nres is.
    loops::Loops loops( current_rigid_core.invert( nres ) );
    loops.grow_all_loops( nres, random_grow_loops_by );
    tr.Info << "Enlarged loops: " << std::endl;
    tr.Info << loops << std::endl;
    current_rigid_core = loops.invert( nres );
  }

  return current_rigid_core;
}

std::string RigidChunkCM::get_name() const {
  return "RigidChunkCM";
}

moves::MoverOP RigidChunkCM::clone() const {
  return new RigidChunkCM( *this );
}

} // abscript
} // abinitio
} // protocols
