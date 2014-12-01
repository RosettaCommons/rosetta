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
#include <protocols/environment/EnvExcn.hh>
#include <protocols/environment/ProtectedConformation.hh>

#include <protocols/environment/claims/CutBiasClaim.hh>
#include <protocols/environment/claims/JumpClaim.hh>
#include <protocols/environment/claims/XYZClaim.hh>
#include <protocols/environment/claims/EnvLabelSelector.hh>

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
#include <utility/excn/Exceptions.hh>

#include <numeric/random/random.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/xyz.functions.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

// tracer
#include <basic/Tracer.hh>

// C++ Headers
#include <algorithm>

// ObjexxFCL Headers

//Req'd on WN32
#include <basic/datacache/WriteableCacheableMap.hh>

static thread_local basic::Tracer tr( "protocols.abinitio.abscript.RigidChunkCM", basic::t_info );

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
  return protocols::moves::MoverOP( new RigidChunkCM );
}

std::string
RigidChunkCMCreator::mover_name() {
  return "RigidChunkCM";
}

RigidChunkCM::RigidChunkCM():
  Parent(),
  selector_( NULL )
{}

RigidChunkCM::RigidChunkCM(
  core::pack::task::residue_selector::ResidueSelectorCOP selector,
  core::pose::Pose const& template_pose
):
  Parent(),
  template_( core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose(template_pose) ) ) ),
  selector_( selector )
{
  loops::Loops regions_in;
  regions_in.add_loop( loops::Loop( 1, templ().total_residue() ) );
  rigid_core( regions_in );
}

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
                                 basic::datacache::DataMap& datamap,
                                 protocols::filters::Filters_map const&,
                                 protocols::moves::Movers_map const& movermap,
                                 core::pose::Pose const& ){

  using namespace basic::options;

  //This option must be provided, but is not always accessed. This is here to fail fast.
  tag->hasOption( "name" );

  std::string const SELECTOR = "selector";
  if( tag->hasOption( SELECTOR ) ){
    set_selector( datamap.get_ptr< core::pack::task::residue_selector::ResidueSelector const >( "ResidueSelector", tag->getOption<std::string>( SELECTOR ) ) );
  }

  std::string const APPLY_TO_TEMPLATE = "apply_to_template";
  utility::vector1< moves::MoverOP > apply_movers;
  if( tag->hasOption( APPLY_TO_TEMPLATE ) ){
    std::string const apply_to_templates = tag->getOption< std::string >( APPLY_TO_TEMPLATE );
    utility::vector1< std::string > movernames = utility::string_split( apply_to_templates, ',' );
    for( utility::vector1< std::string >::const_iterator movername = movernames.begin();
        movername != movernames.end(); ++movername ){
      apply_movers.push_back( movermap.find( *movername )->second );
    }
  }

  if( tag->hasOption("template") ){
    std::string file = tag->getOption< std::string >( "template" );
    if( file == "INPUT" ){
      if( tag->hasOption( APPLY_TO_TEMPLATE ) ){
        std::ostringstream ss;
        ss << "In " << this->get_name() << " the option '" << APPLY_TO_TEMPLATE
        << "' is not combinable with input templates." << std::endl;
        throw utility::excn::EXCN_RosettaScriptsOption( ss.str() );
      }
      template_ = NULL;
    } else {
      core::pose::PoseOP p( new core::pose::Pose() );
      core::import_pose::pose_from_pdb( *p,
                                       *core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ),
                                       file );

      for( Size i = 1; i <= apply_movers.size(); ++i ){
        apply_movers[i]->apply( *p );
        tr.Debug << "RigidChunkCM named " << tag->getOption< std::string >( "name" ) << " applied " << apply_movers[i]->get_name() << std::endl;
      }
      template_ = p;
    }
  } else {
    throw utility::excn::EXCN_BadInput( "RigidChunkCM requires a template pdb with coordinates for rigid regions.");
  }

  loops::Loops regions_in;
  if( tag->hasOption("region_file") ){
    regions_in = read_rigid_core( tag->getOption< std::string >( "region_file" ) );
  } else if ( tag->hasOption( "region" ) ){
    runtime_assert( false ); // unimplemented
  } else if ( tag->hasOption("loop_file") ){
    runtime_assert( false ); // unimplemented
  } else if ( selector() ){
    tr.Debug << "Using ALL residues in '" << tag->getOption< std::string >( SELECTOR )
             << "' in the absence of a specification (e.g. 'region_file')." << std::endl;
    regions_in.add_loop( loops::Loop( 1, templ().total_residue() ) );
  } else {
    std::string const err = "RigidChunkCM did not recieve any information about which chunks to fix. Use options 'region_file' to specify.";
    throw utility::excn::EXCN_BadInput( err );
  }

  rigid_core( select_parts( regions_in, tag->getOption( "random_grow_by",
                                                        (core::Size)( option[ OptionKeys::loops::random_grow_loops_by ]() ) ) ) );

}

void RigidChunkCM::configure(
  core::pose::Pose const& in_p,
  utility::vector1< bool > const sim_selection
) {

  if( !template_ ){
    tr.Debug << "Building template from broker-time pose." << std::endl;
    template_ = core::pose::PoseCOP( core::pose::PoseOP( new core::pose::Pose( in_p ) ) );
  }

  utility::vector1< bool > templ_selection( templ().total_residue() );
  rigid_core().transfer_to_residue_vector( templ_selection, true );

  if( templ_selection.index( true ) == 0 ){
    std::ostringstream ss;
    ss << this->get_name() << " reports that its input loops file (aka rigid core file) contained no residues."
       << " Check your input." << std::endl;
    throw utility::excn::EXCN_BadInput( ss.str() );
  } else if( sim_selection.index( true ) == 0 ) {
    std::ostringstream ss;
    ss << this->get_name() << " reports that its selector (used on the simulation to determine where to insert residues)"
    << " returned an empty selection. Check your input." << std::endl;
    throw utility::excn::EXCN_BadInput( ss.str() );
  }

  core::Size templ_pos = 1;
  core::Size sim_pos = 1;
  std::map< core::Size, core::Size > templ_target;
  std::map< core::Size, core::Size > sim_origin;
  while( templ_pos <= templ().total_residue() &&
         sim_pos <= in_p.total_residue() ){
    if( sim_selection[ sim_pos ] && templ_selection[ templ_pos ] ){
      // this position is selected in both, it's a correspondence.
      templ_target[ templ_pos ] = sim_pos;
      sim_origin[ sim_pos ] = templ_pos;
      templ_pos += 1;
      sim_pos += 1;
    } else if( sim_selection[ sim_pos ] ){
      // this position is only selected in the simulation. Advance the template one.
      templ_pos += 1;
    } else if( templ_selection[ templ_pos ] ){
      // this position is only selected in the template. Advance the simulation one.
      sim_pos += 1;
    } else { // neither
      templ_pos += 1;
      sim_pos += 1;
    }
  }

  this->templ_target( templ_target );
  this->sim_origin( sim_origin );

  for( core::Size i = templ_pos; i <= templ_selection.size(); ++i ){
    if( templ_selection[i] ){
      tr.Warning << "[WARNING] " << this->get_name() << " reports that "
                 << std::count( templ_selection.begin() + i - 1, templ_selection.end(), true )
                 << " residues beginning at " << templ().residue( i ).name3() << i
                 << "in the template do not fit in the selection given. "
                 << "This may or may not be problematic." << std::endl;
      break;
    }
  }
  for( core::Size i = sim_pos; i <= sim_selection.size(); ++i ){
    if( sim_selection[i] ){
      tr.Warning << "[WARNING] " << this->get_name() << " reports that "
                 << std::count( sim_selection.begin() + i - 1, sim_selection.end(), true )
                 << " residues beginning at " << in_p.residue( i ).name3() << i
                 << " in the template do not fit in the selection given. "
                 << "This may or may not be problematic." << std::endl;
      break;
    }
  }
}

claims::EnvClaims RigidChunkCM::yield_claims( core::pose::Pose const& in_p,
                                              basic::datacache::WriteableCacheableMapOP ){
  using namespace claims;
  EnvClaims claims;

  utility::vector1< bool > selection( in_p.total_residue(), false );
  if( selector() ){
    selection = selector()->apply( in_p );
  } else {
    selection = utility::vector1< bool >( in_p.total_residue(), true );
  }
  configure( in_p, selection );

  ClaimingMoverOP this_ptr = utility::pointer::static_pointer_cast< ClaimingMover > ( get_self_ptr() );

  loops::Loops simulation_regions( selection );

  std::pair< core::Size, core::Size > prev_region = std::make_pair( 0, 0 );
  for ( loops::Loops::const_iterator loop_it = simulation_regions.begin();
        loop_it != simulation_regions.end(); ++loop_it ) {

    // Calculate the area that needs to be claimed.
    core::Size excl_region_end = loop_it->stop();
    while( sim_origin().find( excl_region_end ) == sim_origin().end() ){
      assert( excl_region_end != loop_it->start() );
      assert( excl_region_end > 0 );
      excl_region_end -= 1;
    }
    std::pair< core::Size, core::Size > const excl_region = std::make_pair( loop_it->start(),
                                                                            excl_region_end );

    XYZClaimOP xyz_claim( new XYZClaim( this_ptr, "BASE", excl_region ) );

    xyz_claim->strength( EXCLUSIVE, EXCLUSIVE );
    xyz_claim->set_relative( true ); //we don't care where it's located in space; just that the relative locations are ok.
    claims.push_back( xyz_claim );
    tr.Debug << this->get_name() << ": built EXCLUSIVE XYZClaim for " << excl_region.first << "-" << excl_region.second
             << " in " << "BASE" << std::endl;

    // For initialization, we need to some of FDP's fixing, which requires some control outside the region
    std::pair< core::Size, core::Size > const supp_region = std::make_pair( std::max( Size( 1 ), excl_region.first-1 ),
                                                                            std::min( in_p.total_residue(), excl_region.second+1 ) );
    XYZClaimOP support_claim( new XYZClaim( this_ptr, "BASE", supp_region ) );
    support_claim->strength( MUST_CONTROL, MUST_CONTROL );
    claims.push_back( support_claim );
    tr.Debug << this->get_name() << ": built support XYZClaim for " << supp_region.first << "-" << supp_region.second
             << " in " << "BASE" << std::endl;

    claims.push_back( protocols::environment::claims::EnvClaimOP( new CutBiasClaim( this_ptr, "BASE", excl_region, 0.0 ) ) );

    if( prev_region.first != 0 && prev_region.second != 0 ){
      core::Size const jump_start = prev_region.first + ( ( prev_region.second - prev_region.first ) / 2 );
      core::Size const jump_end   = excl_region.first + ( ( excl_region.second - excl_region.first ) / 2 );

      JumpClaimOP j_claim( new JumpClaim( this_ptr,
                                          "RigidChunkJump"+utility::to_string( claims.size()/4 ),
                                          LocalPosition( "BASE", jump_start ),
                                          LocalPosition( "BASE", jump_end ) ) );
      j_claim->strength( EXCLUSIVE, EXCLUSIVE );
      j_claim->physical( false );

      claims.push_back( j_claim );
    }
    prev_region = excl_region;
  }

  return claims;
}
//
//bool missing_density( core::pose::Pose const& pose, loops::Loops const& core, core::Size region_offset ){
//  bool missing_density = false;
//
//  //sanity check: no missing density in backbon in any of the rigid_core residues?
//  for ( loops::Loops::const_iterator it = core.begin(); it!=core.end(); ++it ) {
//    for ( core::Size pos = it->start(); pos <=it->stop(); ++pos ) {
//      // Do we really have Sidechains ?
//      // check this my making sure that no SC atom is more than 20A (?) away from CA
//      numeric::xyzVector< core::Real> ca_pos = pose.residue( pos ).atom("CA").xyz();
//      numeric::xyzVector< core::Real> n_pos = pose.residue( pos ).atom("N").xyz();
//      numeric::xyzVector< core::Real> o_pos = pose.residue( pos ).atom("O").xyz();
//      if ( ( n_pos - ca_pos).length() > 20 || ( ( n_pos - o_pos ).length() > 20 ) ) {
//        tr.Error << "missing backbone in rigid-chunk at " << pos << std::endl;
//        missing_density = true;
//      }
//    }
//  }
//
//  return missing_density;
//}

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
  core::Size const offset = atom.rsd() - ref_atom.rsd();
  if ( has_par1 && ref_has_par1 ) {
    par1O=pose.conformation().atom_tree().atom( atom ).parent()->id();

    std::string const & aname( pose.residue( par1O.rsd() ).atom_name( par1O.atomno() ));
    core::Size const rsd_num = par1O.rsd() - offset;

    ref_par1O=core::id::AtomID( ref_pose.residue( rsd_num ).atom_index( aname ), rsd_num );
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
    core::Size const rsd_num = par2O.rsd() - offset;

    ref_par2O=core::id::AtomID( ref_pose.residue( rsd_num ).atom_index( aname ), rsd_num );
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

class AtomPack {
public:
  AtomPack( core::Size upper_residue,
            core::conformation::Residue const& prev_rsd,
            core::conformation::Residue const& rsd ) :
    bbM1 ( prev_rsd.mainchain_atom( prev_rsd.n_mainchain_atoms()-2 ),  upper_residue-1 ),
    bb0  ( prev_rsd.mainchain_atom( prev_rsd.n_mainchain_atoms()-1 ),  upper_residue-1 ),
    bb1  ( prev_rsd.mainchain_atom( prev_rsd.n_mainchain_atoms()   ),  upper_residue-1 ),
    bb2  (      rsd.mainchain_atom( 1 ),  upper_residue   ),
    bb3  (      rsd.mainchain_atom( 2 ),  upper_residue   ),
    bb4  (      rsd.mainchain_atom( 3 ),  upper_residue   )
  {}
  core::id::AtomID bbM1;
  core::id::AtomID bb0;
  core::id::AtomID bb1;
  core::id::AtomID bb2;
  core::id::AtomID bb3;
  core::id::AtomID bb4;
};

void fix_mainchain_connect( core::pose::Pose& pose,
                            core::Size global_upper,
                            core::pose::Pose const& ref_pose,
                            core::Size local_upper ) {

  //These guys just hold the various AtomIDs required to correclty address the angles that need to be repaired.
  AtomPack const local = AtomPack( local_upper,
                                   ref_pose.residue( local_upper-1),
                                   ref_pose.residue( local_upper ) );
  AtomPack const global = AtomPack( global_upper,
                                    pose.residue( global_upper-1),
                                    pose.residue( global_upper ) );

  core::conformation::Residue const & ref_resi = ref_pose.residue( local_upper );
  tr.Trace << "mainchain torsion: ref: " << ref_resi.mainchain_torsion( 1 ) << " atom-tree: "
           << ref_pose.conformation().torsion_angle( local.bb1, local.bb2, local.bb3, local.bb4 ) << std::endl;

  core::conformation::Residue const & resi = pose.residue( global_upper );
  tr.Trace << "mainchain torsion (before): conf: " << resi.mainchain_torsion( 1 ) << " atom-tree: "
           << pose.conformation().torsion_angle( global.bb1, global.bb2, global.bb3, global.bb4 ) << std::endl;

  pose.conformation().set_bond_length( global.bb1, global.bb2,
                                       ref_pose.conformation().bond_length( local.bb1, local.bb2 ) );
  pose.conformation().set_bond_angle ( global.bb0, global.bb1, global.bb2,
                                       ref_pose.conformation().bond_angle( local.bb0, local.bb1, local.bb2 ) );
  pose.conformation().set_bond_angle ( global.bb1, global.bb2, global.bb3,
                                       ref_pose.conformation().bond_angle( local.bb1, local.bb2, local.bb3 ) );
  pose.conformation().set_torsion_angle( global.bbM1, global.bb0, global.bb1, global.bb2,
                                         ref_pose.conformation().torsion_angle( local.bbM1, local.bb0, local.bb1, local.bb2 ) );
  pose.conformation().set_torsion_angle( global.bb0, global.bb1, global.bb2, global.bb3,
                                         ref_pose.conformation().torsion_angle( local.bb0, local.bb1, local.bb2, local.bb3 ) );
  pose.conformation().set_torsion_angle( global.bb1, global.bb2, global.bb3, global.bb4,
                                         ref_pose.conformation().torsion_angle( local.bb1, local.bb2, local.bb3, local.bb4 ) );

  core::conformation::Residue const & new_resi = pose.residue( global_upper ); //this should trigger update of coords and torsions
  tr.Trace << "mainchain torsion (after): conf: " << new_resi.mainchain_torsion( 1 ) << " atom-tree: "
           << pose.conformation().torsion_angle( global.bb1, global.bb2, global.bb3, global.bb4 ) << std::endl;

  core::conformation::Residue const & prev_rsd( ref_pose.residue( local_upper-1 ) );
  core::conformation::Residue const &      rsd( ref_pose.residue( local_upper ) );

  if ( prev_rsd.has( "O" ) ) {
    core::id::AtomID ref_atomO( prev_rsd.atom_index( "O" ), local_upper-1 );
    core::id::AtomID atomO( pose.residue_type( global_upper-1 ).atom_index( "O" ), global_upper-1 );
    fix_internal_coords_of_siblings( pose, ref_pose, atomO, ref_atomO );
  }
  if ( rsd.has( "H" ) ) {
    core::id::AtomID ref_atomH( rsd.atom_index( "H" ), local_upper );
    core::id::AtomID atomH( new_resi.atom_index( "H" ), global_upper );
    runtime_assert( new_resi.has( "H" ) );
    fix_internal_coords_of_siblings( pose, ref_pose, atomH, ref_atomH );
  }

  if ( tr.Trace.visible() ) {
    bool ideal1( core::pose::is_ideal_position( local_upper, ref_pose ) );
    if ( ideal1 && !core::pose::is_ideal_position( global_upper, pose ) ) {
      tr.Warning << " pose in RigidChunkClaimer is not ideal at position " << global_upper << " although template pose was ideal there " << std::endl;
    }

    bool ideal2( core::pose::is_ideal_position( local_upper-1, ref_pose ) );
    if ( ideal2 && !core::pose::is_ideal_position( global_upper-1, pose ) ) {
      tr.Warning << " pose in RigidChunkClaimer is not ideal at position " << global_upper-1 << " although template pose was ideal there " << std::endl;
    }
  }
}
//////////////////////////////// END MAGIC FPD CODE ////////////////////////////////////////////

void RigidChunkCM::initialize( Pose& pose ){

  DofUnlock activation( pose.conformation(), passport() );

//  if ( missing_density( pose, rigid_core(), region_offset_ ) ) {
//    throw utility::excn::EXCN_BadInput( " missing density in backbone of rigid-chunk. Check your LOOP definitions.");
//  }

  core::pose::Pose reference( pose );

  for ( Size sim_pos = 1; sim_pos <= pose.total_residue(); ++sim_pos ) {

    if( sim_origin().find( sim_pos ) != sim_origin().end() ){
      Size const templ_pos = sim_origin().at( sim_pos );
      Residue const templ_res = templ().residue( templ_pos );

      try {
        tr.Trace << "Replacing simulation " << pose.residue( sim_pos ).name3() << sim_pos
                 << " with template " << templ_res.name3() << templ_pos << std::endl;

        pose.replace_residue( sim_pos, templ_res , false );

      } catch ( EXCN_Env_Security_Exception& e ) {
        std::ostringstream ss;

        ss << this->get_name() << " failed trying to replace the simulation "
           << pose.residue( sim_pos ).name3() << sim_pos << "(" << pose.residue( sim_pos ).natoms() << " atoms) with the template "
           << templ().residue( templ_pos ).name3() << templ_pos <<  " ( "
           << templ().residue( templ_pos ).natoms() << " atoms ).  ";
        if( pose.residue( sim_pos ).name3() != templ().residue( templ_pos ).name3() ){
          ss << "The problem is probably that the sequence shift is off. Surrounding sequences: " << std::endl
             << "template:   ";

          for( Size k = std::max( Size( 1 ), sim_origin().at( sim_pos ) - 5 ); k <= std::min( templ().total_residue(), sim_origin().at( sim_pos ) + 5 ); ++k ){ ss << templ().residue( k ).name1(); }
          ss << std::endl
             << "simulation: ";
          for( Size k = std::max( Size( 1 ), sim_pos - 5 ); k <= std::min( pose.total_residue(), sim_pos + 5 ); ++k ){ ss << pose.residue( k ).name1(); }
          ss << std::endl;
        } else if( ( pose.is_centroid() && !templ().is_centroid() ) ||
                   ( pose.is_fullatom() && !templ().is_fullatom() ) ){
          ss << "The problem is that the input is fullatom and it's being a applied to a centroid pose (or vice versa)--"
             << "note the differing atom counts. If so, the 'apply_to_template' option to apply a "
             << "SwitchResidueTypeSetMover to the template. There could also be a problem with differing variant." << std::endl;
        } else {
          ss << "The problem is that there are a differing number of atoms between the two residues, and replace_residue "
             << "doesn't know how to handle the difference (differing degrees of freedom and all). This is usually from "
             << "differing variants." << std::endl
             << "Simulation: " << utility::to_string( pose.residue( sim_pos ).type().properties().get_list_of_variants() ) << std::endl
             << "Template: " << utility::to_string( templ().residue( templ_pos ).type().properties().get_list_of_variants() ) << std::endl;
        }

        throw utility::excn::EXCN_BadInput( ss.str() );
      }
    }
  }

  if ( tr.Debug.visible() ) {
    tr.Debug << pose.fold_tree() << std::endl;
    tr.Debug << templ().fold_tree() << std::endl;
  }

  // correct mainchain connections

  for( Size sim_pos = 1; sim_pos <= pose.total_residue(); ++sim_pos ){

    // If the residue before this is not in the region
    if( sim_origin().find( sim_pos - 1 ) == sim_origin().end() &&
        sim_origin().find( sim_pos ) != sim_origin().end() ){
      core::Size const templ_pos = sim_origin().at( sim_pos );

      bool lower_connect = ( sim_pos > 1
                            && !pose.residue( sim_pos ).is_lower_terminus()
                            && !pose.fold_tree().is_cutpoint( sim_pos - 1 ) );
      if ( lower_connect && templ_pos - 1 < 1 ){
        tr.Debug << "fixing lower connection for " << sim_pos << " using non-template values." << std::endl;
        fix_mainchain_connect( pose, sim_pos, reference, sim_pos );
      } else if ( lower_connect ) {
        tr.Debug << "fixing lower connection for " << sim_pos << std::endl;
        fix_mainchain_connect( pose, sim_pos, templ(), templ_pos );
      } else {
        tr.Debug << "NOT fixing lower connection for " << sim_pos << std::endl;
      }

    }
    // If the residue after this is not in the region
    if( sim_origin().find( sim_pos + 1 ) == sim_origin().end() &&
        sim_origin().find( sim_pos ) != sim_origin().end() ){
      core::Size const templ_pos = sim_origin().at( sim_pos );

      bool upper_connect = ( sim_pos < pose.total_residue()
                            && !pose.residue( sim_pos ).is_upper_terminus()
                            && !pose.fold_tree().is_cutpoint( sim_pos ) );

      if ( upper_connect && templ_pos + 1 > templ().total_residue() ){
        tr.Debug << "fixing upper connection for " << sim_pos << " using non-template values." << std::endl;
        fix_mainchain_connect( pose, sim_pos+1, reference, sim_pos+1 );
      } else if ( upper_connect ) {
        tr.Debug << "fixing upper connection for " << sim_pos << std::endl;
        fix_mainchain_connect( pose, sim_pos+1, templ(), templ_pos+1 );
      } else {
        tr.Debug << "NOT fixing upper connection for " << sim_pos << std::endl;
      }
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
      if ( numeric::random::rg().uniform() >= it->skip_rate() )  {
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

void RigidChunkCM::set_selector(
  core::pack::task::residue_selector::ResidueSelectorCOP selector
) {
  if( Parent::state_check( __FUNCTION__, ( selector.get() == selector_.get() ) ) ){
    selector_ = selector;
  }
}

core::pack::task::residue_selector::ResidueSelectorCOP RigidChunkCM::selector() const {
  return selector_;
}

std::string RigidChunkCM::get_name() const {
  std::ostringstream ss;
  ss << "RigidChunkCM(";

  if( sim_origin().empty() ) {
    for ( loops::Loops::const_iterator region = rigid_core().begin();
         region != rigid_core().end(); ++region ) {
      ss << "template:" << region->start() << "-" << region->stop() << ",";
    }
  } else {
    ss << "simulation :";
    for( std::map< Size, Size >::const_iterator pos = sim_origin().begin(); pos != sim_origin().end(); ++pos ){
      if( sim_origin().find( pos->first - 1 ) == sim_origin().end() ){
        ss << pos->first << "-";
      }
      if( sim_origin().find( pos->first + 1 ) == sim_origin().end() ){
        ss << pos->first << ",";
      }
    }
  }
  ss << ")";

  return ss.str();
}

void RigidChunkCM::passport_updated() {
  if( !has_passport() ){
    templ_target_.clear();
    sim_origin_.clear();
  }
}

moves::MoverOP RigidChunkCM::clone() const {
  return moves::MoverOP( new RigidChunkCM( *this ) );
}

void RigidChunkCM::rigid_core(
  loops::Loops core_in
) {
  rigid_core_ = core_in;
}

loops::Loops const& RigidChunkCM::rigid_core() const {
  return rigid_core_;
}


} // abscript
} // abinitio
} // protocols
