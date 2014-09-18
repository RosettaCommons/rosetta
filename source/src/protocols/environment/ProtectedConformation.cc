// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/ProtectedConformation.cc
/// @author Justin R. Porter

// Unit Headers
#include <protocols/environment/ProtectedConformation.hh>

// Package headers
#include <core/environment/DofPassport.hh>
#include <protocols/environment/Environment.hh>
#include <protocols/environment/EnvExcn.hh>

// Project headers
#include <core/conformation/symmetry/SymmetricConformation.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/BondedAtom.hh>

#include <core/id/TorsionID.hh>

#include <core/kinematics/AtomTree.hh>

// tracer
#include <boost/bind/bind.hpp>
#include <boost/ref.hpp>
#include <basic/Tracer.hh>

// C++ Headers
#include <algorithm>

// ObjexxFCL Headers

static thread_local basic::Tracer tr( "protocols.environment.ProtectedConformation", basic::t_info );

namespace protocols {
namespace environment {

using namespace core; 

using core::environment::DofPassport;
using core::environment::DofPassportCOP;
using core::conformation::Conformation;
using core::environment::SequenceAnnotation;
using core::environment::SequenceAnnotationCOP;

std::string get_mover_name( std::stack< DofPassportCOP > const& unlocks ){
  if( unlocks.empty() ){
    return "UNKNOWN/NONE (no active passports)";
  } else {
    return unlocks.top()->mover();
  }
}

ProtectedConformation::ProtectedConformation( EnvironmentCAP env, Conformation const& src ):
  Conformation( src ),
  env_( env ),
  environment_exists_( true )
{
  environment()->pconf_creation( this );
}

ProtectedConformation::ProtectedConformation( ProtectedConformation const& src ):
  Conformation( src ),
  unlocks_( src.unlocks_ ),
  annotations_( src.annotations_ ),
  env_( src.environment() ),
  environment_exists_( true )
{
  environment()->pconf_creation( this );
}

ProtectedConformation::~ProtectedConformation() {
  if( environment_exists_ ){
    environment()->pconf_destruction( this );
  }
}

EnvironmentCAP ProtectedConformation::environment() const {
  if( environment_exists_ ){
    return env_;
  } else {
    throw utility::excn::EXCN_NullPointer( "ProtectedEnvironment asked for Environment pointer after environment has already been destroyed/gone out of scope." );
  }
}

core::conformation::ConformationOP ProtectedConformation::clone() const {
  return new ProtectedConformation( *this );
}

/////////////////////////////////////////////////////////////////////////////
//// SECURITY OVERLOADS 	///////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
Conformation& ProtectedConformation::operator=( Conformation const& src ){
  ProtectedConformationCOP conf = dynamic_cast< ProtectedConformation const* >( &src );

  unlocks_ = conf->unlocks_;
  annotations_ = conf->annotations_;
  env_ = conf->environment();
  return Parent::operator=( src );
}

void ProtectedConformation::set_torsion( core::id::TorsionID const & id, core::Real const setting ){
  if( verify( id ) ){
    Parent::set_torsion( id, setting );
  } else {
    fail_verification( utility::to_string(id) );
  }
}

void ProtectedConformation::set_jump( int const jump_number, core::kinematics::Jump const& new_jump ){
  set_jump( jump_atom_id( jump_number ), new_jump );
}

void ProtectedConformation::set_jump( core::id::AtomID const& id, core::kinematics::Jump const& new_jump ){
  if( verify_jump( id ) ){
    Parent::set_jump( id, new_jump );
  } else {
    fail_verification( utility::to_string(id) );
  }
}

void ProtectedConformation::set_secstruct( Size const seqpos, char const setting ){
  if( verify_backbone( seqpos ) ){
    return Parent::set_secstruct( seqpos, setting );
  }
  fail_verification( "secondary structure assignment at position "+utility::to_string( seqpos ) );
}

void bond_lengths( std::map< core::id::DOF_ID, core::Real >& dofs,
                   core::conformation::ConformationOP conf,
                   core::Size seqpos ){
  using namespace core::id;
  using namespace core::kinematics::tree;

  //Collect all internal (and one external upstream) bond length.
  for( core::Size i = 1; i <= conf->residue( seqpos ).natoms(); ++i ){
    AtomID a_id( i, seqpos );
    core::kinematics::tree::Atom const& atom = conf->atom_tree().atom( a_id );

    if( !atom.is_jump() ){
      AtomCOP parent = atom.parent();
      if( parent && !parent->is_jump() ){
          dofs[ DOF_ID( a_id, THETA ) ] = conf->bond_length( atom.id(),
                                                             parent->id() );
      }
    }

    //Collect the (probably one) dependent downstream bond length.
    for( core::Size j = 1; j <= atom.n_children(); ++j ){
      AtomCOP child = atom.child( j );
      if( !child->is_jump() &&
          child->id().rsd() != a_id.rsd() ){
            dofs[ DOF_ID( child->id(), THETA ) ] = conf->bond_length( atom.id(),
                                                                      child->id() );
      }
    }
  }
}

void bond_angles( std::map< core::id::DOF_ID, core::Real >& dofs,
                  core::conformation::ConformationOP conf,
                  core::Size seqpos ){
  using namespace core::id;
  using namespace core::kinematics::tree;

  //Collect all internal (and a couple upstream external) bond angles.
  for( core::Size i = 1; i <= conf->residue( seqpos ).natoms(); ++i ){
    AtomID a_id( i, seqpos );
    core::kinematics::tree::Atom const& atom = conf->atom_tree().atom( a_id );
    if( !atom.is_jump() ){
      AtomCOP parent = atom.parent();
      if( parent && !parent->is_jump() ){
        if( parent->parent() ){
          dofs[ DOF_ID( a_id, THETA ) ] = conf->bond_angle( atom.id(),
                                                                   parent->id(),
                                                                   parent->parent()->id() );
        }
      }
    }

    //Collect dependent downstream bond angles.
    for( core::Size j = 1; j <= atom.n_children(); ++j ){
      AtomCOP child = atom.child( j );
      if( !child->is_jump() && child->id().rsd() != a_id.rsd() ){
        for( core::Size k = 1; k <= child->n_children(); ++k ){
          if( !child->child(k)->is_jump() ){
            AtomCOP subchild = child->child(k);
            dofs[ DOF_ID( subchild->id(), THETA ) ] = conf->bond_angle( atom.id(),
                                                                               child->id(),
                                                                               subchild->id() );
          }
        }
      }
    }

  }
}

void bond_torsions( std::map< core::id::DOF_ID, core::Real >& dofs,
                    core::conformation::ConformationOP conf,
                    core::Size seqpos ){
  using namespace core::id;

  dofs[ conf->dof_id_from_torsion_id( TorsionID( seqpos, BB, psi_torsion ) ) ] =
    conf->torsion( TorsionID( seqpos, BB, psi_torsion ) );
  dofs[ conf->dof_id_from_torsion_id( TorsionID( seqpos, BB, psi_torsion ) ) ] =
    conf->torsion( TorsionID( seqpos, BB, psi_torsion ) );
  dofs[ conf->dof_id_from_torsion_id( TorsionID( seqpos, BB, omega_torsion ) ) ] =
    conf->torsion( TorsionID( seqpos, BB, omega_torsion ) );

  if( conf->is_fullatom() ){
    for( core::Size i = 1; i <= conf->residue( seqpos ).chi().size(); ++i ){
      TorsionID id( seqpos, CHI, i );
      dofs[ conf->dof_id_from_torsion_id( id ) ] = conf->torsion( id );
    }
  }
}

void jump_dofs( std::map< core::id::DOF_ID, core::Real >& dofs,
                core::conformation::ConformationOP conf,
                core::Size seqpos ){

  using namespace core::id;
  using namespace core::kinematics::tree;


  for( core::Size i = 1; i <= conf->residue( seqpos ).natoms(); ++i ){
    core::kinematics::tree::Atom const& atom = conf->atom_tree().atom( AtomID( i, seqpos ) );
    if( atom.is_jump() ){
      dofs[ DOF_ID( atom.id(), RB1 ) ] =  conf->dof( DOF_ID( atom.id(), RB1 ) );
      dofs[ DOF_ID( atom.id(), RB2 ) ] =  conf->dof( DOF_ID( atom.id(), RB2 ) );
      dofs[ DOF_ID( atom.id(), RB3 ) ] =  conf->dof( DOF_ID( atom.id(), RB3 ) );
      dofs[ DOF_ID( atom.id(), RB4 ) ] =  conf->dof( DOF_ID( atom.id(), RB4 ) );
      dofs[ DOF_ID( atom.id(), RB5 ) ] =  conf->dof( DOF_ID( atom.id(), RB5 ) );
      dofs[ DOF_ID( atom.id(), RB6 ) ] =  conf->dof( DOF_ID( atom.id(), RB6 ) );
    }

    for( core::Size j = 1; j <= atom.n_children(); ++j ){
      AtomCOP child = atom.child(j);
      if( child->is_jump() ){
        dofs[ DOF_ID( child->id(), RB1 ) ] =  conf->dof( DOF_ID( child->id(), RB1 ) );
        dofs[ DOF_ID( child->id(), RB2 ) ] =  conf->dof( DOF_ID( child->id(), RB2 ) );
        dofs[ DOF_ID( child->id(), RB3 ) ] =  conf->dof( DOF_ID( child->id(), RB3 ) );
        dofs[ DOF_ID( child->id(), RB4 ) ] =  conf->dof( DOF_ID( child->id(), RB4 ) );
        dofs[ DOF_ID( child->id(), RB5 ) ] =  conf->dof( DOF_ID( child->id(), RB5 ) );
        dofs[ DOF_ID( child->id(), RB6 ) ] =  conf->dof( DOF_ID( child->id(), RB6 ) );
      }
    }
  }

}

/// @brief The plan here is to collect all the DoFs in the residue at 'seqpos' and on either side.
std::map< core::id::DOF_ID, core::Real > collect_dofs( core::Size const seqpos, core::conformation::ConformationCOP conf ){
  using namespace core::id;

  typedef std::map< core::id::DOF_ID, core::Real > DofMap;

  DofMap old_dofs;

  for( core::Size i = std::max( (Size) 1, seqpos-1 ); i <= std::min( conf->size(), seqpos+1 ); ++i ) {
    for( core::Size j = 1; j <= conf->residue(i).natoms(); ++j ){
      AtomID a_id = AtomID( j, i );

      if( conf->atom_tree().atom( a_id ).is_jump() ) {
        // If the atom's a jump atom, store the RB DoFs
        for( core::Size rbi = core::id::RB1; rbi <= core::id::RB6; ++rbi ){
          old_dofs[ DOF_ID( a_id, core::id::DOF_Type( rbi ) ) ] =
          conf->dof( DOF_ID( a_id, core::id::DOF_Type( rbi ) ) );
        }
      } else {
        // If the atom's a normal bonded atom, store the bond length/angle/torsions.
        DOF_ID d    ( AtomID( j, i ), core::id::D );
        DOF_ID theta( AtomID( j, i ), core::id::THETA );
        DOF_ID phi  ( AtomID( j, i ), core::id::PHI );

        old_dofs[ d ] = conf->dof( d );
        old_dofs[ theta ]  = conf->dof( theta );
        old_dofs[ phi ] = conf->dof( phi );
      }
    }
  }

  return old_dofs;
}

std::string dof_id_to_string( core::id::DOF_ID const& id, ProtectedConformation const& conf ){
  std::ostringstream ss;
  if( id.type() == core::id::PHI ){
    ss << "Torsion Angle ";
  } else if ( id.type() == core::id::THETA ){
    ss << "Bond Angle ";
  } else if ( id.type() == core::id::D ) {
    ss << "Bond Length ";
  } else if ( id.type() >= core::id::RB1 &&
             id.type() <= core::id::RB6 ){
    ss << "RB" << id.type() - core::id::RB1 + 1;
  }
  ss << " on " << conf.residue( id.rsd() ).name3() << id.rsd() << " owned by atom " << id.atomno() << "("
  << conf.residue( id.rsd() ).atom_name( (int) id.atomno() ) << ")";
  return ss.str();
}

template< typename Param >
void ProtectedConformation::replace_residue_sandbox( Size const seqpos, Residue const & new_rsd, Param p ){
  // TODO: accommodate design?

  if( new_rsd.natoms() != this->residue( seqpos ).natoms() ){
    std::ostringstream ss;
    ss << "Residue replacement of residue " << residue( seqpos ).name3() << seqpos << " ("
       << residue( seqpos ).natoms() << " atoms) failed because the input residue ("
       << new_rsd.name3() << ", " << new_rsd.natoms() << " atoms) differed in the number of atoms from the current residue.";
    throw EXCN_Env_Security_Exception( ss.str(), get_mover_name( unlocks_ ), environment() );
  }

  typedef std::map< core::id::DOF_ID, core::Real > DofMap;
  DofMap old_dofs = collect_dofs( seqpos, this );
  Residue const old_rsd( residue( seqpos ) );

  //do the actual replacement
  Parent::replace_residue( seqpos, new_rsd, p );

  DofMap new_dofs = collect_dofs( seqpos, this );;

  try{
    // Check to see that the passport was respected.
    if( old_dofs.size() != new_dofs.size() ){
      std::ostringstream ss;
      ss << "residue replacement of residue " << residue( seqpos ).name3() << seqpos << " with a "
         << new_rsd.name3() << " (" << new_rsd.natoms() << " failed"
         << ": dofs were created or destroyed during replacement (old size: " << old_dofs.size()
         << "; new size: " << new_dofs.size() << ")" << std::endl;
      throw EXCN_Env_Security_Exception( ss.str(), get_mover_name( unlocks_ ), environment() );
    }

    for( DofMap::iterator it = old_dofs.begin();
        it != old_dofs.end(); ++it ){

      core::id::DOF_ID const& dof = it->first;

      // Use a linear difference for D and RB1-6 and angular difference for PHI and THETA
      core::Real const angular_delta = std::cos( old_dofs[ dof ] ) - std::cos( new_dofs[ dof ] );
      core::Real const linear_delta = old_dofs[ dof ] - new_dofs[ dof ];
      core::Real const delta = dof.type() <= core::id::THETA ? angular_delta : linear_delta;

      if( std::abs( delta ) > 1e-6 && !verify( dof ) ){
        std::ostringstream ss;
        ss << dof_id_to_string( dof, *this ) << " in replace_residue (dof delta " << delta << ")";
        fail_verification( ss.str() );
      }
    }
  } catch( EXCN_Env_Security_Exception const& e ){
    //If we found an error, reset all the dofs to the previous values.
    Parent::replace_residue( seqpos, old_rsd, false );

    // re-throw the exception
    throw e;
  }
}

void ProtectedConformation::replace_residue( Size const seqpos, Residue const & new_rsd, bool const orient_backbone ) {
  replace_residue_sandbox( seqpos, new_rsd, orient_backbone );
}

void ProtectedConformation::replace_residue( Size const seqpos, core::conformation::Residue const& new_rsd,
                                             utility::vector1< std::pair< std::string, std::string > > const & atom_pairs ){
  // This is ok, since the parent function contains a call to the above function.
  Parent::replace_residue( seqpos, new_rsd, atom_pairs );
  //replace_residue_sandbox( seqpos, new_rsd, atom_pairs );
}

void ProtectedConformation::set_stub_transform( core::id::StubID const& stub_id1,
                                                core::id::StubID const& stub_id2,
                                                core::kinematics::RT const& target_rt ){
  int direction = 0;
  if( verify_jump( atom_tree().get_jump_atom_id( stub_id1, stub_id2, direction ) ) ){
    return Parent::set_stub_transform( stub_id1, stub_id2, target_rt );
  }
  fail_verification( "Stubs (AtomID="+utility::to_string( atom_tree().get_jump_atom_id( stub_id1, stub_id2, direction ) )
                     +", Stub1="+utility::to_string(stub_id1)+", Stub2="+utility::to_string(stub_id2) +")" );
}

void ProtectedConformation::set_dof( DOF_ID const& id, core::Real const setting ){
  if( verify( id ) ){
    return Parent::set_dof( id, setting);
  } else {
    fail_verification( dof_id_to_string( id, *this ) );
  }
}

void ProtectedConformation::set_torsion_angle( AtomID const & atom1,
                                               AtomID const & atom2,
                                               AtomID const & atom3,
                                               AtomID const & atom4,
                                               core::Real const setting,
                                               bool quiet){
  DOF_ID id = atom_tree().torsion_angle_dof_id( atom1, atom2, atom3, atom4, quiet );
  if( verify( id ) ){
    return Parent::set_torsion_angle(atom1, atom2, atom3, atom4, setting, quiet);
  } else {
    fail_verification( dof_id_to_string( id, *this )+" (via set_torsion_angle)" );
  }
}

void ProtectedConformation::set_bond_angle( AtomID const & atom1,
                                            AtomID const & atom2,
                                            AtomID const & atom3,
                                            core::Real const setting ){
  core::Real offset = 0.0;
  DOF_ID id = atom_tree().bond_angle_dof_id( atom1, atom2, atom3, offset );
  if( verify( id ) ){
    return Parent::set_bond_angle(atom1, atom2, atom3, setting);
  } else {
    fail_verification( dof_id_to_string( id, *this ) + " (via set_bond_angle)" );
  }
}

void ProtectedConformation::set_bond_length( AtomID const & atom1,
                                             AtomID const & atom2,
                                             core::Real const setting ){
  DOF_ID id = atom_tree().bond_length_dof_id( atom1, atom2 );
  if( verify( id ) ){
    return Parent::set_bond_length(atom1, atom2, setting);
  } else {
    fail_verification( dof_id_to_string( id, *this )+" (via set_bond_angle)" );
  }
}

void ProtectedConformation::insert_fragment( core::id::StubID const & /* instub_id */,
    FragRT const & /* outstub_transforms */,
    FragXYZ const & /* frag_xyz */ )
{
  // TODO: implement
  runtime_assert( false );
}


// ANNOTATION-RELATED METHODS
void ProtectedConformation::attach_annotation( SequenceAnnotationCOP annotations ){
  annotations_ = annotations;
}

SequenceAnnotationCOP ProtectedConformation::resolver() const {
  return annotations_;
}

// MISC OVERRIDES
bool ProtectedConformation::same_type_as_me( Conformation const & other, bool recurse ) const {
  if( dynamic_cast< ProtectedConformation const* >( &other ) ){
    if( recurse ){
      return other.same_type_as_me( *this, false );
    } else {
      return true;
    }
  } else {
    return false;
  }
}


// Security housekeeping
void ProtectedConformation::set_environment( EnvironmentCAP env ){
  env_ = env;
}

void ProtectedConformation::push_passport( DofPassportCOP passport ){
  unlocks_.push( passport );
}

DofPassportCOP ProtectedConformation::pop_passport(){
  DofPassportCOP popped = unlocks_.top();
  unlocks_.pop();
  return popped;
}

bool ProtectedConformation::has_passport() const {
  return !unlocks_.empty();
}

// VERIFICATION CONVIENENCE METHODS

bool ProtectedConformation::verify_backbone( Size const& seqpos ){
  return ( verify( core::id::TorsionID( seqpos, core::id::BB, core::id::phi_torsion ) ) &&
           verify( core::id::TorsionID( seqpos, core::id::BB, core::id::psi_torsion ) ) &&
           verify( core::id::TorsionID( seqpos, core::id::BB, core::id::omega_torsion ) ) );
}

bool ProtectedConformation::verify_jump( core::id::AtomID const& a_id ){
  bool allow = true;

  for( Size rb_i = core::id::RB1; rb_i <= core::id::RB6; ++rb_i ){
    allow &= verify( core::id::DOF_ID( a_id, core::id::DOF_Type(rb_i) ) );
  }

  return allow;
}

bool ProtectedConformation::verify( TorsionID const& id ){
  core::id::DOF_ID const d_id = dof_id_from_torsion_id( id );
  return !d_id.valid() || verify( dof_id_from_torsion_id( id ) );
}

bool ProtectedConformation::verify( core::id::DOF_ID const& id ){
  return ( !unlocks_.empty() && unlocks_.top()->dof_access( *environment(), id ) );
}

void ProtectedConformation::fail_verification( std::string const& str ){
  throw EXCN_Env_Security_Exception( str, *this, get_mover_name( unlocks_ ), environment() );
}

// ALWAYS-FAILING SECURITY OVERLOADS
void ProtectedConformation::fold_tree( FoldTree const & ){
  fail_verification( "direct setting of foldtree" );
}

void ProtectedConformation::chain_endings( utility::vector1< Size > const & ) {
  fail_verification( "chain endings list" );
}

void ProtectedConformation::insert_chain_ending( Size const ) {
  fail_verification( "chain endings insertion" );
}

void ProtectedConformation::delete_chain_ending( Size const ) {
  fail_verification( "chain endings deletion" );
}

void ProtectedConformation::reset_chain_endings() {
  fail_verification( "chain endings reset" );
}

void ProtectedConformation::chains_from_termini() {
  fail_verification( "chain status from termini" );
}

void ProtectedConformation::append_residue_by_jump( core::conformation::Residue const&, Size const,
                                                   std::string const&, std::string const&,
                                                   bool const ) {
  fail_verification( "residue-by-jump append" );
}

void ProtectedConformation::append_polymer_residue_after_seqpos( core::conformation::Residue const&,
                                                                Size const, bool const ) {
  fail_verification( "residue-after-seqpos append" );
}

void ProtectedConformation::safely_append_polymer_residue_after_seqpos( core::conformation::Residue const&,
                                                                       Size const, bool const ) {
  fail_verification( "residue-after-seqpos append" );
}

void ProtectedConformation::prepend_polymer_residue_before_seqpos( core::conformation::Residue const&,
                                                                  Size const, bool const ) {
  fail_verification( "residue-before-seqpos prepend" );
}

void ProtectedConformation::safely_prepend_polymer_residue_before_seqpos( core::conformation::Residue const&,
                                                                         Size const, bool const ) {
  fail_verification( "residue-before-seqpos prepend" );
}

void ProtectedConformation::delete_polymer_residue( Size const ) {
  fail_verification( "polymer residue delete" );
}

void ProtectedConformation::delete_residue_slow( Size const ) {
  fail_verification( "slow delete residue" );
}

void ProtectedConformation::delete_residue_range_slow( Size const, Size const ) {
  fail_verification( "slow delete residue range" );
}

void ProtectedConformation::declare_chemical_bond( Size const, std::string const&, Size const,
                                                  std::string const& ) {
  fail_verification( "chemical bond declaration" );
}

void ProtectedConformation::insert_conformation_by_jump( Conformation const&, Size const, Size const,
                                                        Size const, Size const, std::string const&,
                                                        std::string const& ) {
  fail_verification( "conformation insertion (seriously?)" );
}

void ProtectedConformation::rebuild_polymer_bond_dependent_atoms( Size const ) {
  fail_verification( "polymer bond dependent atom rebuild"  );
}

void ProtectedConformation::insert_ideal_geometry_at_polymer_bond( Size const ) {
  fail_verification( "insert_ideal_geometry_at_polymer_bond" );
}

void ProtectedConformation::insert_ideal_geometry_at_residue_connection( Size const,
                                                                         Size const ) {
  fail_verification( "insert_ideal_geometry_at_residue_connection" );
}

void ProtectedConformation::set_polymeric_connection( Size, Size ){
  fail_verification( "set_polymeric_connection" );
}

void ProtectedConformation::fix_disulfides( utility::vector1< std::pair<Size, Size> > ){
  fail_verification( "fix_disulfides" );
}

void ProtectedConformation::set_xyz( AtomID const &, core::PointPosition const & )
{
  fail_verification( "set_xyz" );
}

void ProtectedConformation::batch_set_xyz( utility::vector1<AtomID> const &,
                                           utility::vector1< core::PointPosition > const & ){
  fail_verification( "batch_set_xyz" );
}

void ProtectedConformation::reset_move_data(){
  fail_verification( "reset_move_data" );
}

void ProtectedConformation::clear(){
  fail_verification( "clear" );
}

void ProtectedConformation::fill_missing_atoms( core::id::AtomID_Mask ){
  fail_verification( "fill_missing_atoms" );
}

} // environment
} // protocols
