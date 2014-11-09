// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/claims/TorsionClaim.cc
/// @brief Implementation file for TorsionClaim.
/// @author Justin R. Porter

// Unit Headers
#include <protocols/environment/claims/TorsionClaim.hh>

// Package Headers
#include <core/environment/LocalPosition.hh>

#include <core/pack/task/residue_selector/ResidueSelector.fwd.hh>

#include <protocols/environment/claims/EnvLabelSelector.hh>
#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/ClaimStrength.hh>

#include <protocols/environment/ClaimingMover.hh>
#include <protocols/environment/ProtectedConformation.hh>

// Project Headers
#include <core/id/TorsionID.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>

// C++ headers

// option key includes


static thread_local basic::Tracer tr( "protocols.environment.TorsionClaim", basic::t_info );

namespace protocols {
namespace environment {
namespace claims {

using core::environment::LocalPosition;
using core::environment::LocalPositions;
using core::conformation::Conformation;

TorsionClaim::TorsionClaim( ClaimingMoverOP owner,
                            core::pack::task::residue_selector::ResidueSelectorCOP selector ) :
  EnvClaim( owner ),
  selector_( selector ),
  c_str_( MUST_CONTROL ),
  i_str_( DOES_NOT_CONTROL ),
  claim_sidechain_( false ),
  claim_backbone_( true )
{}


TorsionClaim::TorsionClaim( ClaimingMoverOP owner,
                            utility::tag::TagCOP tag,
                            basic::datacache::DataMap& datamap ):
  EnvClaim( owner ),
  c_str_( Parent::parse_ctrl_str( tag->getOption< std::string >( "control_strength" ) ) ),
  i_str_( Parent::parse_ctrl_str( tag->getOption< std::string >( "initialization_strength", "DOES_NOT_CONTROL" ) ) )
{
  claim_backbone_ = tag->getOption< bool >( "backbone", false );
  claim_sidechain_ = tag->getOption< bool >( "sidechain", false );
  if( !claim_sidechain_ && !claim_backbone_ ){
    utility::tag::TagCOP tag_parent = tag->getParent().lock();
    tr.Warning << type() << " owned by mover named "
               << tag_parent->getOption< std::string >( "name" )
               << " was not configured to claim neither sidechains nor backbone torsions."
               << " Are you sure this is what you wanted?" << std::endl;
  }

  selector_ = datamap.get_ptr< core::pack::task::residue_selector::ResidueSelector const >( "ResidueSelector", tag->getOption<std::string>( "selector" ) );
}

TorsionClaim::TorsionClaim( ClaimingMoverOP owner,
                            LocalPosition const & local_pos):
  EnvClaim( owner ),
  selector_( core::pack::task::residue_selector::ResidueSelectorCOP( core::pack::task::residue_selector::ResidueSelectorOP( new EnvLabelSelector( local_pos ) ) ) ),
  c_str_( MUST_CONTROL ),
  i_str_( DOES_NOT_CONTROL ),
  claim_sidechain_( false ),
  claim_backbone_( true )
{}

TorsionClaim::TorsionClaim( ClaimingMoverOP owner,
                            std::string const & label,
                            std::pair< core::Size, core::Size > const & range ):
  EnvClaim( owner ),
  c_str_( MUST_CONTROL ),
  i_str_( DOES_NOT_CONTROL ),
  claim_sidechain_( false ),
  claim_backbone_( true )
{
  LocalPositions local_positions = LocalPositions();

  for( Size i = range.first; i <= range.second; ++i){
    local_positions.push_back( core::environment::LocalPositionOP( new LocalPosition( label, i ) ) );
  }

  selector_ = core::pack::task::residue_selector::ResidueSelectorCOP( core::pack::task::residue_selector::ResidueSelectorOP( new EnvLabelSelector( local_positions ) ) );
}

TorsionClaim::TorsionClaim( ClaimingMoverOP owner,
                            LocalPositions const & positions ):
  EnvClaim( owner ),
  selector_( core::pack::task::residue_selector::ResidueSelectorCOP( core::pack::task::residue_selector::ResidueSelectorOP( new EnvLabelSelector( positions ) ) ) ),
  c_str_( MUST_CONTROL ),
  i_str_( DOES_NOT_CONTROL ),
  claim_sidechain_( false ),
  claim_backbone_( true )
{}


DOFElement TorsionClaim::wrap_dof_id( core::id::DOF_ID const & id ) const {
  DOFElement e = Parent::wrap_dof_id( id );

  e.c_str = ctrl_strength();
  e.i_str = init_strength();

  return e;
}

void TorsionClaim::insert_dof_element( core::conformation::Conformation const & conf,
                                       DOFElements& elements,
                                       core::Size seqpos,
                                       core::id::TorsionType type,
                                       core::Size torsion_number) const {
  using namespace core::id;

  TorsionID const t_id( seqpos, type, torsion_number );
  assert( t_id.valid() );

  DOF_ID d_id = conf.dof_id_from_torsion_id( t_id );
  if( d_id.valid() ){
    elements.push_back( wrap_dof_id( d_id ) );
  }
}

void TorsionClaim::yield_elements( core::pose::Pose const & pose, DOFElements& elements ) const {

  using namespace core::id;

  utility::vector1< bool > subset = selector_->apply( pose );

  Conformation const & conf = pose.conformation();

  for( Size seqpos = 1; seqpos <= pose.total_residue(); ++seqpos ){
    if( !subset[seqpos] ){
      // If this residue isn't selected by the selector, we're not supposed to do anything.
      // Continue on to the next residue.
      continue;
    }

    if( conf.residue( seqpos ).type().is_protein() ) {
      if( claim_backbone() ){
        insert_dof_element( conf, elements, seqpos, BB, phi_torsion );
        insert_dof_element( conf, elements, seqpos, BB, psi_torsion );
        insert_dof_element( conf, elements, seqpos, BB, omega_torsion );
      }
      if( claim_sidechain() ){
        for( Size i = 1; i <= conf.residue( seqpos ).nchi(); ++i ){
          insert_dof_element( conf, elements, seqpos, CHI, i );
        }
      }
    } else if( conf.residue( seqpos ).is_virtual_residue() ) {
      tr.Debug << *this << " ignoring seqpos " << seqpos << ", because it's a virtual residue." << std::endl;
    } else if( conf.residue( seqpos ).is_NA() ){
      tr.Fatal << "[FATAL] The Broker is not currently built to deal with nucleic acids. Go in to TorsionClaim::"
         << __FUNCTION__ << " in " << __FILE__ << ":" << __LINE__ << " and build the appropriate DOF_IDs." << std::endl;
      utility_exit();
    } else {
      // Hey there! If you got here because you're using sugars or something,
      // all you have to do is put an else if for your case (see the virtual check above),
      // and correctly generate the DOF_IDs for all your torsions. Then call wrap_dof_id( dof_id )
      // to build a DOFElement out of it, and put it into elements. See insert_dof_element above. JRP

      std::ostringstream ss;
      ss << "[ERROR] TorsionClaim::" << __FUNCTION__ << " owned by " << owner()->get_name()
         << " doesn't know how to handle residue " << seqpos << " with name3 "
         << conf.residue( seqpos ).name3() << " because it is not protein."
         << "If you're using non-protein residue types (sugars, ligands, RNA/DNA, etc.), check out "
         << __FILE__ << ":" << __LINE__ << " to tell the TorsionClaim how to turn torsion "
         << "numbers in to DOF_IDs" << std::endl;
      throw utility::excn::EXCN_BadInput( ss.str() );
    }
  }
}

void TorsionClaim::strength( ControlStrength const& c_str, ControlStrength const& i_str ){
  if( c_str > EXCLUSIVE || c_str < DOES_NOT_CONTROL ){
    throw utility::excn::EXCN_RangeError( "Sampling ControlStrengths are limited to values between DOES_NOT_CONTROL and EXCLUSIVE" );
  } else {
    c_str_ = c_str;
  }

  if( i_str > EXCLUSIVE || i_str < DOES_NOT_CONTROL ){
    throw utility::excn::EXCN_RangeError( "Initialization ControlStrengths are limited to values between DOES_NOT_INITIALIZE and EXCLUSIVE" );
  } else {
    i_str_ = i_str;
  }
}

ControlStrength const& TorsionClaim::ctrl_strength() const {
  return c_str_;
}


ControlStrength const& TorsionClaim::init_strength() const {
  return i_str_;
}

EnvClaimOP TorsionClaim::clone() const {
  return EnvClaimOP( new TorsionClaim( *this ) );
}

std::string TorsionClaim::type() const{
  return "Torsion";
}

void TorsionClaim::show( std::ostream& os ) const {
  os << type() << " owned by a '" << owner()->get_name() << "' with control strength "
     << ctrl_strength() << " and init strength " << init_strength() << " with selector " << selector()->get_name();
}

} //claims
} //environment
} //protocols
