// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/claims/EnvClaim.cc
/// @brief  Abstract class for implementing an EnvClaim (e.g. TorsionClaim).
/// @author Justin R. Porter

// Unit Headers
#include <protocols/environment/claims/EnvClaim.hh>

// Package Headers
#include <protocols/environment/ClaimingMover.hh>
#include <core/environment/SequenceAnnotation.hh>

#include <protocols/environment/claims/ClaimStrength.hh>

#include <protocols/environment/claims/CutBiasClaim.hh>
#include <protocols/environment/claims/JumpClaim.hh>
#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/claims/VirtResClaim.hh>
#include <protocols/environment/claims/XYZClaim.hh>

// Project Headers

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

#include <utility/tag/Tag.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

//// C++ headers
#include <iostream>

// option key includes


static basic::Tracer tr("protocols.environment.EnvClaim",basic::t_info);

namespace protocols {
namespace environment {
namespace claims {

EnvClaimOP EnvClaim::make_claim( std::string const& name,
                                 ClaimingMoverOP owner,
                                 utility::tag::TagCOP tag,
                                 basic::datacache::DataMap& datamap ) {
  if      ( name == "CutBiasClaim" ) return new CutBiasClaim( owner, tag, datamap );
  else if ( name == "JumpClaim" )    return new JumpClaim( owner, tag, datamap );
  else if ( name == "TorsionClaim" ) return new TorsionClaim( owner, tag, datamap );
  else if ( name == "VrtResClaim" )  return new VirtResClaim( owner, tag, datamap );
  else if ( name == "XYZClaim" )     return new XYZClaim( owner, tag, datamap );
  else throw utility::excn::EXCN_RosettaScriptsOption( "'" + name + "' is not a known EnvClaim type." );
}

bool EnvClaim::is_claim( std::string const& name ) {
  if      ( name == "CutBiasClaim" ) return true;
  else if ( name == "JumpClaim" )    return true;
  else if ( name == "TorsionClaim" ) return true;
  else if ( name == "VrtResClaim" )  return true;
  else if ( name == "XYZClaim" )     return true;
  else return false;
}

EnvClaim::EnvClaim( ClaimingMoverOP owner ):
  ReferenceCount(),
  claim_source_( owner )
{}

/// @details Auto-generated virtual destructor
EnvClaim::~EnvClaim() {}

void EnvClaim::show( std::ostream& os ) const {
    os << "owned by, " << owner()->type() << ";";
}

DOFElement EnvClaim::wrap_dof_id( core::id::DOF_ID const& id ) const {
  DOFElement e;
  e.id = id;

  return e;
}

ClaimingMoverOP EnvClaim::owner() const {
  return claim_source_;
}

ControlStrength EnvClaim::parse_ctrl_str( std::string const& str ) const {
  std::string lower = str;
  boost::algorithm::to_lower( lower );

  if( lower == "does_not_control" ){
    return DOES_NOT_CONTROL;
  } else if( lower == "can_control" ){
    return CAN_CONTROL;
  } else if( lower == "must_control" ){
    return MUST_CONTROL;
  } else if( lower == "exclusive" ){
    return EXCLUSIVE;
  } else {
    throw utility::excn::EXCN_BadInput( "The initialization strength '" + str +
                                        "' is not recognized." );
  }
}

void EnvClaim::annotate( core::pose::Pose const& pose, core::environment::SequenceAnnotationOP ann ) const {
  BOOST_FOREACH( AnnotatingSelectors::value_type pair, selector_list_ ){
    std::string const& label = pair.first;
    ResidueSelectorCOP selector = pair.second;

    utility::vector1< bool > subset( pose.total_residue(), false );
    selector->apply( pose, subset );

    utility::vector1< Size > trues;
    for( Size i = 1; i <= subset.size(); ++i ){
      if( subset[i] )
        trues.push_back( i );
    }

    try{
      ann->add_seq_label( label, trues );
    } catch ( utility::excn::EXCN_KeyError& e ) {
      std::ostringstream ss;
      ss << "While " << *this << " was annotating the pose for broking, the " << selector->get_name()
         << "Selector produced a conflicting residue selection.";
      e.add_msg( ss.str() );
      throw e;
    }
  }
}

void EnvClaim::queue_for_annotation( std::string const& label, ResidueSelectorCOP selector ) {
  if( selector_list_.find( label ) == selector_list_.end() ||            // no selector with label 'label' exists
      selector_list_.find( label )->second.get() == selector.get() ){    // selector with label 'label' *is* 'selector'
    selector_list_[ label ] = selector;
  } else {
    std::ostringstream ss;
    ss << "In claim " << *this << ", the label '" << label << "' was to be used for the selection of a "
       << selector->get_name() << "Selector, but that label was already used by an existing "
       << selector_list_[label]->get_name() << "Selector.";
    if( selector->get_name() == selector_list_[label]->get_name() )
      ss << " They aren't the same ResidueSelector.";  // Make sure error is clear that selector-sameness has be checked already.
    throw utility::excn::EXCN_BadInput( ss.str() );
  }
}

extern std::ostream& operator<<( std::ostream& os, EnvClaim const& claim ) {
  claim.show( os );
  return os;
}

extern std::ostream& operator<<( std::ostream& os, EnvClaims const& claims ) {
  for ( EnvClaims::const_iterator it = claims.begin(); it != claims.end(); ++it ) {
    if ( *it ) {
      os << **it << "\n";
    } else {
      os << "No-Claim\n";
    }
  }
  return os;
}

} //claims
} //environment
} //protocols
