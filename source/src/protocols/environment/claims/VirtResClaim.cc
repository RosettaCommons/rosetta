// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file VirtResClaim
/// @brief Claims access to a torsional angle.
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/claims/VirtResClaim.hh>

// Package Headers
#include <core/environment/LocalPosition.hh>
#include <core/environment/DofPassport.hh>

#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/BrokerElements.hh>

#include <protocols/environment/ProtectedConformation.hh>
#include <core/environment/SequenceAnnotation.hh>
#include <protocols/environment/ClaimingMover.hh>

// Project Headers
#include <core/id/TorsionID.hh>

#include <core/pose/util.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

// C++ headers

// option key includes

static basic::Tracer tr("protocols.environment.VirtResClaim",basic::t_info);

namespace protocols {
namespace environment {
namespace claims {

VirtResClaim::VirtResClaim( ClaimingMoverOP owner,
                            LocalPosition parent,
                            std::string const& jump_label,
                            std::string const& vrt_label ):
  EnvClaim( owner ),
  vrt_label_( vrt_label ),
  parent_( parent ),
  j_claim( owner, jump_label, parent, LocalPosition( vrt_label, 1 ), LocalPosition( vrt_label, 0 ) )
{}

void VirtResClaim::yield_elements( FoldTreeSketch const&, ResidueElements& elements ) const{
  ResidueElement e;

  e.label = vrt_label();

  elements.push_back( e );
}

void VirtResClaim::yield_elements( FoldTreeSketch const& fts, JumpElements& elements ) const {
  j_claim.yield_elements( fts, elements );
}

EnvClaimOP VirtResClaim::clone() const {
  return new VirtResClaim( *this );
}

std::string const& VirtResClaim::jump_label() const {
  return j_claim.label();
}

std::string const& VirtResClaim::vrt_label() const {
  return vrt_label_;
}

std::string VirtResClaim::str_type() const{
  return "VirtRes";
}

void VirtResClaim::show( std::ostream& os ) const {
  os << str_type() << " '" << vrt_label() << "with jump '"
     << jump_label() << "' owned by a " << owner()->get_name();
}

} //claims
} //environment
} //protocols
