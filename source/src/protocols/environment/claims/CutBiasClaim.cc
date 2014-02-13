// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CutBiasClaim
/// @brief Claims access to a torsional angle.
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/claims/CutBiasClaim.hh>

// Package Headers
#include <core/environment/LocalPosition.hh>
#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/ClaimingMover.hh>

// Project Headers
#include <core/fragment/SecondaryStructure.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

// C++ headers

// option key includes


static basic::Tracer tr("protocols.environment.CutBiasClaim",basic::t_info);

namespace protocols {
namespace environment {
namespace claims {

using core::environment::LocalPosition;
using core::environment::LocalPositions;

CutBiasClaim::CutBiasClaim( ClaimingMoverOP owner,
                            std::string const& label,
                            core::fragment::SecondaryStructure const& ss_in ):
  Parent( owner ),
  label_( label ),
  ss_( ss_in )
{}

void CutBiasClaim::yield_elements( FoldTreeSketch const&, CutBiasElements& elements ) const {

  ObjexxFCL::FArray1D_float const& loop_frac = ss_.loop_fraction();

  for( Size i = 1; i <= ss_.total_residue(); ++i ){
    CutBiasElement e;

    e.p = LocalPosition( label_, i );
    e.bias = loop_frac( (int) i );

    elements.push_back( e );
  }

}

EnvClaimOP CutBiasClaim::clone() const {
  return new CutBiasClaim( *this );
}

std::string CutBiasClaim::str_type() const{
  return "CutBias";
}

void CutBiasClaim::show( std::ostream& os ) const {
  os << str_type() << " owned by a " << owner()->get_name() << " with ss " << ss_;
}

} //claims
} //environment
} //protocols
