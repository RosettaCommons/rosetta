// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file CrossPeakList.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_Assignment_hh
#define INCLUDED_protocols_noesy_assign_Assignment_hh


// Unit Headers
#include <protocols/noesy_assign/PeakAssignment.fwd.hh>
#include <protocols/noesy_assign/CrossPeakList.hh>
#include <protocols/noesy_assign/ResonanceList.hh>
#include <core/types.hh>

// Package Headers
#include <core/id/NamedAtomID.fwd.hh>

// Project Headers

// Utility headers
#include <utility/exit.hh>
// #include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
// #include <numeric/numeric.functions.hh>
// #include <basic/prof.hh>
//#include <basic/Tracer.hh>
// #include <basic/options/option.hh>
// #include <basic/options/keys/abinitio.OptionKeys.gen.hh>
// #include <basic/options/keys/run.OptionKeys.gen.hh>
//#include <basic/options/keys/templates.OptionKeys.gen.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <list>
#include <map>

namespace protocols {
namespace noesy_assign {


///@brief fast access to assignments that are stored in CrossPeak -- similar to FragID
class PeakAssignment : public utility::pointer::ReferenceCount  {
public:
  PeakAssignment( CrossPeakAP, core::Size assign_spin1, core::Size assign_spin2 );

  core::Size spin_id( core::Size select ) {
    runtime_assert( select == 1 || select == 2 );
    return select == 1 ? spin_assign_index1_ : spin_assign_index2_;
  }

  ///@brief return resonance_id, i.e., pointer into Resonance list that will resolve in assigned atom
  core::Size resonance_id( core::Size select ) const {
    runtime_assert( select == 1 || select == 2 );
    return select == 1 ? resonance1_ : resonance2_;
  }



 ///@brief returns atom 1 or 2 of the assigned cross-peak.  --- might throw Exception if atom not found in ResonanceList
  core::id::NamedAtomID const& atom( ResonanceList const& resonances, core::Size iatom ) const{
    return resonances[ resonance_id( iatom ) ].atom();
  }

  ///@brief returns residue number of a1 or a2 of the assigned cross-peak, --- might throw Exception if atom not found
  core::Size resid( ResonanceList const& resonances, core::Size iatom ) const {
    return resonances[ resonance_id( iatom ) ].resid();
  }


  bool operator==( PeakAssignment const& other ) const {
    return crosspeak_.get() == other.crosspeak_.get() && assignment_index_ == other.assignment_index_;
  }

  bool is_symmetric_partner_of( PeakAssignment const& other ) const {
    return resonance_id( 1 ) == other.resonance_id( 2 ) && resonance_id( 2 ) == other.resonance_id( 1 );
  }

  void invalidate_assignment();

  CrossPeak& crosspeak() { return *crosspeak_; };
  CrossPeak const& crosspeak() const { return *crosspeak_; }

  core::Size assignment_index() const { return assignment_index_; }

private:
  CrossPeakAP crosspeak_;
  core::Size spin_assign_index1_; //points to assignment of spin1 and spin2
  core::Size spin_assign_index2_; //points to assignment of spin1 and spin2
  core::Size resonance1_; //the spin1 and spin2 resonance (short cut to avoid lengthy access thru CrossPeak)
  core::Size resonance2_; //the spin1 and spin2 resonance (short cut to avoid lengthy access thru CrossPeak)



  // should these live in inherited class or in PeakAssignmentWeights ?
  Real chem_shift_overlap_; //Ck
  Real symmetry_compliance_; //Tk
  Real covalent_compliance_; //Ok
  Real decoy_compatibility_; //Dk
  Real Network_anchoring_; //Nk
};

extern PeakAssignment const BOGUS_ASSIGNMENT;


} //namespace noesy_assign
} //namespace protocols

#endif
