// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/noesy_assign/PeakAssignment.hh
/// @author Oliver Lange

#ifndef INCLUDED_protocols_noesy_assign_PeakAssignment_HH
#define INCLUDED_protocols_noesy_assign_PeakAssignment_HH


// Unit Headers
#include <protocols/noesy_assign/PeakAssignment.fwd.hh>
// Package Headers

#include <protocols/noesy_assign/CrossPeak.hh>

#include <protocols/noesy_assign/ResonanceList.hh>
#include <protocols/noesy_assign/Resonance.hh>


// Project Headers
// AUTO-REMOVED #include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>
#include <core/pose/Pose.fwd.hh>
// AUTO-REMOVED #include <core/scoring/func/Func.hh>

#include <core/id/NamedAtomID.fwd.hh>
#include <core/types.hh>

#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.fwd.hh>
#include <core/scoring/constraints/AmbiguousNMRConstraint.fwd.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
// AUTO-REMOVED #include <cstdlib>
#include <string>
// AUTO-REMOVED #include <list>
#include <map>
#include <iostream>

#include <core/scoring/func/Func.hh>

//Auto Headers

namespace protocols {
namespace noesy_assign {


///@brief fast access to assignments that are stored in CrossPeak -- similar to FragID
class PeakAssignment : public utility::pointer::ReferenceCount {
public:
	///@brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~PeakAssignment();
	typedef core::scoring::constraints::AmbiguousNMRDistanceConstraintOP NmrConstraintOP;
	typedef core::scoring::constraints::AmbiguousNMRDistanceConstraint NmrConstraint;

  PeakAssignment( CrossPeak *, core::Size assign_spin1, core::Size assign_spin2 );

  core::Size spin_id( core::Size select ) const {
    runtime_assert( select == 1 || select == 2 );
    return select == 1 ? spin_assign_index1_ : spin_assign_index2_;
  }

  ///@brief return resonance_id, i.e., pointer into Resonance list that will resolve in assigned atom
  core::Size resonance_id( core::Size select ) const {
    runtime_assert( select == 1 || select == 2 );
    return select == 1 ? resonance1_ : resonance2_;
  }

	bool has_label( core::Size select ) const {
		return  crosspeak_->has_label( select );
	}

	bool has_proton( core::Size select ) const {
		return  crosspeak_->has_proton( select );
	}

	core::Size float_ambiguity() const;

	///@brief return resonance_id, i.e., pointer into Resonance list that will resolve in assigned atom
  core::Size label_resonance_id( core::Size select ) const;

 ///@brief returns atom 1 or 2 of the assigned cross-peak.  --- might throw Exception if atom not found in ResonanceList
  core::id::NamedAtomID const& atom( core::Size iatom ) const{
    return resonances()[ resonance_id( iatom ) ].atom();
  }

  core::id::NamedAtomID const& label_atom( core::Size iatom ) const{
    return resonances()[ label_resonance_id( iatom ) ].atom();
  }

	Resonance const& resonance( core::Size iatom ) const {
		return resonances()[ resonance_id( iatom ) ];
	}

	Resonance const& label_resonance( core::Size iatom ) const {
		return resonances()[ label_resonance_id( iatom ) ];
	}

	CALIBRATION_ATOM_TYPE calibration_atom_type( core::Size iatom ) const{
    return resonances()[ resonance_id( iatom ) ].calibration_atom_type();
  }

	core::scoring::constraints::ConstraintOP create_constraint(
  	core::pose::Pose const& pose,
		core::scoring::func::FuncOP = NULL  ) const;


	NmrConstraintOP create_constraint(
  	core::pose::Pose const& pose,
		core::Size ifloat, //if float ambiguity is present enumerate all possible constraints with 1<=ifloat <=float_ambiguity()
		core::scoring::func::FuncOP = NULL  ) const;

  ///@brief returns residue number of a1 or a2 of the assigned cross-peak, --- might throw Exception if atom not found
  core::Size resid( core::Size iatom ) const {
    return resonances()[ resonance_id( iatom ) ].resid();
  }

  bool operator==( PeakAssignment const& other ) const {
    return crosspeak_ == other.crosspeak_ && spin_id( 1 ) == other.spin_id( 1 ) && spin_id( 2 ) == other.spin_id( 2 );
  }

  bool is_symmetric_partner_of( PeakAssignment const& other ) const;

	ResonanceList const& resonances() const; //can't inline due to circularity
  //void invalidate_assignment();

  CrossPeak& crosspeak() { return *crosspeak_; };
  CrossPeak const& crosspeak() const { return *crosspeak_; }

  //core::Size assignment_index() const { return assignment_index_; }
	void set_symmetry( core::Real setting = 0 ) {
		symmetry_compliance_ = setting;
	}

	core::Real symmetry_compliance() const {
		return symmetry_compliance_;
	}

	core::Real chemshift_compliance() const {
		return chemshift_overlap_; //Ck
	}

	void set_network_anchoring( core::Real setting ) {
		network_anchoring_ = setting;
		//		network_anchoring_per_residue_ = reswise_setting;
	}

	void set_network_anchoring_per_residue( core::Real reswise_setting ) {
		network_anchoring_per_residue_ = reswise_setting;
	}

	core::Real network_anchoring() const {
		return network_anchoring_;
	}

	core::Real network_anchoring_per_residue() const {
		return network_anchoring_per_residue_;
	}

	core::Real decoy_compatibility() const {
		return decoy_compatibility_;
	}

	void set_decoy_compatibility( core::Real setting ) {
		decoy_compatibility_ = setting;
	}

	void update_chemshiftscore_from_peak();
	void update_upperdistance_score();

	///@brief this is not normalized
	core::Real peak_volume() const;

	///@brief only correct if peak_volumes have been update in CrossPeaks.
	core::Real normalized_peak_volume() const;

	void show( std::ostream& os ) const;

	void dump_weights( std::ostream& os ) const;

	void set_native_distance_viol( core::Real setting ) {
		native_distance_viol_= setting;
	}
	core::Real native_distance_viol() {
		return native_distance_viol_;
	}

private:
	void update_resonances_from_peak();

	core::scoring::constraints::AmbiguousNMRConstraintOP create_float_constraint(
  	core::pose::Pose const& pose,
		core::scoring::func::FuncOP = NULL  ) const;

  CrossPeak * crosspeak_;
  core::Size spin_assign_index1_; //points to assignment of spin1 and spin2
  core::Size spin_assign_index2_; //points to assignment of spin1 and spin2
  core::Size resonance1_; //the spin1 and spin2 resonance (short cut to avoid lengthy access thru CrossPeak)
  core::Size resonance2_; //the spin1 and spin2 resonance (short cut to avoid lengthy access thru CrossPeak)


  // should these live in inherited class or in PeakAssignmentWeights ?
  core::Real chemshift_overlap_; //Ck
	core::Real symmetry_compliance_; //Tk
  bool covalent_compliance_; //Ok
  core::Real decoy_compatibility_; //Dk
  core::Real network_anchoring_; //Nk
	core::Real network_anchoring_per_residue_; //N_{AB}
	core::Real native_distance_viol_;

};

std::ostream& operator<<( std::ostream& os, PeakAssignment const& pa );

extern PeakAssignment const BOGUS_ASSIGNMENT;

}
}
#endif
