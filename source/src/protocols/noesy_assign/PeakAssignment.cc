// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange

// Unit Headers
#include <protocols/noesy_assign/PeakAssignment.hh>
#include <protocols/noesy_assign/CrossPeak.hh>

#include <protocols/noesy_assign/ResonanceList.hh>

// Package Headers
#include <protocols/noesy_assign/PeakAssignmentParameters.hh>

// Project Headers
#include <protocols/noesy_assign/util.hh>
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>
#include <core/scoring/constraints/AmbiguousNMRConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>

// Utility headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>

// Utility headers


//// C++ headers
#include <cstdlib>
#include <string>

#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.noesy_assign.assignments" );

namespace protocols {
namespace noesy_assign {

/// @details Auto-generated virtual destructor
PeakAssignment::~PeakAssignment() {}
using namespace core;

PeakAssignment::PeakAssignment( CrossPeak * cp, core::Size assign_spin1, core::Size assign_spin2 )
: crosspeak_( cp ),
	spin_assign_index1_( assign_spin1 ),
	spin_assign_index2_( assign_spin2 ),
	chemshift_overlap_( 1.0 ),
	symmetry_compliance_( 0 ),
	covalent_compliance_( false ),
	decoy_compatibility_( 1.0 ),
	network_anchoring_( 1.0 ),
	network_anchoring_per_residue_( 200 ), //so we are not eliminated by default
	native_distance_viol_( -1 )
{
	if ( cp )  update_resonances_from_peak();
}

ResonanceList const& PeakAssignment::resonances() const {
	return crosspeak_->resonances();
}

void PeakAssignment::dump_weights( std::ostream& os ) const {
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	using namespace ObjexxFCL::format;
	os << F( 5, 2, chemshift_overlap_ ) << " "
		<< F( 5, 2, std::max( symmetry_compliance_*params.symmetry_compliance_weight_, 1.0 ) ) << " "
		<< F( 5, 2, covalent_compliance_ ? params.covalent_compliance_weight_ : 1.0 ) << " "
		<< F( 5, 2, decoy_compatibility_ ) << " "
		<< F( 5, 2, network_anchoring_ ) << " "
		<< F( 5, 2, network_anchoring_per_residue_ ) << " ";
	os << F( 5, 2, native_distance_viol_ ) << " ";
}

Real PeakAssignment::peak_volume() const {
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	return chemshift_overlap_
		* std::min( params.smax_,
		std::max( symmetry_compliance_*params.symmetry_compliance_weight_, 1.0 )
		* ( covalent_compliance_ ? params.covalent_compliance_weight_ : 1.0 )
		* network_anchoring_ )
		* decoy_compatibility_;
}

Real PeakAssignment::normalized_peak_volume() const {
	if ( crosspeak_->cumulative_peak_volume() > 0.0 ) return peak_volume() / crosspeak_->cumulative_peak_volume();
	return 0.0;
}

void PeakAssignment::update_resonances_from_peak() {
	resonance1_ = crosspeak_->proton( 1 ).assignment( spin_id( 1 ) );
	resonance2_ = crosspeak_->proton( 2 ).assignment( spin_id( 2 ) );
}

void PeakAssignment::update_chemshiftscore_from_peak() {
	PeakAssignmentParameters const& params( *PeakAssignmentParameters::get_instance() );
	Real const& weight( params.chemshift_overlap_weight_ );
	Real sum( 0.0 );
	for ( Size d=1; d<=crosspeak_->dimension(); d++ ) {
		CrossPeak::Spin const& spin( crosspeak_->spin( d ) );
		Resonance const& assigned_resonance( resonances()[ spin.assignment( spin_id( d>2 ? d-2 : d ) ) ] );
		Real s = 1.0/weight*assigned_resonance.pmatch( spin.freq(), crosspeak_->tolerance( d ), crosspeak_->folder( d ) );
		//  Real diff( spin.freq()-crosspeak_->fold_resonance( assigned_resonance.freq(), d ) );
		//  Real s( diff/weight/std::max( crosspeak_->tolerance( d ), assigned_resonance.tolerance() ) );
		sum += s*s;
	}
	chemshift_overlap_ = exp( -0.5*sum );
}

/// @details WARNING WARNING WARNING THREAD UNSAFE covalent_compliance RELIES ON THREAD-UNSAFE SINGLETON CovalentCompliance
void PeakAssignment::update_upperdistance_score(/*dmax*/ ) {
	core::id::NamedAtomID const& atom1( atom( 1 ) );
	core::id::NamedAtomID const& atom2( atom( 2 ) );
	covalent_compliance_ = covalent_compliance( atom1, atom2 );
}

core::Size PeakAssignment::float_ambiguity() const {
	Size float_ambiguity=1;
	for ( Size iatom=1; iatom<=2; ++iatom ) {
		if ( has_proton( iatom ) ) float_ambiguity *= resonance( iatom ).ambiguity();
	}
	return float_ambiguity;
}

PeakAssignment::NmrConstraintOP PeakAssignment::create_constraint(
	pose::Pose const& pose,
	core::Size ifloat,
	core::scoring::func::FuncOP func
) const {
	basic::ProfileThis doit( basic::NOESY_ASSIGN_PA_GEN_CST );

	core::Size ifloat_1( (ifloat-1)/resonance( 2 ).ambiguity() + 1 );
	core::Size ifloat_2( (ifloat-1)%resonance( 2 ).ambiguity() + 1 );
	core::Size resonance1 = resonances()[ resonance_id( 1 ) ].float_label( ifloat_1 );
	core::Size resonance2 = resonances()[ resonance_id( 2 ) ].float_label( ifloat_2 );

	bool flip = resonances()[ resonance1 ].resid() > resonances()[ resonance2 ].resid();
	core::id::NamedAtomID const& atom1( resonances()[ flip ? resonance2 : resonance1 ].atom() );
	core::id::NamedAtomID const& atom2( resonances()[ flip ? resonance1 : resonance2 ].atom() );

	tr.Debug << "create constraint for atom: " << atom1 << " " << atom2 << " from resonances id: "
		<< resonance1 << " " << resonance2 << std::endl;

	using namespace core::scoring::constraints;
	if ( !func ) {
		func = core::scoring::func::FuncOP( new BoundFunc( 1.5,
			5.5,
			1.0,
			"VC "+ObjexxFCL::string_of( normalized_peak_volume(), 3 )
			) );
	}
	return PeakAssignment::NmrConstraintOP( new NmrConstraint( atom1, atom2, pose, func ) ); //later figure out what the func should be...
}

core::scoring::constraints::ConstraintOP PeakAssignment::create_constraint(
	pose::Pose const& pose,
	core::scoring::func::FuncOP func
) const {
	if ( float_ambiguity() > 1 ) {
		return create_float_constraint( pose, func );
	}
	return create_constraint( pose, 1, func );
}

core::scoring::constraints::AmbiguousNMRConstraintOP PeakAssignment::create_float_constraint(
	pose::Pose const& pose,
	core::scoring::func::FuncOP func
) const {
	/// create Ambiguous constraint with BoundFunc describing the whole thing...
	using namespace core::scoring::constraints;
	AmbiguousNMRConstraintOP full_cst( new AmbiguousNMRConstraint( func ) );
	core::Size const n_float = float_ambiguity();

	// add individual constraints --- their potential does not matter ... it is ignored in evaluation...
	for ( core::Size ifloat = 1; ifloat <= n_float; ++ifloat ) {
		AmbiguousNMRDistanceConstraintOP new_cst( create_constraint( pose, ifloat ) );
		full_cst->add_individual_constraint( new_cst );
	}
	return full_cst;
}

/// @brief return resonance_id, i.e., pointer into Resonance list that will resolve in assigned atom
core::Size PeakAssignment::label_resonance_id( core::Size select ) const {
	runtime_assert( select == 1 || select == 2 );
	runtime_assert( crosspeak_ );
	runtime_assert( crosspeak_->has_label( select ) );
	return crosspeak_->label( select ).assignment( spin_id( select ) );
}

PeakAssignment const BOGUS_ASSIGNMENT( NULL, 0, 0 );

bool PeakAssignment::is_symmetric_partner_of( PeakAssignment const& other ) const {
	bool match=true;

	//because I also allow now the peaks that have not turned around residue numbers we have to exclude assignments with the same peak
	if ( crosspeak_->same_peak( *other.crosspeak_ ) ) {
		//  tr.Debug << " SAME " << std::endl;
		return false;
	}
	for ( Size select = 1; select <=2 && match; ++select ) {
		core::Size self = select;
		core::Size partner = select % 2 + 1;
		if ( has_label( self ) && other.has_label( partner ) ) {
			match = match && label_resonance_id( self ) == other.label_resonance_id( partner );
		}
		//  if ( has_proton( self ) && other.has_proton( partner ) ) {
		match = match && resonance_id( self ) == other.resonance_id( partner );
		//  } else {
		//   match = match && label_atom_
	}
	// return match;
	//we also have to match 1-1 and 2-2 since this can be the case when comparing
	// CHH and CCH spectra.
	if ( crosspeak_->exp_hash() == other.crosspeak().exp_hash() ) return match;
	if ( match ) return true;
	match = true;
	for ( Size select = 1; select <=2 && match; ++select ) {
		core::Size self = select;
		core::Size partner = select;
		if ( has_label( self ) && other.has_label( partner ) ) {
			match = match && label_resonance_id( self ) == other.label_resonance_id( partner );
		}
		//  if ( has_proton( self ) && other.has_proton( partner ) ) {
		match = match && resonance_id( self ) == other.resonance_id( partner );
		//  }
	}
	return match;
	// return resonance_id( 1 ) == other.resonance_id( 2 ) && resonance_id( 2 ) == other.resonance_id( 1 );
}

void PeakAssignment::show( std::ostream& os ) const {
	for ( Size select =1; select <= 2; ++select ) {
		if ( has_proton( select ) ) {
			os << atom( select ) << "   ";
		} else {
			os << "[ " << atom( select ) << "] ";
		}
		if ( has_label( select ) ) {
			os << label_atom( select ) << "   ";
		}
	}
}

std::ostream& operator<<( std::ostream& os, PeakAssignment const& pa ) {
	pa.show( os );
	return os;
}

// void PeakAssignment::invalidate_assignment() {
//   if ( crosspeak_ ) crosspeak_->invalidate_assignment( assignment_index_ );
// }

}
}
