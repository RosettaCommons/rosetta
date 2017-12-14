// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/NonResidueTypeConstraint.cc
///
/// @brief
/// @author Sarel Fleishman


#include <core/scoring/constraints/NonResidueTypeConstraint.hh>

#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>

#include <core/id/SequenceMapping.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility>
#include <utility/vector1.hh>


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


static basic::Tracer non_residue_type_constraint_tracer( "core.scoring.constraints.NonResidueTypeConstraint" );


NonResidueTypeConstraint::NonResidueTypeConstraint(
	core::pose::Pose const & pose,
	Size seqpos,
	core::Real favor_non_native_bonus
):
	Constraint( core::scoring::res_type_constraint ),
	seqpos_( seqpos ),
	rsd_type_name3_( pose.residue_type( seqpos ).name3() ),
	favor_non_native_bonus_( favor_non_native_bonus )
{}


NonResidueTypeConstraint::NonResidueTypeConstraint(
	core::pose::Pose const &,
	Size seqpos,
	std::string const & AAname,
	core::Real favor_non_native_bonus
):
	Constraint( core::scoring::res_type_constraint ),
	seqpos_( seqpos ),
	rsd_type_name3_( AAname ),
	favor_non_native_bonus_( favor_non_native_bonus )
{}

NonResidueTypeConstraint::NonResidueTypeConstraint(
	Size seqpos,
	std::string const & aa_in,
	std::string const & name3_in,
	core::Real bonus_in
):
	Constraint( core::scoring::res_type_constraint ),
	seqpos_( seqpos ),
	AAname( aa_in ),
	rsd_type_name3_( name3_in ),
	favor_non_native_bonus_( bonus_in )
{}


NonResidueTypeConstraint::~NonResidueTypeConstraint() = default;

ConstraintOP
NonResidueTypeConstraint::clone() const
{
	return ConstraintOP( new NonResidueTypeConstraint( *this ) );
}

utility::vector1< core::Size >
NonResidueTypeConstraint::residues() const {
	utility::vector1< core::Size > pos_list(1, seqpos_); // length 1 containing "all" seqpos_ values
	return pos_list;
}

void
NonResidueTypeConstraint::show( std::ostream & out ) const {
	out << "NonResidueTypeConstraint; ";
	out << "seqpos: " << seqpos_;
	out << "; AAname: "<< AAname;
	out << "; rsd_type_name3: "<< rsd_type_name3_;
	out << "; favor_non_native_bonus: "<< favor_non_native_bonus_;
}

ConstraintOP
NonResidueTypeConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	core::Size newseqpos = seqmap[ seqpos_ ];
	if ( newseqpos != 0 ) {
		return ConstraintOP( new NonResidueTypeConstraint( newseqpos, AAname, rsd_type_name3_, favor_non_native_bonus_ ) );
	} else {
		return nullptr;
	}
}


// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
NonResidueTypeConstraint::score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const
{
	Real const weight(weights[ this->score_type() ] );
	if ( weight == 0 ) return; // what's the point?

	conformation::Residue const & rsd( xyz_func.residue(seqpos_) );
	if ( rsd.type().name3() == rsd_type_name3_ ) {
		emap[ this->score_type() ] -= favor_non_native_bonus_;
	}
	// no match, don't adjust score
}

core::Real
NonResidueTypeConstraint::dist( core::scoring::func::XYZ_Func const & xyz ) const {
	conformation::Residue const & rsd( xyz.residue(seqpos_) );
	return rsd.type().name3() == rsd_type_name3_;
}

void
NonResidueTypeConstraint::fill_f1_f2(
	AtomID const & ,//atom,
	func::XYZ_Func const &,
	Vector & ,//F1,
	Vector & ,//F2,
	EnergyMap const & //weights
) const
{
	// Do nothing.
	// Derivative of this restraint is effectively zero
	// so we just "add zero" to F1 and F2.
}

bool NonResidueTypeConstraint::operator == ( Constraint const & rhs ) const
{
	if ( ! same_type_as_me( rhs ) ) return false;
	if ( ! rhs.same_type_as_me( *this ) ) return false;

	auto const & rhs_nrtc( static_cast< NonResidueTypeConstraint const & > ( rhs ) );
	if ( seqpos_ != rhs_nrtc.seqpos_ ) return false;
	if ( AAname != rhs_nrtc.AAname ) return false;
	if ( rsd_type_name3_ != rhs_nrtc.rsd_type_name3_ ) return false;
	return favor_non_native_bonus_ == rhs_nrtc.favor_non_native_bonus_;
}

bool NonResidueTypeConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< NonResidueTypeConstraint const *  > ( & other );
}


} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::constraints::NonResidueTypeConstraint::NonResidueTypeConstraint() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::NonResidueTypeConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( seqpos_ ) ); // Size
	arc( CEREAL_NVP( AAname ) ); // std::string
	arc( CEREAL_NVP( rsd_type_name3_ ) ); // std::string
	arc( CEREAL_NVP( favor_non_native_bonus_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::NonResidueTypeConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	arc( seqpos_ ); // Size
	arc( AAname ); // std::string
	arc( rsd_type_name3_ ); // std::string
	arc( favor_non_native_bonus_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::NonResidueTypeConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::NonResidueTypeConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_NonResidueTypeConstraint )
#endif // SERIALIZATION
