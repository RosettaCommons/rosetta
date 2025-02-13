// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/ResidueTypeConstraint.cc
///
/// @brief
/// @author Sarel Fleishman


#include <core/scoring/constraints/ResidueTypeConstraint.hh>

#include <core/conformation/Residue.hh>
//  -- REALLY?
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>

#include <core/id/SequenceMapping.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/Pose.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {


static basic::Tracer TR( "core.scoring.constraints.ResidueTypeConstraint" );

ResidueTypeConstraint::ResidueTypeConstraint():
	Constraint( core::scoring::res_type_constraint )
{}

ResidueTypeConstraint::ResidueTypeConstraint(
	core::pose::Pose const & pose,
	Size seqpos,
	core::Real favor_native_bonus
):
	Constraint( core::scoring::res_type_constraint ),
	seqpos_( seqpos ),
	rsd_type_name3_( pose.residue_type( seqpos ).name3() ),
	favor_native_bonus_( favor_native_bonus )
{}


ResidueTypeConstraint::ResidueTypeConstraint(
	core::pose::Pose const &, //pose,
	Size seqpos,
	std::string const & AA_name,
	core::Real favor_native_bonus
):
	Constraint( core::scoring::res_type_constraint ),
	seqpos_( seqpos ),
	rsd_type_name3_( AA_name ),
	favor_native_bonus_( favor_native_bonus )
{}

ResidueTypeConstraint::ResidueTypeConstraint(
	core::pose::Pose const &, //pose,
	Size seqpos,
	std::string const & AA_name,
	core::Real favor_native_bonus,
	bool base_name
):
	Constraint( core::scoring::res_type_constraint ),
	seqpos_( seqpos ),
	rsd_type_name3_( AA_name ),
	favor_native_bonus_( favor_native_bonus ),
	base_name_( base_name )
{}

ResidueTypeConstraint::ResidueTypeConstraint(
	Size seqpos,
	std::string const & aa_in,
	std::string const & name3_in,
	core::Real bonus_in
):
	Constraint( core::scoring::res_type_constraint ),
	seqpos_( seqpos ),
	AAname(std::move( aa_in )),
	rsd_type_name3_(std::move( name3_in )),
	favor_native_bonus_( bonus_in )
{}


ResidueTypeConstraint::~ResidueTypeConstraint() = default;

ConstraintOP
ResidueTypeConstraint::clone() const
{
	return ConstraintOP( new ResidueTypeConstraint( *this ) );
}

utility::vector1< core::Size >
ResidueTypeConstraint::residues() const {
	utility::vector1< core::Size > pos_list(1, seqpos_); // length 1 containing "all" seqpos_ values
	return pos_list;
}

void
ResidueTypeConstraint::show( std::ostream & out ) const {
	out << "ResidueTypeConstraint; ";
	out << "seqpos: " << seqpos_;
	out << "; AAname: "<< AAname;
	out << "; rsd_type_name3: "<< rsd_type_name3_;
	out << "; favor_native_bonus: "<< favor_native_bonus_;
}

ConstraintOP
ResidueTypeConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	core::Size newseqpos = seqmap[ seqpos_ ];
	if ( newseqpos != 0 ) {
		return utility::pointer::make_shared< ResidueTypeConstraint >( newseqpos, AAname, rsd_type_name3_, favor_native_bonus_ );
	} else {
		return nullptr;
	}
}

bool
ResidueTypeConstraint::operator == ( Constraint const & other_cst ) const
{
	if ( ! same_type_as_me( other_cst ) ) return false;
	if ( ! other_cst.same_type_as_me( *this ) ) return false;

	auto const & other( static_cast< ResidueTypeConstraint const & > (other_cst) );

	if ( seqpos_ != other.seqpos_ ) return false;
	if ( AAname != other.AAname ) return false;
	if ( rsd_type_name3_ != other.rsd_type_name3_ ) return false;
	if ( favor_native_bonus_ != other.favor_native_bonus_ ) return false;
	if ( score_type() != other.score_type() ) return false;

	return true;
}

bool
ResidueTypeConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< ResidueTypeConstraint const * > ( &other );
}

ConstraintOP
ResidueTypeConstraint::remapped_clone( pose::Pose const&, pose::Pose const&, id::SequenceMappingCOP smap ) const {

	core::Size newseqpos = seqpos_;
	if ( smap ) {
		newseqpos = (*smap)[ seqpos_ ];
		if ( newseqpos == 0 ) return nullptr;
	}

	return utility::pointer::make_shared< ResidueTypeConstraint >(newseqpos, AAname, rsd_type_name3_, favor_native_bonus_);
}

bool
ResidueTypeConstraint::get_base_name_active( ) const
{
	return base_name_;
}

void
ResidueTypeConstraint::set_base_name_active( bool base_name )
{
	base_name_ = base_name;
}

// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
ResidueTypeConstraint::score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const
{
	Real const weight(weights[ this->score_type() ] );
	if ( weight == 0 ) return; // what's the point?

	conformation::Residue const & rsd( xyz_func.residue(seqpos_) );
	if ( base_name_ ) {
		if ( rsd.type().base_name() == rsd_type_name3_ ) {
			emap[ this->score_type() ] -= favor_native_bonus_;
		}
	} else {
		if ( rsd_type_name3_.length() > 3 ) {
			TR.Warning <<
				"Warning: base_name_ is not set (set_base_name_active(true)), but rsd_type_name3_ is greater than 3: "
				<< rsd_type_name3_ << std::endl;
		}
		if ( rsd.type().name3() == rsd_type_name3_ ) {
			emap[ this->score_type() ] -= favor_native_bonus_;
		}
	}
	// no match, don't adjust score
}

core::Real
ResidueTypeConstraint::dist( core::scoring::func::XYZ_Func const & xyz ) const {
	conformation::Residue const & rsd( xyz.residue(seqpos_) );
	if ( base_name_ ) {
		return rsd.type().base_name() == rsd_type_name3_;
	} else {
		if ( rsd_type_name3_.length() > 3 ) {
			TR.Warning <<
				"Warning: base_name_ is not set (set_base_name_active(true)), but rsd_type_name3_ is greater than 3: "
				<< rsd_type_name3_ << std::endl;
		}
		return rsd.type().name3() == rsd_type_name3_;
	}
}

void
ResidueTypeConstraint::fill_f1_f2(
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

core::Real
ResidueTypeConstraint::get_favor_native_bonus() const{
	return favor_native_bonus_;
}

std::string
ResidueTypeConstraint::get_rsd_type_name3() const{
	return rsd_type_name3_;
}




} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::ResidueTypeConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( seqpos_ ) ); // Size
	arc( CEREAL_NVP( AAname ) ); // std::string
	arc( CEREAL_NVP( rsd_type_name3_ ) ); // std::string
	arc( CEREAL_NVP( favor_native_bonus_ ) ); // core::Real
	arc( CEREAL_NVP( base_name_ ) ); // bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::ResidueTypeConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	arc( seqpos_ ); // Size
	arc( AAname ); // std::string
	arc( rsd_type_name3_ ); // std::string
	arc( favor_native_bonus_ ); // core::Real
	arc( base_name_  ); // bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::ResidueTypeConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::ResidueTypeConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_ResidueTypeConstraint )
#endif // SERIALIZATION
