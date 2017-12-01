// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/MembraneSpanConstraint.cc
///
/// @brief constarint the TM spans specified by user to stay in the membrane
/// @details penalizes the pose for TM span residues that are far from the membrane mid-plane (Z=0)
/// mostly useful for hetero fold and dock, where sampling space is large, and Rosetta seems to like extracting the TM spans from the membrane
//
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)


#include <core/scoring/constraints/MembraneSpanConstraint.hh>

#include <core/conformation/Residue.hh>
//  -- REALLY?
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>

#include <core/id/SequenceMapping.hh>
#include <core/pose/util.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/Pose.hh>

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#include <cereal/access.hpp>
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>

#endif // SERIALIZATION
// Membrane headers
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/Exceptions.hh>


namespace core {
namespace scoring {
namespace constraints {


static basic::Tracer TR( "core.scoring.constraints.MembraneSpanConstraint" );

MembraneSpanConstraint::MembraneSpanConstraint():
	Constraint( core::scoring::membrane_span_constraint )
{
}

MembraneSpanConstraint::MembraneSpanConstraint(
	core::pose::Pose const & pose
):
	Constraint( core::scoring::membrane_span_constraint )
{
	utility::vector1< utility::pointer::shared_ptr< core::conformation::membrane::Span > > const spans = pose.conformation().membrane_info()->spanning_topology()->get_spans();

	span_res_.reserve( spans.size() );

	for ( core::Size i = 1; i <= spans.size(); ++i ) {
		core::Size res_i = spans[ i ]->start() + ( spans[ i ]->end() - spans[ i ]->start() ) / 2;
		span_res_.push_back( pose.residue( res_i ).get_self_ptr() );
	}
}

MembraneSpanConstraint::~MembraneSpanConstraint() {}

ConstraintOP
MembraneSpanConstraint::clone() const
{
	return ConstraintOP( new MembraneSpanConstraint( *this ) );
}

utility::vector1< core::Size >
MembraneSpanConstraint::residues() const {
	utility::vector1< core::Size > pos_list( 1 );
	for ( auto const & res : span_res_ ) {
		pos_list.push_back( res->seqpos() ); // length 1 containing "all" seqpos_ values
	}
	return pos_list;
}

void
MembraneSpanConstraint::show( std::ostream & out ) const {
	out << "MembraneSpanConstraint" << std::endl;
	out << "i    z    score" << std::endl;
	for ( auto const & res : span_res_ ) {
		out << res->seqpos() << " " << res->nbr_atom_xyz().z() << " " << score_res_z( res->nbr_atom_xyz().z() ) << std::endl;
	}
}

ConstraintOP
MembraneSpanConstraint::remap_resid( core::id::SequenceMapping const & ) const
{
	return NULL;
}

bool
MembraneSpanConstraint::operator == ( Constraint const & other_cst ) const
{
	if ( ! same_type_as_me( other_cst ) ) return false;
	if ( ! other_cst.same_type_as_me( *this ) ) return false;

	MembraneSpanConstraint const & other( static_cast< MembraneSpanConstraint const & > (other_cst) );

	if ( span_res_ != other.span_res_ ) return false;
	if ( score_type() != other.score_type() ) return false;

	return true;
}

bool
MembraneSpanConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< MembraneSpanConstraint const * > ( &other );
}

ConstraintOP
MembraneSpanConstraint::remapped_clone( pose::Pose const& pose, pose::Pose const&, id::SequenceMappingCOP ) const {
	return ConstraintOP( new MembraneSpanConstraint( pose ) );
}


// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
MembraneSpanConstraint::score( func::XYZ_Func const &, EnergyMap const &, EnergyMap & emap ) const
{
	emap[ this->score_type() ] += calc_score();
	//show( TR );
}

core::Real
MembraneSpanConstraint::calc_score() const
{
	core::Real score( 0 );
	for ( auto const & res : span_res_ ) {
		score += score_res_z(res->nbr_atom_xyz().z());
	}
	return score;
}

core::Real
MembraneSpanConstraint::score_res_z( core::Real z ) const
{
	return std::pow( z / 8, 4 );
}

core::Real
MembraneSpanConstraint::dist( core::scoring::func::XYZ_Func const & ) const {
	return 0.0;
}

void
MembraneSpanConstraint::fill_f1_f2(
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

} // namespace constraints
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::MembraneSpanConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( span_res_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::MembraneSpanConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	utility::vector1< core::conformation::ResidueOP > local_span_res;
	arc( local_span_res );
	span_res_ = local_span_res;
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::MembraneSpanConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::MembraneSpanConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_MembraneSpanConstraint )
#endif // SERIALIZATION
