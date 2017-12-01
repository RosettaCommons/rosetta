// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/MembraneSpanTermZConstraint.cc
///
/// @brief constraint the Z axis distance between TM span edges
/// @details in hetero fold and dock (w/w- design) rosetta commonly shifts the chains along the Z axis
/// if you prefer them not shofting, this constraint penalizes the pose by the distance on the Z axis between the span edges
/// @author Jonathan Weinstein (jonathan.weinstein@weizmann.ac.il)


#include <core/scoring/constraints/MembraneSpanTermZConstraint.hh>

#include <core/conformation/Residue.hh>
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
#include <cereal/types/string.hpp>
#include <cereal/types/utility.hpp>

#endif // SERIALIZATION
// Membrane headers
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/Exceptions.hh>


namespace core {
namespace scoring {
namespace constraints {


static basic::Tracer TR( "core.scoring.constraints.MembraneSpanTermZConstraint" );

MembraneSpanTermZConstraint::MembraneSpanTermZConstraint():
	Constraint( core::scoring::membrane_span_term_z_constraint )
{
}

MembraneSpanTermZConstraint::MembraneSpanTermZConstraint(
	core::pose::Pose const & pose
):
	Constraint( core::scoring::membrane_span_term_z_constraint )
{
	utility::vector1< utility::pointer::shared_ptr< core::conformation::membrane::Span > > const spans = pose.conformation().membrane_info()->spanning_topology()->get_spans();

	in_res_.reserve( spans.size() );
	out_res_.reserve( spans.size() );

	for ( core::Size i = 1; i <= spans.size(); ++i ) {
		core::Size start = spans[ i ]->start();
		core::Size end = spans[ i ]->end();

		in_res_.push_back( pose.residue( start ).get_self_ptr() );
		out_res_.push_back( pose.residue( end ).get_self_ptr() );
	}
}

MembraneSpanTermZConstraint::~MembraneSpanTermZConstraint() {}

ConstraintOP
MembraneSpanTermZConstraint::clone() const
{
	return ConstraintOP( new MembraneSpanTermZConstraint( *this ) );
}

utility::vector1< core::Size >
MembraneSpanTermZConstraint::residues() const {
	utility::vector1< core::Size > pos_list( 1 );
	for ( auto const & res : in_res_ ) {
		pos_list.push_back( res->seqpos() ); // length 1 containing "all" seqpos_ values
	}
	for ( auto const & res : out_res_ ) {
		pos_list.push_back( res->seqpos() ); // length 1 containing "all" seqpos_ values
	}
	return pos_list;
}

void
MembraneSpanTermZConstraint::show( std::ostream & out ) const {
	out << "MembraneSpanTermZConstraint" << std::endl;
	out << "in facing span termini" << std::endl;
	for ( auto const & res1 : in_res_ ) {
		for ( auto const & res2 : in_res_ ) {
			if ( res1->seqpos() != res2->seqpos() ) {
				out << "in: res1 " << res1->seqpos() << " Z " << res1->nbr_atom_xyz().z() << " res2 " << res2->seqpos() << " " << res2->nbr_atom_xyz().z() << " score " << score_delta( res1->nbr_atom_xyz().z(), res2->nbr_atom_xyz().z() )<< std::endl;
			}
		}
	}
	out << "out facing span termini" << std::endl;
	for ( auto const & res1 : out_res_ ) {
		for ( auto const & res2 : out_res_ ) {
			if ( res1->seqpos() != res2->seqpos() ) {
				out << "out: res1 " << res1->seqpos() << " Z " << res1->nbr_atom_xyz().z() << " res2 " << res2->seqpos() << " " << res2->nbr_atom_xyz().z() << " score " << score_delta( res1->nbr_atom_xyz().z(), res2->nbr_atom_xyz().z() )<< std::endl;
			}
		}
	}
}

ConstraintOP
MembraneSpanTermZConstraint::remap_resid( core::id::SequenceMapping const & ) const
{
	return NULL;
}

bool
MembraneSpanTermZConstraint::operator == ( Constraint const & other_cst ) const
{
	if ( ! same_type_as_me( other_cst ) ) return false;
	if ( ! other_cst.same_type_as_me( *this ) ) return false;

	MembraneSpanTermZConstraint const & other( static_cast< MembraneSpanTermZConstraint const & > (other_cst) );

	if ( in_res_ != other.in_res_ ) return false;
	if ( out_res_ != other.out_res_ ) return false;
	if ( score_type() != other.score_type() ) return false;

	return true;
}

bool
MembraneSpanTermZConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< MembraneSpanTermZConstraint const * > ( &other );
}

ConstraintOP
MembraneSpanTermZConstraint::remapped_clone( pose::Pose const& pose, pose::Pose const&, id::SequenceMappingCOP ) const {
	return ConstraintOP( new MembraneSpanTermZConstraint( pose ) );
}


// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
MembraneSpanTermZConstraint::score( func::XYZ_Func const &, EnergyMap const &, EnergyMap & emap ) const
{
	emap[ this->score_type() ] += calc_score();
	//show( TR );
}

core::Real
MembraneSpanTermZConstraint::calc_score() const
{
	core::Real score;
	score = 0;

	for ( auto const & res1 : in_res_ ) {
		for ( auto const & res2 : in_res_ ) {
			if ( res1->seqpos() != res2->seqpos() ) {
				score += score_delta( res1->nbr_atom_xyz().z(),  res2->nbr_atom_xyz().z());
			}
		}
	}
	for ( auto const & res1 : out_res_ ) {
		for ( auto const & res2 : out_res_ ) {
			if ( res1->seqpos() != res2->seqpos() ) {
				score += score_delta( res1->nbr_atom_xyz().z(),  res2->nbr_atom_xyz().z());
			}
		}
	}
	score /= 2;
	return score;
}

core::Real
MembraneSpanTermZConstraint::score_delta( core::Real a, core::Real b ) const
{
	return std::pow( a - b, 2 );
}

core::Real
MembraneSpanTermZConstraint::dist( core::scoring::func::XYZ_Func const & ) const {
	return 0.0;
}

void
MembraneSpanTermZConstraint::fill_f1_f2(
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
core::scoring::constraints::MembraneSpanTermZConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( in_res_ ) ); // Size
	arc( CEREAL_NVP( out_res_ ) ); // Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::MembraneSpanTermZConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	utility::vector1< core::conformation::ResidueOP > local_in_res;
	arc ( local_in_res );
	in_res_ = local_in_res;
	utility::vector1< core::conformation::ResidueOP > local_out_res;
	arc( local_out_res );
	out_res_ = local_out_res;
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::MembraneSpanTermZConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::MembraneSpanTermZConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_MembraneSpanTermZConstraint )
#endif // SERIALIZATION
