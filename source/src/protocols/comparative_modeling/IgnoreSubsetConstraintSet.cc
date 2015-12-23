// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file IgnoreSubsetConstraintSet.cc
/// @brief
/// @details
/// @author James Thompson


// Package Headers

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/types.hh>
#include <protocols/comparative_modeling/IgnoreSubsetConstraintSet.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <basic/Tracer.hh>

//// C++ headers
#include <string>
#include <set>

#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "protocols.comparative_modeling" );

using core::scoring::constraints::ConstraintSet;
using core::scoring::constraints::ConstraintSetOP;

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace comparative_modeling {

using namespace core;

IgnoreSubsetConstraintSet::IgnoreSubsetConstraintSet(
	std::set< int > residues_to_ignore,
	ConstraintSet const & other
) :
	ConstraintSet( other ),
	ignore_list_( residues_to_ignore )
{}

/// @copy constructor. Does nothing.
IgnoreSubsetConstraintSet::IgnoreSubsetConstraintSet(
	IgnoreSubsetConstraintSet const & other
) :
	ConstraintSet( other ),
	ignore_list_( other.ignore_list() )
{}

ConstraintSet &
IgnoreSubsetConstraintSet::operator = ( ConstraintSet const & rhs )
{
	if ( this != & rhs ) {
		IgnoreSubsetConstraintSet const * iscs_rhs = dynamic_cast< IgnoreSubsetConstraintSet const * > ( & rhs );
		if ( ! iscs_rhs ) {
			throw utility::excn::EXCN_Msg_Exception( "IgnoreSubsetConstraintSet handed a non IgnoreSubsetConstraintSet in operator =" );
		}
		ConstraintSet::operator = ( rhs );
		ignore_list_ = iscs_rhs->ignore_list_;
	}
	return *this;
}

ConstraintSetOP
IgnoreSubsetConstraintSet::clone() const {
	return ConstraintSetOP( new IgnoreSubsetConstraintSet( *this ) );
}

void IgnoreSubsetConstraintSet::detached_copy( ConstraintSet const & src ) {
	IgnoreSubsetConstraintSet const * iscs_src = dynamic_cast< IgnoreSubsetConstraintSet const * > ( & src );
	if ( ! iscs_src ) {
		throw utility::excn::EXCN_Msg_Exception( "IgnoreSubsetConstraintSet handed a non IgnoreSubsetConstraintSet in detatched_copy" );
	}
	deep_copy( src );
	ignore_list_ = iscs_src->ignore_list_;
}


ConstraintSetOP
IgnoreSubsetConstraintSet::detached_clone() const
{
	IgnoreSubsetConstraintSetOP newset( new IgnoreSubsetConstraintSet(*this) );
	newset->detached_copy( *this );
	return newset;
}

bool
IgnoreSubsetConstraintSet::same_type_as_me( ConstraintSet const & other, bool recurse /* = true */ ) const
{
	IgnoreSubsetConstraintSet const * iscs_other = dynamic_cast< IgnoreSubsetConstraintSet const * > ( &other );
	if ( ! iscs_other ) return false;
	if ( recurse ) {
		return other.same_type_as_me( *this, false );
	}
	return true;
}


/*void
IgnoreSubsetConstraintSet::eval_atom_derivative_for_residue_pairs (
id::AtomID const & atom_id,
pose::Pose const & pose,
scoring::ScoreFunction const &,
scoring::EnergyMap const & weights,
Vector & F1,
Vector & F2
) const
{
using scoring::constraints::ResidueConstraints;

Size const resi( atom_id.rsd() );

for ( ResidueConstraints::const_iterator
it= residue_pair_constraints_begin( resi ), ite = residue_pair_constraints_end( resi );
it != ite; ++it ) {
Size const resj( it->first );
if ( !ignore(resi) && !ignore(resj) ) {
it->second->eval_atom_derivative(
atom_id, pose.conformation(), weights, F1, F2
);
} // if ( !ignore)
}

} // eval_atom_derivative_for_residue_pairs
*/


void
IgnoreSubsetConstraintSet::residue_pair_energy(
	Residue const & rsd1,
	Residue const & rsd2,
	Pose const & pose,
	core::scoring::ScoreFunction const & scorefxn,
	core::scoring::EnergyMap & emap
) const
{
	int const pos1( rsd1.seqpos() ), pos2( rsd2.seqpos() );
	if ( ignore(pos1) || ignore(pos2) ) {
		return; //cast avoids warning
	}
	ConstraintSet::residue_pair_energy( rsd1, rsd2, pose, scorefxn, emap );
}

void
IgnoreSubsetConstraintSet::setup_for_minimizing_for_residue(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::kinematics::MinimizerMapBase const & minmap,
	core::scoring::ResSingleMinimizationData & res_data_cache
) const
{
	if ( ignore( rsd.seqpos() ) ) return;
	ConstraintSet::setup_for_minimizing_for_residue( rsd, pose, sfxn, minmap, res_data_cache );
}


void
IgnoreSubsetConstraintSet::setup_for_minimizing_for_residue_pair(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	core::scoring::ScoreFunction const & sfxn,
	core::kinematics::MinimizerMapBase const & minmap,
	core::scoring::ResSingleMinimizationData const & res1_data_cache,
	core::scoring::ResSingleMinimizationData const & res2_data_cache,
	core::scoring::ResPairMinimizationData & respair_data_cache
) const
{
	if ( ignore( rsd1.seqpos() ) || ignore( rsd2.seqpos() ) ) return;
	ConstraintSet::setup_for_minimizing_for_residue_pair(
		rsd1, rsd2, pose, sfxn, minmap, res1_data_cache, res2_data_cache, respair_data_cache );
}

bool
IgnoreSubsetConstraintSet::ignore( int const pos ) const {
	if ( ignore_list_.find(pos) == ignore_list_.end() ) {
		return false;
	}

	return true;
}

void
IgnoreSubsetConstraintSet::ignore_residue( int const pos ) {
	ignore_list_.insert( pos );
}

std::set< int >
IgnoreSubsetConstraintSet::ignore_list() const {
	return ignore_list_;
}

} // comparative_modeling
} // protocols


#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
protocols::comparative_modeling::IgnoreSubsetConstraintSet::IgnoreSubsetConstraintSet() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::comparative_modeling::IgnoreSubsetConstraintSet::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::constraints::ConstraintSet >( this ) );
	arc( CEREAL_NVP( ignore_list_ ) ); // std::set<int>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::comparative_modeling::IgnoreSubsetConstraintSet::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::constraints::ConstraintSet >( this ) );
	arc( ignore_list_ ); // std::set<int>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::comparative_modeling::IgnoreSubsetConstraintSet );
CEREAL_REGISTER_TYPE( protocols::comparative_modeling::IgnoreSubsetConstraintSet )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_comparative_modeling_IgnoreSubsetConstraintSet )
#endif // SERIALIZATION
