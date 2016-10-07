// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ScoreFunction.cc
/// @brief  ScoreFunction class definition. A easy way to deal with the low-res score function problem in replica docking
/// @author Zhe Zhang


// Unit headers
#include <core/scoring/MinScoreScoreFunction.hh>

// Package headers

// // Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/Jump.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/Energies.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>

// Utility headers
#include <utility/options/OptionCollection.hh>
#include <utility/options/keys/OptionKeyList.hh>

// C++ headers
#include <limits>


static THREAD_LOCAL basic::Tracer tr( "core.scoring.MinScoreScoreFunction" );

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace scoring {

///////////////////////////////////////////////////////////////////////////////
MinScoreScoreFunction::MinScoreScoreFunction():
	MinScoreScoreFunction( basic::options::option )
{}

MinScoreScoreFunction::MinScoreScoreFunction( utility::options::OptionCollection const & options ):
	ScoreFunction( options ),
	// RM 2013-11-18: Make it work like regular scorefunction if min not specified.
	min_score_( -1 * std::numeric_limits< core::Real >::max() )
{}

///////////////////////////////////////////////////////////////////////////////
ScoreFunctionOP
MinScoreScoreFunction::clone() const
{
	MinScoreScoreFunctionOP newscorefxn( new MinScoreScoreFunction );
	newscorefxn->assign( *this );
	return newscorefxn;
}

///////////////////////////////////////////////////////////////////////////////
void
MinScoreScoreFunction::assign( ScoreFunction const & src )
{
	ScoreFunction::assign( src );
	min_score_ = -1 * std::numeric_limits< core::Real >::max();
}

void
MinScoreScoreFunction::assign( MinScoreScoreFunction const & src )
{
	ScoreFunction::assign( src );
	// MinScoreScoreFunction specific values
	min_score_ = src.min_score_;
}

///////////////////////////////////////////////////////////////////////////////

/// @details delegate to the option-collection accepting version of this constructor
MinScoreScoreFunction::MinScoreScoreFunction( ScoreFunction const & src, core::Real min_score ) :
	MinScoreScoreFunction( src, min_score, basic::options::option )
{}

MinScoreScoreFunction::MinScoreScoreFunction( ScoreFunction const & src, core::Real min_score, utility::options::OptionCollection const & options  ) :
	ScoreFunction( options )
{
	ScoreFunction::assign( src );
	min_score_ = min_score;
}

/// @details delegate to the option-collection accepting version of this constructor
MinScoreScoreFunction::MinScoreScoreFunction( core::Real min_score ) :
	MinScoreScoreFunction( min_score, basic::options::option )
{
}

MinScoreScoreFunction::MinScoreScoreFunction( core::Real min_score, utility::options::OptionCollection const & options ) :
	ScoreFunction( options )
{
	min_score_ = min_score;
}

///////////////////////////////////////////////////////////////////////////////

// to start out, just thinking fullatom energies
//
// NOTE: no freakin rotamer trials inside scoring!
Real
MinScoreScoreFunction::operator()( pose::Pose & pose ) const
{
	ScoreFunction::operator()( pose ); //score -- but without atom_pair_constraints..
	// is probably cheaper to not apply a completely new scorefunction...
	EnergyMap cst_free_weights( weights() );
	Real cst_weight = cst_free_weights.get( atom_pair_constraint );
	cst_free_weights[ atom_pair_constraint ] = 0;
	Real uncst_energy = pose.energies().total_energies().dot( cst_free_weights );
	Real min_energy = uncst_energy < min_score_ ? min_score_ : uncst_energy;
	tr.Debug << "uncst_energy: " << uncst_energy << " min_energy: " << min_energy << std::endl;
	pose.energies().total_energies()[ total_score ] = min_energy + pose.energies().total_energies()[ atom_pair_constraint ]*cst_weight;
	pose::setPoseExtraScore( pose, "min_score", uncst_energy );
	return pose.energies().total_energies()[ total_score ];
}

void
MinScoreScoreFunction::list_options_read( utility::options::OptionKeyList & options_read )
{
	ScoreFunction::list_options_read( options_read );
}

///////////////////////////////////////////////////////////////////////////////
} // namespace scoring
} // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::MinScoreScoreFunction::save( Archive & arc ) const {
	arc( cereal::base_class< core::scoring::ScoreFunction >( this ) );
	arc( CEREAL_NVP( min_score_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::MinScoreScoreFunction::load( Archive & arc ) {
	arc( cereal::base_class< core::scoring::ScoreFunction >( this ) );
	arc( min_score_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::MinScoreScoreFunction );
CEREAL_REGISTER_TYPE( core::scoring::MinScoreScoreFunction )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_MinScoreScoreFunction )
#endif // SERIALIZATION
