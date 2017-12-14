// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/EnergyLandscapeEvaluator.cc
/// @brief Evaluates a set of score/rms points
/// @author Tom Linsky (tlinsky@gmail.com)

#include <protocols/toolbox/EnergyLandscapeEvaluator.hh>

// Protocol headers

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility>

// C headers
#include <cmath>

static basic::Tracer TR( "protocols.toolbox.EnergyLandscapeEvaluator" );

namespace protocols {
namespace toolbox {

ScoreRmsPoint::ScoreRmsPoint( core::Real const score_val, core::Real rms_val ):
	utility::pointer::ReferenceCount(),
	score_( score_val ),
	rms_( rms_val )
{}

core::Real
ScoreRmsPoint::score() const
{
	return score_;
}

core::Real
ScoreRmsPoint::rms() const
{
	return rms_;
}

std::ostream &
operator<<( std::ostream & os, ScoreRmsPoint const & score )
{
	os << score.score_ << " " << score.rms_;
	return os;
}

ScoreRmsPoints::ScoreRmsPoints( ScoreRmsPoint const & bg_val ):
	utility::vector1< ScoreRmsPoint >(),
	bg_( bg_val )
{}

ScoreRmsPoint const &
ScoreRmsPoints::bg() const
{
	return bg_;
}

EnergyLandscapeEvaluator::EnergyLandscapeEvaluator():
	utility::pointer::ReferenceCount()
{

}

EnergyLandscapeEvaluator::~EnergyLandscapeEvaluator()= default;

RotamerBoltzmannWeightEvaluator::RotamerBoltzmannWeightEvaluator( core::Real const temperature, bool const include_bg_point_in_sum ):
	EnergyLandscapeEvaluator(),
	temperature_( temperature ),
	include_bg_( include_bg_point_in_sum )
{}

RotamerBoltzmannWeightEvaluator::~RotamerBoltzmannWeightEvaluator() = default;

EnergyLandscapeEvaluatorOP
RotamerBoltzmannWeightEvaluator::clone() const
{
	return EnergyLandscapeEvaluatorOP( new RotamerBoltzmannWeightEvaluator( *this ) );
}

core::Real
RotamerBoltzmannWeightEvaluator::compute( ScoreRmsPoints const & points ) const
{
	core::Real boltz_sum = 0.0;
	if ( include_bg_ ) {
		boltz_sum += 1.0; // exp(0) == 1.0
	}
	for ( auto p=points.begin(); p!=points.end(); ++p ) {
		boltz_sum += exp( ( points.bg().score() - p->score() ) / temperature_ );
	}

	if ( boltz_sum < 1e-6 ) return 0.0;

	return ( 1 / boltz_sum );
}


MulliganPNearEvaluator::MulliganPNearEvaluator( core::Real const temperature, core::Real const lambda ):
	EnergyLandscapeEvaluator(),
	temperature_( temperature ),
	lambda_sq_( lambda*lambda )
{}

MulliganPNearEvaluator::~MulliganPNearEvaluator() = default;

EnergyLandscapeEvaluatorOP
MulliganPNearEvaluator::clone() const
{
	return EnergyLandscapeEvaluatorOP( new MulliganPNearEvaluator( *this ) );
}

core::Real
MulliganPNearEvaluator::compute( ScoreRmsPoints const & points ) const
{
	core::Real numerator = 1.0; // exp( -(0*0)/lambda_sq )*exp( 0/temp ) == 1.0
	core::Real denominator = 1.0; // exp( 0/temp ) == 1.0
	for ( auto p=points.begin(); p!=points.end(); ++p ) {
		core::Real const energy = exp( ( points.bg().score() - p->score() ) / temperature_ );
		core::Real const nearness = exp( -(p->rms()*p->rms())/(lambda_sq_) );
		numerator += energy * nearness;
		denominator += energy;
	}
	TR.Debug << "PNearEvaluator: numerator: " << numerator << " denominator: " << denominator << std::endl;
	if ( denominator < 1e-6 ) return 0.0;
	return numerator / denominator;
}

/*
////////////////  Modified DDG Evaluator

ModifiedDdgEvaluator::ModifiedDdgEvaluator( ScoreRmsPoints const & ddg_values ):
EnergyLandscapeEvaluator(),
ddg_values_( ddg_values )
{}

ModifiedDdgEvaluator::~ModifiedDdgEvaluator()
{}

EnergyLandscapeEvaluatorOP
ModifiedDdgEvaluator::clone() const
{
return EnergyLandscapeEvaluatorOP( new ModifiedDdgEvaluator( *this ) );
}

core::Real
ModifiedDdgEvaluator::compute( ScoreRmsPoints const & points ) const
{
return 0.0;
}
*/

} //protocols
} //toolbox

