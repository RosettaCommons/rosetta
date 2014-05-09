// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief TempWeightedMetropolisHastingsMover methods implemented
/// @author Zhe Zhang


// Unit Headers
#include <devel/replica_docking/TempWeightedMetropolisHastingsMover.hh>
#include <devel/replica_docking/TempWeightedMetropolisHastingsMoverCreator.hh>
#include <devel/replica_docking/TempInterpolatorFactory.hh>
#include <devel/replica_docking/TempInterpolator.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
// C++ Headers
#include <algorithm>
#include <numeric/random/random.hh>
#include <iterator>

using basic::T;
using basic::Error;
using basic::Warning;

static basic::Tracer tr( "protocols.canonical_sampling.MetropolisHastingsMover" );
static numeric::random::RandomGenerator RG(29783409);
using namespace protocols::canonical_sampling;

namespace devel {
namespace replica_docking {


std::string
TempWeightedMetropolisHastingsMoverCreator::keyname() const {
	return TempWeightedMetropolisHastingsMoverCreator::mover_name();
}

protocols::moves::MoverOP
TempWeightedMetropolisHastingsMoverCreator::create_mover() const {
	return new TempWeightedMetropolisHastingsMover;
}

std::string
TempWeightedMetropolisHastingsMoverCreator::mover_name() {
	return "TempWeightedMetropolisHastings";
}

TempWeightedMetropolisHastingsMover::TempWeightedMetropolisHastingsMover() {
	last_temp_level_in_random_mover_=99999;
}

TempWeightedMetropolisHastingsMover::TempWeightedMetropolisHastingsMover(
	TempWeightedMetropolisHastingsMover const & other
) :	Parent(other)
{
	std::copy( other.overall_weights_.begin(), other.overall_weights_.end(), std::back_inserter( overall_weights_ ) );
	std::copy( other.weight_controllers_.begin(), other.weight_controllers_.end(), std::back_inserter( weight_controllers_ ) );
}

std::string
TempWeightedMetropolisHastingsMover::get_name() const {
	return "TempWeightedMetropolisHastingsMover";
}

protocols::moves::MoverOP
TempWeightedMetropolisHastingsMover::clone() const {
	return new devel::replica_docking::TempWeightedMetropolisHastingsMover(*this);
}

protocols::moves::MoverOP
TempWeightedMetropolisHastingsMover::fresh_instance() const {
	return new TempWeightedMetropolisHastingsMover;
}

void
TempWeightedMetropolisHastingsMover::add_mover(
	ThermodynamicMoverOP th_mover,
	core::Real weight,
	utility::tag::TagCOP const& subtag
) {
	tr.Debug << "add_mover( mover, weight, tag) in TempWeightedMetropolisHastings " << std::endl;
	utility::vector0< utility::tag::TagCOP > const wc_tags( subtag->getTags() );
	TempInterpolatorBaseOP weight_controller=NULL;

	for ( utility::vector0< utility::tag::TagCOP >::const_iterator subtag_it = wc_tags.begin(); subtag_it != wc_tags.end(); ++subtag_it) {
		utility::tag::TagCOP const wc_tag = *subtag_it;
		if ( wc_tag->getName() == "Interp" && wc_tag->getOption< std::string >( "key" )=="weight" ) {
			weight_controller = TempInterpolatorFactory::get_instance()->new_tempInterpolator( wc_tag , tempering()->temperature_level() );
			tr.Debug <<"Interpolate weight_constroller" << std::endl;
			continue;
		}
	}

	if ( !weight_controller ) {
		weight_controller = new devel::replica_docking::TempFixValue(1);
		tr.Debug << "Fix Value weight_controller " << std::endl;
	}

	Parent::add_mover( th_mover, weight );
	weight_controllers_.push_back(weight_controller);
	overall_weights_.push_back( weight );
}

ThermodynamicMoverOP
TempWeightedMetropolisHastingsMover::random_mover() const {
	tr.Trace << "random_mover" << std::endl;
	core::Size temp_level( tempering()->temperature_level() );
	bool has_changed( temp_level != last_temp_level_in_random_mover_ );
	tr.Trace << "current temp_level " << temp_level << " last temp_level " << last_temp_level_in_random_mover_ << " has_changed " << has_changed << std::endl;
	if ( has_changed ) {
		current_weighted_sampler_.clear();
		Weights::const_iterator weight_it=overall_weights_.begin();
		Interpolators::const_iterator interp_it=weight_controllers_.begin(); // interp_it is a pointer point to InterpolatorOP, so by *interp_it you get the real OP
		runtime_assert( overall_weights_.size() == weight_controllers_.size() );
		tr.Debug << "For temp_level " << temp_level <<" the temp_weighted weights and overall weights are: " << std::endl;
		for ( ; interp_it != weight_controllers_.end(); ++weight_it, ++interp_it ) {
			current_weighted_sampler_.add_weight( (*interp_it)->get_value( temp_level )*(*weight_it) );
			tr.Debug << (*interp_it)->get_value( temp_level ) << " , " << *weight_it << std::endl;
		}
		last_temp_level_in_random_mover_ = temp_level ;
		runtime_assert( current_weighted_sampler_.size() == overall_weights_.size() );
		runtime_assert( current_weighted_sampler_.size() > 0 );
	}
	return mover_by_index( current_weighted_sampler_.random_sample(RG) );
}

} //moves
} //protocols

