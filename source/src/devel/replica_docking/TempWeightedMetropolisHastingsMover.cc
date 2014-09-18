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

static thread_local basic::Tracer tr( "devel.replica_docking.TempWeightedMetropolisHastingsMover" );
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
	//	std::copy( other.weight_controllers_.begin(), other.weight_controllers_.end(), std::back_inserter( weight_controllers_ ) );
 	std::copy( other.weight_contro_1_.begin(), other.weight_contro_1_.end(), std::back_inserter( weight_contro_1_ ) );
 	std::copy( other.weight_contro_2_.begin(), other.weight_contro_2_.end(), std::back_inserter( weight_contro_2_ ) );
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
	tr.Debug << "add mover " << th_mover->get_name() << " in TempWeightedMetropolisHastings " << std::endl;
	Parent::add_mover( th_mover, weight );
	overall_weights_.push_back( weight );
	utility::vector0< utility::tag::TagCOP > const wc_tags( subtag->getTags() );
	tr.Debug << "wc_tags.size() " << wc_tags.size() << std::endl;

	bool added_1 = false;
	bool added_2 = false;
	for ( utility::vector0< utility::tag::TagCOP >::const_iterator subtag_it = wc_tags.begin(); subtag_it != wc_tags.end(); ++subtag_it ) {
		utility::tag::TagCOP const wc_tag = *subtag_it;
		TempInterpolatorBaseOP wc=NULL;
		if ( wc_tag->getName() == "Interp" && wc_tag->getOption< std::string >("key") == "weight" ) {
			core::Size dim=wc_tag->getOption< core::Size >( "dim",1 );
			tr.Debug << "dim" << dim << " nlevels_per_dim" << tempering()->nlevels_per_dim(dim) << std::endl;

			if ( tempering()->nlevels_per_dim(dim) > 1 ) {
				wc=TempInterpolatorFactory::get_instance()->new_tempInterpolator( wc_tag, tempering()->nlevels_per_dim(dim) );
				tr.Debug << "weight controller for " << dim << "th dimension is generated with " << tempering()->nlevels_per_dim(dim) << " levels" << std::endl;
			} else {
				wc=new devel::replica_docking::TempFixValue(1);
				tr.Debug << "Fix value weight controller used for " << dim << "th dimension" << std::endl;
			}

			if ( dim==1 ) {
				weight_contro_1_.push_back( wc );
				added_1 = true;
			} else {
				weight_contro_2_.push_back( wc );
				added_2 = true;
			}
			tr.Debug << "added tempInterpolator to " << dim << "th dimension weight controllers" << std::endl;
		}
	}
	tr.Debug << "added_1 " << added_1 << " added_2 " << added_2 << std::endl;
	if ( !added_1 ) {
		weight_contro_1_.push_back( new devel::replica_docking::TempFixValue( 1 ) );
		tr.Debug << "added fix value interpolator to 1st dimension " << std::endl;
	}
	if ( !added_2 ) {
		weight_contro_2_.push_back( new devel::replica_docking::TempFixValue( 1 ) );
		tr.Debug << "added fix value interpolator to 2nd dimension " << std::endl;
	}

	tr.Debug << "done adding mover " << th_mover->get_name() << std::endl;
}


ThermodynamicMoverOP
TempWeightedMetropolisHastingsMover::random_mover() const {
	tr.Trace << "random_mover" << std::endl;
	core::Size temp_level( tempering()->temperature_level() );
	GridCoord grid_coord( tempering()->level_2_grid_coord( temp_level ) );
	if ( grid_coord.size()==1 ) grid_coord.push_back( 1 ); // because of the weight_contro_2_
	bool has_changed( temp_level != last_temp_level_in_random_mover_ );
	tr.Trace << "current temp_level " << temp_level << " last temp_level " << last_temp_level_in_random_mover_ << " has_changed " << has_changed << std::endl;
	if ( has_changed ) {
		current_weighted_sampler_.clear();
		Weights::const_iterator weight_it=overall_weights_.begin();
		Interpolators::const_iterator interp_it1=weight_contro_1_.begin(); // interp_it is a pointer point to InterpolatorOP, so by *interp_it you get the real OP
		Interpolators::const_iterator interp_it2=weight_contro_2_.begin();
		runtime_assert( overall_weights_.size() == weight_contro_1_.size() );
		runtime_assert( overall_weights_.size() == weight_contro_2_.size() );

		tr.Debug << "For temp_level " << temp_level <<" the temp_weighted weights and overall weights are: " << std::endl;
		for ( ; interp_it1 != weight_contro_1_.end(); ++weight_it, ++interp_it1, ++interp_it2 ) {
			current_weighted_sampler_.add_weight( (*weight_it) * (*interp_it1)->get_value( grid_coord[1] ) * (*interp_it2)->get_value( grid_coord[2] ) );
			tr.Debug << (*interp_it1)->get_value( grid_coord[1] ) << " , " << (*interp_it2)->get_value( grid_coord[2] ) << " , "  << *weight_it << std::endl;
		}
		last_temp_level_in_random_mover_ = temp_level ;
		runtime_assert( current_weighted_sampler_.size() == overall_weights_.size() );
		runtime_assert( current_weighted_sampler_.size() > 0 );
	}
	return mover_by_index( current_weighted_sampler_.random_sample(numeric::random::rg()) ); // instead of using weighted_sampler_, we use current_weighted_sampler, which is a combination of overall_weights and tempweighted-weights.
}

} //moves
} //protocols

