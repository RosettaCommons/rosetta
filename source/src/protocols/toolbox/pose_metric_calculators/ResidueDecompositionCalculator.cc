// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/PoseMetricCalculators/ResidueDecompositionCalculator.cc
/// @brief  ResidueDecompositionCalculator class
/// @author Colin A. Smith

// Unit headers
#include <protocols/toolbox/pose_metric_calculators/ResidueDecompositionCalculator.hh>

#include <core/pose/Pose.hh>
#include <basic/MetricValue.hh>

#include <utility/string_util.hh>

#include <utility/assert.hh>

#include <utility/vector1.hh>


using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/set.hpp>
#include <cereal/types/string.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {

ResidueDecompositionCalculator::ResidueDecompositionCalculator():
	StructureDependentCalculator()
{
}

ResidueDecompositionCalculator::ResidueDecompositionCalculator(
	ResidueDecompositionCalculator const & calculator
):
	StructureDependentCalculator(),
	residue_decomposition_(calculator.residue_decomposition()),
	residue_set_numbers_(calculator.residue_set_numbers())
{
}

std::string
ResidueDecompositionCalculator::print(
	std::string const & key
) const
{
	if ( key == "num_sets" ) {
		return utility::to_string( residue_decomposition_.size() );
	} else if ( key == "residue_decomposition" ) {
		std::ostringstream sstream;
		sstream << "[";
		for ( core::Size i = 1; i <= residue_decomposition_.size(); ++i ) {
			if ( i > 1 ) sstream << ", ";
			sstream << "Set " << i << " (" << set_names_[i] << "): [";
			for ( std::set<core::Size>::const_iterator iter = residue_decomposition_[i].begin();
					iter != residue_decomposition_[i].end(); ++iter ) {
				if ( iter != residue_decomposition_[i].begin() ) sstream << ", ";
				sstream << *iter;
			}
			sstream << "]";
		}
		sstream << "]";
		return sstream.str();
	} else if ( key == "residue_set_numbers" ) {
		std::ostringstream sstream;
		sstream << "[";
		for ( core::Size i = 1; i <= residue_set_numbers_.size(); ++i ) {
			if ( i > 1 ) sstream << ", ";
			sstream << i << ": " << residue_set_numbers_[i];
		}
		sstream << "]";
		return sstream.str();
	} else if ( key == "set_names" ) {
		std::ostringstream sstream;
		sstream << "[";
		for ( core::Size i = 1; i <= set_names_.size(); ++i ) {
			if ( i > 1 ) sstream << ", ";
			sstream << set_names_[i];
		}
		sstream << "]";
		return sstream.str();
	}

	basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";
}

void
ResidueDecompositionCalculator::lookup(
	std::string const & key,
	basic::MetricValueBase * valptr
) const
{
	if ( key == "num_sets" ) {
		core::Size dummy_size(0);
		basic::check_cast( valptr, &dummy_size, "num_sets expects to return a Size" );
		(static_cast<basic::MetricValue<core::Size> *>(valptr))->set( residue_decomposition_.size() );

	} else if ( key == "residue_decomposition" ) {
		basic::check_cast( valptr, &residue_decomposition_, "decomposition expects to return a utility::vector1<std::set<core::Size> >" );
		(static_cast<basic::MetricValue<utility::vector1<std::set<core::Size> > > *>(valptr))->set( residue_decomposition_ );

	} else if ( key == "residue_set_numbers" ) {
		basic::check_cast( valptr, &residue_set_numbers_, "set_numbers expects to return a utility::vector1< core::Size >" );
		(static_cast<basic::MetricValue<utility::vector1<core::Size> > *>(valptr))->set( residue_set_numbers_ );

	} else if ( key == "set_names" ) {
		basic::check_cast( valptr, &set_names_, "set_names expects to return a utility::vector1< std::string >" );
		(static_cast<basic::MetricValue<utility::vector1<std::string> > *>(valptr))->set( set_names_ );

	} else {
		basic::Error() << "ResidueDecompositionCalculator cannot compute the requested metric " << key << std::endl;
		utility_exit();
	}
}

void
ResidueDecompositionCalculator::residue_decomposition_to_set_numbers(
	core::pose::Pose const & this_pose
)
{
	// reset residue_set_numbers_ to zeros
	residue_set_numbers_.assign(this_pose.size(), 0);

	for ( core::Size i = 1; i <= residue_decomposition_.size(); ++i ) {
		for ( std::set<core::Size>::iterator iter = residue_decomposition_[i].begin(); iter != residue_decomposition_[i].end(); ++iter ) {
			if ( *iter <= residue_set_numbers_.size() ) residue_set_numbers_[*iter] = i;
		}
	}
}

void
ResidueDecompositionCalculator::residue_set_numbers_to_decomposition()
{
	// reset residue_decomposition_ to the correct number of empty sets
	core::Size max_set_number = 0;
	for ( core::Size i = 1; i <= residue_set_numbers_.size(); ++i ) {
		core::Size const residue_set_number(residue_set_numbers_[i]);
		if ( residue_set_number > max_set_number ) max_set_number = residue_set_number;
	}
	residue_decomposition_.resize(max_set_number);
	for ( core::Size i = 1; i <= residue_decomposition_.size(); ++i ) {
		residue_decomposition_[i].clear();
	}

	for ( core::Size i = 1; i <= residue_set_numbers_.size(); ++i ) {
		core::Size const residue_set_number(residue_set_numbers_[i]);
		if ( residue_set_number > 0 ) residue_decomposition_[residue_set_number].insert(i);
	}
}


} // PoseMetricCalculators
} // toolbox
} // protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::ResidueDecompositionCalculator::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( CEREAL_NVP( residue_decomposition_ ) ); // utility::vector1<std::set<core::Size> >
	arc( CEREAL_NVP( residue_set_numbers_ ) ); // utility::vector1<core::Size>
	arc( CEREAL_NVP( set_names_ ) ); // utility::vector1<std::string>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::toolbox::pose_metric_calculators::ResidueDecompositionCalculator::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( residue_decomposition_ ); // utility::vector1<std::set<core::Size> >
	arc( residue_set_numbers_ ); // utility::vector1<core::Size>
	arc( set_names_ ); // utility::vector1<std::string>
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::toolbox::pose_metric_calculators::ResidueDecompositionCalculator );
CEREAL_REGISTER_TYPE( protocols::toolbox::pose_metric_calculators::ResidueDecompositionCalculator )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_ResidueDecompositionCalculator )
#endif // SERIALIZATION
