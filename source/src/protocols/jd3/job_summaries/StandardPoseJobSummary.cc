// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/job_summaries/StandardPoseJobSummary.cc
/// @brief A JobSummary for arbitrary SimpleMetrics.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

// Unit headers
#include <protocols/jd3/job_summaries/StandardPoseJobSummary.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/SimpleMetricData.hh>
#include <core/simple_metrics/util.hh>

// Basic/Utility headers

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace jd3 {
namespace job_summaries {

using namespace protocols::jd3;
using namespace core::simple_metrics;

//Constructor
StandardPoseJobSummary::StandardPoseJobSummary()
{}

StandardPoseJobSummary::StandardPoseJobSummary( core::pose::Pose const & pose ){
	extract_summary(pose);
}

//Destructor
StandardPoseJobSummary::~StandardPoseJobSummary() = default;

core::Real
StandardPoseJobSummary::energy() const {
	return energy_;
}

SimpleMetricDataCOP
StandardPoseJobSummary::metric_data() const {
	return metric_data_;
}

void
StandardPoseJobSummary::extract_summary(core::pose::Pose const & pose){
	extract_energy(pose);
	extract_simple_metric_data( pose );
}

void
StandardPoseJobSummary::extract_energy(core::pose::Pose const & pose){
	if ( pose.energies().energies_updated() ) {
		energy_ = pose.energies().total_energy();
	}
}

void
StandardPoseJobSummary::extract_simple_metric_data(const core::pose::Pose &pose){
	if ( has_sm_data(pose) ) {
		metric_data_ = SimpleMetricDataOP( new SimpleMetricData( * get_sm_data(pose) ) );
	}
}

void
StandardPoseJobSummary::set_energy(core::Real energy){
	energy_ = energy;
}

} //protocols
} //jd3
} //job_summaries

#ifdef    SERIALIZATION
template< class Archive >
void
protocols::jd3::job_summaries::StandardPoseJobSummary::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::JobSummary >( this ) );
	arc( energy_ ); // core::Real
	core::simple_metrics::SimpleMetricDataOP metric_data;
	arc( metric_data);
	metric_data_ = metric_data;
}

template< class Archive >
void
protocols::jd3::job_summaries::StandardPoseJobSummary::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::JobSummary >( this ) );
	arc( CEREAL_NVP( energy_ ) ); // core::Real
	arc( CEREAL_NVP( metric_data_ ) );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::job_summaries::StandardPoseJobSummary );
CEREAL_REGISTER_TYPE( protocols::jd3::job_summaries::StandardPoseJobSummary )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_job_summaries_StandardPoseJobSummary )

#endif // SERIALIZATION


