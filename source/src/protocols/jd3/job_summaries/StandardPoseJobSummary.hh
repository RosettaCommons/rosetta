// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/job_summaries/StandardPoseJobSummary.hh
/// @brief A JobSummary for arbitrary Energies and SimpleMetrics.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_jd3_job_summaries_StandardPoseJobSummary_HH
#define INCLUDED_protocols_jd3_job_summaries_StandardPoseJobSummary_HH

// Unit headers
#include <protocols/jd3/job_summaries/StandardPoseJobSummary.fwd.hh>
#include <protocols/jd3/JobSummary.hh>
#include <core/simple_metrics/SimpleMetric.fwd.hh>
#include <core/simple_metrics/SimpleMetricData.fwd.hh>
#include <core/pose/Pose.fwd.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace jd3 {
namespace job_summaries {

///@brief A JobSummary that extracts the energy and SimpleMetricData from a pose.
///  The Job object should score the pose and run any SimpleMetrics desired.

class StandardPoseJobSummary: public protocols::jd3::JobSummary
{
public:

	StandardPoseJobSummary();

	StandardPoseJobSummary(
		core::pose::Pose const & pose
	);

	virtual ~StandardPoseJobSummary();

public:

	///////////////////////////////
	/// Data extraction
	///////////////////////////////


	///@brief
	///  1. Extracts the energy from the Energies object
	///  2. Extract a clone of the SimpleMetricData from the pose.
	virtual void
	extract_summary(
		core::pose::Pose const & pose
	);


public:
	///Data Access

	///@brief Get the energy stored here.
	/// If no energy has been set, the energy is 0.
	///
	core::Real
	energy() const;

	///@brief Get the SimpleMetric results.
	/// If no simple_metrics have been run, this will be a nullptr.
	///
	core::simple_metrics::SimpleMetricDataCOP
	metric_data() const;

public:
	///Data Overrides

	///@brief Set an arbitrary energy to use for the result.
	void
	set_energy(core::Real energy);

private:

	///@brief Extract the energy from a pose
	void
	extract_energy( core::pose::Pose const & pose );

	void
	extract_simple_metric_data( core::pose::Pose const & pose );

private:

	core::Real energy_ = 0.0;
	core::simple_metrics::SimpleMetricDataCOP metric_data_ = nullptr;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //protocols
} //jd3
} //job_summaries

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_job_summaries_StandardPoseJobSummary )
#endif // SERIALIZATION

#endif //protocols_jd3_job_summaries_StandardPoseJobSummary_HH

