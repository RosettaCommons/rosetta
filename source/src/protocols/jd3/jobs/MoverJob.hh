// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/jd3/jobs/MoverJob.hh
/// @brief A base class for any job that takes and uses a mover.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_jd3_jobs_MoverJob_HH
#define INCLUDED_protocols_jd3_jobs_MoverJob_HH

// Unit headers
#include <protocols/jd3/jobs/MoverJob.fwd.hh>
#include <protocols/jd3/job_summaries/StandardPoseJobSummary.fwd.hh>
#include <protocols/jd3/job_results/PoseJobResult.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/simple_metrics/SimpleMetric.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Protocol headers
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>

#include <protocols/jd3/CompletedJobOutput.hh>
#include <protocols/jd3/Job.hh>

//Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace jd3 {
namespace jobs {

///@brief A Generic Mover(andPose)Job.  Can be subclassed for derived jobs.
///
///@details JobSummary is a StandardPoseJobSummary that includes the energy and a SimpleMetricDataOP for any set metrics.
/// If you require anything else in your JobResult, derive from the PoseJobResult and then derive this class.
///
class MoverJob: public protocols::jd3::Job{

public:
	MoverJob();

	~MoverJob();

	///@brief Set the pose this job will run on.
	void
	pose( core::pose::PoseOP setting );


	///@brief Run the job on any set private variables (such as a PoseOP)
	/// Return a completed job output after this Job is done.
	virtual jd3::CompletedJobOutput
	run() override;


	///@brief Set the Mover this job will run.  Use a MoverContainer to use multiple movers.
	void
	set_mover( moves::MoverOP mover );

	///@brief Add a mover to a sequence mover.
	///
	///@details Overrides any single Mover set
	void
	add_mover( moves::MoverOP mover );

	///@brief Set any metrics for this particular job.
	/// The metrics will be stored in the pose so they can eventually be output.
	///
	void
	add_metric(
		core::simple_metrics::SimpleMetricCOP metric,
		std::string const & data_prefix );

	///@brief Set any metrics for this particular job.
	/// The metrics will be stored in the pose so they can eventually be output.
	///
	void
	add_metrics(
		utility::vector1< core::simple_metrics::SimpleMetricCOP > metrics,
		std::string const & data_prefix );
public:

	///@brief Get a modifieable pose
	core::pose::PoseOP
	pose();

	///@brief Geta the const pose
	core::pose::PoseCOP
	pose() const;

protected:

	///@brief PoseJobResult can be derived
	virtual job_results::PoseJobResultOP
	create_job_result( core::pose::PoseOP pose);

	///@brief StandardPoseJobSummary can be derived.
	virtual job_summaries::StandardPoseJobSummaryOP
	create_job_summary( core::pose::Pose const & pose );

	///@brief Get the mover for this job
	moves::MoverOP
	mover();

private:


	moves::MoverOP mover_ = nullptr;

	core::pose::PoseOP pose_;

	std::map< std::string, utility::vector1< core::simple_metrics::SimpleMetricCOP >> metrics_;

	///@brief Used in case the pose has does not have an energies object after running the mover.
	/// This ensures that all poses are scored, and that the JobSummary has an energy of the pose.
	core::scoring::ScoreFunctionOP scorefxn_;
};

} //protocols
} //jd3
} //jobs

#endif //protocols_jd3_jobs_MoverJob_HH
