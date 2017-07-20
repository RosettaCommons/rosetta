// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/MoverAndFullModelJob.hh
/// @brief  The definition of class protocols::jd3::MoverAndFullModelJob, which applies a
/// Mover to a Pose in its run() method.
/// @author Andy Watkins (amw579@stanford.edu)

#ifndef INCLUDED_protocols_jd3_full_model_MoverAndFullModelJob_HH
#define INCLUDED_protocols_jd3_full_model_MoverAndFullModelJob_HH

// Unit headers
#include <protocols/jd3/full_model/MoverAndFullModelJob.fwd.hh>

// Package headers
#include <protocols/jd3/Job.hh>
#include <protocols/jd3/JobResult.hh>
#include <protocols/jd3/JobSummary.hh>

// Project headers
#include <protocols/moves/Mover.fwd.hh>
#include <core/pose/Pose.fwd.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {
namespace full_model {

class MoverAndFullModelJob : public Job
{
public:

	MoverAndFullModelJob();
	~MoverAndFullModelJob() override;

	virtual
	CompletedJobOutput run() override;

	void mover( moves::MoverOP setting );
	void pose( core::pose::PoseOP setting );

	moves::MoverOP mover();
	core::pose::PoseOP pose();
	core::pose::PoseCOP pose() const;

protected:
	/// @brief Factory method so derived classes can return their own result classes,
	/// which themselves should derive from FullModelJobResult.
	virtual FullModelJobResultOP create_job_result();

	virtual EnergyJobSummaryOP create_job_summary();

	/// @brief Method that allows derived classes to tuck data into the result object
	/// as they see fit. Noop in the base class.
	virtual void finalize_job_result( FullModelJobResultOP result );

	/// @brief Method that allows derived classes to tuck data into the job summary
	/// object as they see fit. Noop in the base class.
	virtual void finalize_job_summary( EnergyJobSummaryOP summary );

private:
	moves::MoverOP mover_;
	core::pose::PoseOP pose_;

};

class FullModelJobResult : public JobResult
{
public:
	FullModelJobResult();
	~FullModelJobResult() override;

	JobStatus status() const override;

	core::pose::PoseOP pose();
	void pose( core::pose::PoseOP setting );
private:

	core::pose::PoseOP pose_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

class EnergyJobSummary : public JobSummary
{
public:
	EnergyJobSummary();
	virtual ~EnergyJobSummary();

	core::Real energy() const;
	void energy( core::Real setting );
private:
	core::Real energy_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace full_model
} // namespace jd3
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_full_model_MoverAndFullModelJob )
#endif // SERIALIZATION


#endif //INCLUDED_protocols_jd3_full_model_MoverAndFullModelJob_HH
