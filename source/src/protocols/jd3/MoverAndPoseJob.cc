// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/MoverAndPoseJob.cc
/// @brief  The class method definitions for MoverAndPoseJob and PoseJobResult
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit headers
#include <protocols/jd3/MoverAndPoseJob.hh>

// Project headers
#include <protocols/moves/Mover.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>

//utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace jd3 {

MoverAndPoseJob::MoverAndPoseJob() {}
MoverAndPoseJob::~MoverAndPoseJob() = default;

CompletedJobOutput
MoverAndPoseJob::run()
{
	using namespace protocols::moves;

	mover_->apply( *pose_ );
	PoseJobResultOP job_result = create_job_result();
	job_result->pose( pose_ );

	EnergyJobSummaryOP summary = create_job_summary();
	summary->energy( pose_->energies().total_energy() );
	if ( mover_->get_last_move_status() == MS_SUCCESS ) {
		summary->status( jd3_job_status_success );
	} else if ( mover_->get_last_move_status() == FAIL_RETRY ) {
		summary->status( jd3_job_status_failed_retry );
	} else if ( mover_->get_last_move_status() == FAIL ) { // treat fail like fail retry?
		summary->status( jd3_job_status_failed_retry );
	} else if ( mover_->get_last_move_status() == FAIL_DO_NOT_RETRY ) {
		summary->status( jd3_job_status_failed_do_not_retry );
	} else if ( mover_->get_last_move_status() == FAIL_BAD_INPUT ) {
		summary->status( jd3_job_status_inputs_were_bad );
	}

	finalize_job_result( job_result );
	finalize_job_summary( summary );

	CompletedJobOutput job_output;
	job_output.first  = summary;
	job_output.second = job_result;
	return job_output;
}

void
MoverAndPoseJob::mover( moves::MoverOP setting ) {
	mover_ = setting;
}

void
MoverAndPoseJob::pose( core::pose::PoseOP setting ) {
	pose_ = setting;
}

moves::MoverOP MoverAndPoseJob::mover() {
	return mover_;
}

core::pose::PoseOP MoverAndPoseJob::pose() {
	return pose_;
}

core::pose::PoseCOP MoverAndPoseJob::pose() const {
	return pose_;
}

PoseJobResultOP
MoverAndPoseJob::create_job_result() { return PoseJobResultOP( new PoseJobResult ); }

EnergyJobSummaryOP
MoverAndPoseJob::create_job_summary() { return EnergyJobSummaryOP( new EnergyJobSummary ); }


void MoverAndPoseJob::finalize_job_result( PoseJobResultOP ) {}

void MoverAndPoseJob::finalize_job_summary( EnergyJobSummaryOP ) {}


PoseJobResult::PoseJobResult() {}
PoseJobResult::~PoseJobResult() = default;

JobStatus PoseJobResult::status() const { return jd3_job_status_success; }

core::pose::PoseOP PoseJobResult::pose() { return pose_; }
void PoseJobResult::pose( core::pose::PoseOP setting ) { pose_ = setting; }

EnergyJobSummary::EnergyJobSummary() : energy_( 0.0 ) {}
EnergyJobSummary::~EnergyJobSummary() {}

core::Real
EnergyJobSummary::energy() const
{
	return energy_;
}

void
EnergyJobSummary::energy( core::Real setting )
{
	energy_ = setting;
}

} // namespace jd3
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::PoseJobResult::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::JobResult >( this ) );
	arc( CEREAL_NVP( pose_ ) ); // core::pose::PoseOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::PoseJobResult::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::JobResult >( this ) );
	arc( pose_ ); // core::pose::PoseOP
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::PoseJobResult );
CEREAL_REGISTER_TYPE( protocols::jd3::PoseJobResult )


/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::EnergyJobSummary::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::JobSummary >( this ) );
	arc( CEREAL_NVP( energy_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::EnergyJobSummary::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::JobSummary >( this ) );
	arc( energy_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::EnergyJobSummary );
CEREAL_REGISTER_TYPE( protocols::jd3::EnergyJobSummary )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_MoverAndPoseJob )
#endif // SERIALIZATION
