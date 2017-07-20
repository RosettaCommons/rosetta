// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/full_model/MoverAndFullModelJob.cc
/// @brief  The class method definitions for MoverAndFullModelJob and FullModelJobResult
/// @author Andy Watkins (amw579@stanford.edu)

// Unit headers
#include <protocols/jd3/full_model/MoverAndFullModelJob.hh>

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
namespace full_model {

MoverAndFullModelJob::MoverAndFullModelJob() {}
MoverAndFullModelJob::~MoverAndFullModelJob() = default;

CompletedJobOutput
MoverAndFullModelJob::run()
{
	using namespace protocols::moves;
	mover_->apply( *pose_ );

	CompletedJobOutput job_output;

	// prepare the job status
	if ( mover_->get_last_move_status() == MS_SUCCESS ) {
		job_output.status = jd3_job_status_success;
	} else if ( mover_->get_last_move_status() == FAIL_RETRY ) {
		job_output.status = jd3_job_status_failed_retry;
	} else if ( mover_->get_last_move_status() == FAIL ) { // treat fail like fail retry?
		job_output.status = jd3_job_status_failed_retry;
	} else if ( mover_->get_last_move_status() == FAIL_DO_NOT_RETRY ) {
		job_output.status = jd3_job_status_failed_do_not_retry;
	} else if ( mover_->get_last_move_status() == FAIL_BAD_INPUT ) {
		job_output.status = jd3_job_status_inputs_were_bad;
	}

	// Retrieve all Poses created by the Mover and append them individually to
	// the CompletedJobOutput's job_results vector.
	while ( pose_ ) {
		FullModelJobResultOP job_result = create_job_result();
		job_result->pose( pose_ );

		EnergyJobSummaryOP summary = create_job_summary();
		summary->energy( pose_->energies().total_energy() );

		finalize_job_result( job_result );
		finalize_job_summary( summary );

		job_output.job_results.push_back( std::make_pair( summary, job_result ));

		// keep retrieving Poses from the Mover until the Mover stops returning any.
		// This could perhaps be expensive if the Mover wants to deliver 10K poses and
		// was not anticipating that they would all need to live in memory simultaneously!
		pose_ = mover_->get_additional_output();
	}

	return job_output;
}

void
MoverAndFullModelJob::mover( moves::MoverOP setting ) {
	mover_ = setting;
}

void
MoverAndFullModelJob::pose( core::pose::PoseOP setting ) {
	pose_ = setting;
}

moves::MoverOP MoverAndFullModelJob::mover() {
	return mover_;
}

core::pose::PoseOP MoverAndFullModelJob::pose() {
	return pose_;
}

core::pose::PoseCOP MoverAndFullModelJob::pose() const {
	return pose_;
}

FullModelJobResultOP
MoverAndFullModelJob::create_job_result() { return FullModelJobResultOP( new FullModelJobResult ); }

EnergyJobSummaryOP
MoverAndFullModelJob::create_job_summary() { return EnergyJobSummaryOP( new EnergyJobSummary ); }


void MoverAndFullModelJob::finalize_job_result( FullModelJobResultOP ) {}

void MoverAndFullModelJob::finalize_job_summary( EnergyJobSummaryOP ) {}


FullModelJobResult::FullModelJobResult() {}
FullModelJobResult::~FullModelJobResult() = default;

JobStatus FullModelJobResult::status() const { return jd3_job_status_success; }

core::pose::PoseOP FullModelJobResult::pose() { return pose_; }
void FullModelJobResult::pose( core::pose::PoseOP setting ) { pose_ = setting; }

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

} // namespace full_model
} // namespace jd3
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::full_model::FullModelJobResult::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::JobResult >( this ) );
	arc( CEREAL_NVP( pose_ ) ); // core::pose::PoseOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::full_model::FullModelJobResult::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::JobResult >( this ) );
	arc( pose_ ); // core::pose::PoseOP
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::full_model::FullModelJobResult );
CEREAL_REGISTER_TYPE( protocols::jd3::full_model::FullModelJobResult )


/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::full_model::EnergyJobSummary::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::JobSummary >( this ) );
	arc( CEREAL_NVP( energy_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::full_model::EnergyJobSummary::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::JobSummary >( this ) );
	arc( energy_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::full_model::EnergyJobSummary );
CEREAL_REGISTER_TYPE( protocols::jd3::full_model::EnergyJobSummary )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_full_model_MoverAndFullModelJob )
#endif // SERIALIZATION
