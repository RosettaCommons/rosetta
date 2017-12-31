// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/chunk_library/MoverAndChunkLibraryJob.cc
/// @brief  The class method definitions for MoverAndChunkLibraryJob and ChunkLibraryJobResult
/// @author Andy Watkins (amw579@stanford.edu)

// Unit headers
#include <protocols/jd3/chunk_library/MoverAndChunkLibraryJob.hh>

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
namespace chunk_library {

MoverAndChunkLibraryJob::MoverAndChunkLibraryJob() {}
MoverAndChunkLibraryJob::~MoverAndChunkLibraryJob() = default;

CompletedJobOutput
MoverAndChunkLibraryJob::run()
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

	ChunkLibraryJobResultOP job_result = create_job_result();
	job_result->pose( pose_ );

	lores_pose_ = mover_->get_additional_output();
	job_result->lores_pose( lores_pose_ );
	job_result->native_pose( mover_->get_native_pose() );


	/*  if ( options_->output_filters() && ( rna_fragment_monte_carlo_ != nullptr ) ) {
	s.add_energy(  "lores_early", rna_fragment_monte_carlo_->lores_score_early() );
	if ( options_->minimize_structure() ) s.add_energy( "lores_final", rna_fragment_monte_carlo_->lores_score_final() );
	}
	*/

	// AMW: do we have to do anything for the energies here for lores? don't think so.
	EnergyJobSummaryOP summary = create_job_summary();
	summary->energy( pose_->energies().total_energy() );

	finalize_job_result( job_result );
	finalize_job_summary( summary );

	job_output.job_results.push_back( std::make_pair( summary, job_result ));

	return job_output;
}

void
MoverAndChunkLibraryJob::mover( moves::MoverOP setting ) {
	mover_ = setting;
}

void
MoverAndChunkLibraryJob::pose( core::pose::PoseOP setting ) {
	pose_ = setting;
}

moves::MoverOP MoverAndChunkLibraryJob::mover() {
	return mover_;
}

core::pose::PoseOP MoverAndChunkLibraryJob::pose() {
	return pose_;
}

core::pose::PoseCOP MoverAndChunkLibraryJob::pose() const {
	return pose_;
}

ChunkLibraryJobResultOP
MoverAndChunkLibraryJob::create_job_result() { return ChunkLibraryJobResultOP( new ChunkLibraryJobResult ); }

EnergyJobSummaryOP
MoverAndChunkLibraryJob::create_job_summary() { return EnergyJobSummaryOP( new EnergyJobSummary ); }


void MoverAndChunkLibraryJob::finalize_job_result( ChunkLibraryJobResultOP ) {//job_result ) {
	// align poses
	//rna_fragment_monte_carlo_->align_pose( pose, true /*verbose*/ );
	//
}

void MoverAndChunkLibraryJob::finalize_job_summary( EnergyJobSummaryOP ) {}


ChunkLibraryJobResult::ChunkLibraryJobResult() {}
ChunkLibraryJobResult::~ChunkLibraryJobResult() = default;

JobStatus ChunkLibraryJobResult::status() const { return jd3_job_status_success; }

core::pose::PoseOP ChunkLibraryJobResult::pose() const { return pose_; }
void ChunkLibraryJobResult::pose( core::pose::PoseOP setting ) { pose_ = setting; }

core::pose::PoseOP ChunkLibraryJobResult::lores_pose() const { return lores_pose_; }
void ChunkLibraryJobResult::lores_pose( core::pose::PoseOP setting ) { lores_pose_ = setting; }

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

} // namespace chunk_library
} // namespace jd3
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::chunk_library::ChunkLibraryJobResult::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::JobResult >( this ) );
	arc( CEREAL_NVP( pose_ ) ); // core::pose::PoseOP
	arc( CEREAL_NVP( lores_pose_ ) ); // core::pose::PoseOP
	arc( CEREAL_NVP( native_pose_ ) ); // core::pose::PoseOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::chunk_library::ChunkLibraryJobResult::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::JobResult >( this ) );
	arc( pose_ ); // core::pose::PoseOP
	arc( lores_pose_ ); // core::pose::PoseOP
	core::pose::PoseOP native_pose;
	arc( native_pose );
	native_pose_ = native_pose;
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::chunk_library::ChunkLibraryJobResult );
CEREAL_REGISTER_TYPE( protocols::jd3::chunk_library::ChunkLibraryJobResult )


/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::jd3::chunk_library::EnergyJobSummary::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::JobSummary >( this ) );
	arc( CEREAL_NVP( energy_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::jd3::chunk_library::EnergyJobSummary::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::JobSummary >( this ) );
	arc( energy_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::jd3::chunk_library::EnergyJobSummary );
CEREAL_REGISTER_TYPE( protocols::jd3::chunk_library::EnergyJobSummary )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_jd3_chunk_library_MoverAndChunkLibraryJob )
#endif // SERIALIZATION
