// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HelicalBundlePredictApplication_MPI.cc
/// @brief MPI implementation of HelicalBundlePredictApplication.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef USEMPI

// Project headers:
#include <protocols/helical_bundle_predict/HelicalBundlePredictApplication_MPI.hh>

// Core headers:
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>

// Protocols headers
#include <protocols/helical_bundle_predict/HBP_MoveGenerator.hh>
#include <protocols/helical_bundle_predict/HBP_HelixCoilMoveGenerator.hh>
#include <protocols/helical_bundle_predict/HelicalBundlePredictApplication.hh>

// Basic headers:
#include <basic/Tracer.hh>

// Utility headers:
#include <utility/pointer/memory.hh>

static basic::Tracer TR( "protocols.helical_bundle_predict.HelicalBundlePredictApplication_MPI" );
static basic::Tracer TR_summary( "protocols.helical_bundle_predict.HelicalBundlePredictApplication_MPI_summary" );


namespace protocols {
namespace helical_bundle_predict {

/// @brief Default constructor.
HelicalBundlePredictApplication_MPI::HelicalBundlePredictApplication_MPI():
	protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication(TR, TR_summary),
	centroid_move_generator_(nullptr),
#ifdef MULTI_THREADED
	centroid_move_generator_mutex_(),
	sfxn_fullatom_mutex_(),
#endif
	options_(nullptr)
{}

/// @brief Constructor with options
///
HelicalBundlePredictApplication_MPI::HelicalBundlePredictApplication_MPI(
	int const MPI_rank,
	int const MPI_n_procs,
	core::scoring::ScoreFunctionCOP centroid_sfxn_in,
	core::scoring::ScoreFunctionCOP fullatom_sfxn_in,
	core::Size const total_hierarchy_levels,
	utility::vector1 < core::Size > const & procs_per_hierarchy_level,
	utility::vector1 < core::Size > const &batchsize_per_level,
	std::string const &sort_type,
	bool const select_highest,
	core::Real const &output_fraction,
	std::string const &output_filename,
	core::Real const &lambda,
	core::Real const &kbt,
	bool const compute_rmsd_to_lowest,
	core::Real const & compute_pnear_to_lowest_fract,
	bool const compute_sasa_metrics,
	core::Size const threads_per_worker_proc //Only used in multi-threaded build.
) :
	protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication(
		TR, TR_summary, MPI_rank, MPI_n_procs, centroid_sfxn_in, total_hierarchy_levels,
		procs_per_hierarchy_level, batchsize_per_level, sort_type, select_highest,
		output_fraction, output_filename, lambda, kbt, compute_rmsd_to_lowest,
		compute_pnear_to_lowest_fract, compute_sasa_metrics, threads_per_worker_proc
	),
	centroid_move_generator_(nullptr),
	sfxn_fullatom_(fullatom_sfxn_in->clone()),
#ifdef MULTI_THREADED
	centroid_move_generator_mutex_(),
	sfxn_fullatom_mutex_(),
#endif
	options_(nullptr)
{}

/// @brief Copy constructor.
HelicalBundlePredictApplication_MPI::HelicalBundlePredictApplication_MPI( HelicalBundlePredictApplication_MPI const &src ):
	protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication( src ),
	centroid_move_generator_( src.centroid_move_generator_ == nullptr ? nullptr : src.centroid_move_generator_->clone() ),
#ifdef MULTI_THREADED
	centroid_move_generator_mutex_(),
	sfxn_fullatom_mutex_(),
#endif
	options_( src.options_ == nullptr ? nullptr : src.options_->clone())
{}

/// @brief Destructor.
HelicalBundlePredictApplication_MPI::~HelicalBundlePredictApplication_MPI(){}

/// @brief Clone function: create a copy of this object and return an owning pointer to the copy.
protocols::cyclic_peptide_predict::HierarchicalHybridJDApplicationOP
HelicalBundlePredictApplication_MPI::clone() const {
	return utility::pointer::make_shared< HelicalBundlePredictApplication_MPI >(*this);
}

/// @brief Set the options for this application.
void
HelicalBundlePredictApplication_MPI::set_options(
	protocols::helical_bundle_predict::HelicalBundlePredictApplicationOptionsCOP options
) {
	options_ = options;
}

//////////////////////////////////////////////////////////////// PROTECTED FUNCTIONS ////////////////////////////////////////////////////////////////

/// @brief Get the protocol-specific settings.
/// @details The director reads these from disk and broadcasts them to all other nodes.  This function should be called from all nodes;
/// it figures out which behaviour it should be performing.
/// @note Implements a pure virtual function in the base class.
void
HelicalBundlePredictApplication_MPI::get_protocol_specific_settings() {
	if( options_ == nullptr ) options_ = utility::pointer::make_shared< HelicalBundlePredictApplicationOptions >(); //Initializes from options system, once.
	centroid_move_generator_ = HelicalBundlePredictApplication::create_centroid_move_generator(); //Currently, makes a HBP_HelixCoilMoveGenerator.  Creator function is invoked since I might have different types in the future.
	HBP_HelixCoilMoveGeneratorOP helix_coil_generator( utility::pointer::dynamic_pointer_cast< HBP_HelixCoilMoveGenerator >( centroid_move_generator_ ) );
	runtime_assert( helix_coil_generator != nullptr );
	helix_coil_generator->set_up_user_helix_assignments( options_->helix_assignment_file_contents() );
}

/// @brief Create an instance of the appropriate app class, and carry out N jobs on a single process.
/// @details This code is called in a single thread in multi-threaded mode, and is used in the single-threaded version too.
/// @note Implements a pure virtual function in the base class.
void
HelicalBundlePredictApplication_MPI::derived_worker_carry_out_n_jobs(
	core::Size const njobs_from_above,
	utility::vector1 < protocols::cyclic_peptide_predict::HierarchicalHybridJD_JobResultsSummaryOP > &jobsummaries,
	utility::vector1 < core::io::silent::SilentStructOP > &all_output,
	core::scoring::ScoreFunctionOP sfxn,
	core::pose::PoseCOP native,
	std::string const & //sequence	
) const {

	// First, copy the centroid move generator object.  We lock a mutex to do this in multithreaded mode:
#ifdef MULTI_THREADED
	centroid_move_generator_mutex_.lock();
#endif //MULTI_THREADED
	HBP_MoveGeneratorOP centroid_move_generator_copy( centroid_move_generator_->clone() );
#ifdef MULTI_THREADED
	centroid_move_generator_mutex_.unlock();
#endif //MULTI_THREADED

#ifdef MULTI_THREADED
	sfxn_fullatom_mutex_.lock();
#endif //MULTI_THREADED
	runtime_assert_string_msg( sfxn_fullatom_  != nullptr, "Error in HelicalBundlePredictApplication_MPI::derived_worker_carry_out_n_jobs(): No fullatom scoring function was provided." );
	core::scoring::ScoreFunctionOP sfxn_fullatom_local( sfxn_fullatom_->clone() );
#ifdef MULTI_THREADED
	sfxn_fullatom_mutex_.unlock();
#endif //MULTI_THREADED

	HelicalBundlePredictApplicationOP predict_app( utility::pointer::make_shared< HelicalBundlePredictApplication >( options_, centroid_move_generator_copy, sfxn, sfxn_fullatom_local ) );
	predict_app->set_output( &jobsummaries, &all_output );
	predict_app->set_my_rank( MPI_rank() );
	predict_app->set_nstruct( njobs_from_above );
	predict_app->set_already_completed_job_count( worker_job_count() );
	predict_app->set_native( native );
	predict_app->run();
}

/// @brief Compute the RMSD between a pose and a reference pose.
/// @details Must be implemented by derived classes, since this might be done differently for
/// different classes of molecule.
core::Real
HelicalBundlePredictApplication_MPI::derived_worker_compute_rmsd(
	core::pose::Pose const & pose,
	core::pose::Pose const & reference_pose,
	std::string const &//sequence
) const {

#ifdef MULTI_THREADED
	centroid_move_generator_mutex_.lock();
#endif //MULTI_THREADED
	HBP_MoveGeneratorOP centroid_move_generator_copy( centroid_move_generator_->clone() );
#ifdef MULTI_THREADED
	centroid_move_generator_mutex_.unlock();
#endif //MULTI_THREADED

#ifdef MULTI_THREADED
	sfxn_fullatom_mutex_.lock();
#endif //MULTI_THREADED
	runtime_assert_string_msg( sfxn_fullatom_  != nullptr, "Error in HelicalBundlePredictApplication_MPI::derived_worker_carry_out_n_jobs(): No fullatom scoring function was provided." );
	core::scoring::ScoreFunctionOP sfxn_fullatom_local( sfxn_fullatom_->clone() );
#ifdef MULTI_THREADED
	sfxn_fullatom_mutex_.unlock();
#endif //MULTI_THREADED

	HelicalBundlePredictApplicationOP predict_app(
		utility::pointer::make_shared< HelicalBundlePredictApplication >( options_, centroid_move_generator_copy, sfxn_fullatom_local, sfxn_fullatom_local )
	);
	predict_app->set_my_rank( MPI_rank() );
	predict_app->set_native( reference_pose.clone() );
	return predict_app->align_to_native_pose( *( pose.clone() ) );
}

} //protocols
} //helical_bundle

#endif //USEMPI
