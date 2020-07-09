// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HelicalBundlePredictApplication_MPI.hh
/// @brief MPI implementation of HelicalBundlePredictApplication.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifdef USEMPI

#ifndef INCLUDED_protocols_helical_bundle_predict_HelicalBundlePredictApplication_MPI_hh
#define INCLUDED_protocols_helical_bundle_predict_HelicalBundlePredictApplication_MPI_hh

// Project headers
#include <protocols/helical_bundle_predict/HelicalBundlePredictApplication_MPI.fwd.hh>

// Protocols headers
#include <protocols/helical_bundle_predict/HBP_MoveGenerator.fwd.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJDApplication.hh>
#include <protocols/helical_bundle_predict/HelicalBundlePredictApplication.hh>

// Core headers
#include <core/scoring/ScoreFunction.fwd.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

// Multithreading headers
#ifdef MULTI_THREADED
#include <mutex>
#endif

namespace protocols {
namespace helical_bundle_predict {

/// @brief MPI implementation of HelicalBundlePredictApplication.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class HelicalBundlePredictApplication_MPI : public protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication {

public:

	/// @brief Default constructor.
	HelicalBundlePredictApplication_MPI();

	/// @brief Constructor with options
	///
	HelicalBundlePredictApplication_MPI(
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
	);

	/// @brief Copy constructor.
	HelicalBundlePredictApplication_MPI(HelicalBundlePredictApplication_MPI const & src);

	/// @brief Destructor.
	virtual ~HelicalBundlePredictApplication_MPI();

	/// @brief Clone function: create a copy of this object and return an owning pointer to the copy.
	protocols::cyclic_peptide_predict::HierarchicalHybridJDApplicationOP clone() const override;

	/// @brief Set the options for this application.
	void set_options( protocols::helical_bundle_predict::HelicalBundlePredictApplicationOptionsCOP options );

protected: //Methods

	/// @brief Get the protocol-specific settings.
	/// @details The director reads these from disk and broadcasts them to all other nodes.  This function should be called from all nodes;
	/// it figures out which behaviour it should be performing.
	/// @note Pure virtual.  Must be implemented by derived classes.
	void get_protocol_specific_settings() override;

	/// @brief Create an instance of the appropriate app class, and carry out N jobs on a single process.
	/// @details This code is called in a single thread in multi-threaded mode, and is used in the single-threaded version too.
	/// @note Pure virutal function must be implemented by derived classes.
	void derived_worker_carry_out_n_jobs(
		core::Size const njobs_from_above,
		utility::vector1 < protocols::cyclic_peptide_predict::HierarchicalHybridJD_JobResultsSummaryOP > &jobsummaries,
		utility::vector1 < core::io::silent::SilentStructOP > &all_output,
		core::scoring::ScoreFunctionOP sfxn,
		core::pose::PoseCOP native,
		std::string const &sequence		
	) const override;

	/// @brief Compute the RMSD between a pose and a reference pose.
	/// @details Must be implemented by derived classes, since this might be done differently for
	/// different classes of molecule.
	core::Real
	derived_worker_compute_rmsd(
		core::pose::Pose const & pose,
		core::pose::Pose const & reference_pose,
		std::string const &sequence
	) const override;

private:

	/// @brief The move generator used in centroid move.
	/// @details Initialized once at the start of the protocol (which requires a read from disk).  Cloned to avoid
	/// repeated reads from disk.
	HBP_MoveGeneratorOP centroid_move_generator_;

	/// @brief A scorefunction to use for fullatom refinement
	core::scoring::ScoreFunctionOP sfxn_fullatom_;

#ifdef MULTI_THREADED
	/// @brief A mutex for the centroid move generator, to facilitate threadsafe creation and cloning.
	mutable std::mutex centroid_move_generator_mutex_;
	/// @brief A mutex for the fullatom refinement scorefunction, to facilitate threadsafe creation and cloning.
	mutable std::mutex sfxn_fullatom_mutex_;
#endif //MULTI_THREADED

	/// @brief Options for helical bundle structure prediction.
	HelicalBundlePredictApplicationOptionsCOP options_;

};

} //protocols
} //helical_bundle

#endif //INCLUDED_protocols_helical_bundle_predict_HelicalBundlePredictApplication_MPI_hh

#endif //USEMPI
