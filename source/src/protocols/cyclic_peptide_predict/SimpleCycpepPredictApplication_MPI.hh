// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI.hh
/// @brief Wrapper for SimpleCycpepPredictApplication that allows the app to use hierarchical MPI/pthreads based job distribution.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)

#ifdef USEMPI

#ifndef INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_MPI_hh
#define INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_MPI_hh

#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI.fwd.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJDApplication.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace protocols {
namespace cyclic_peptide_predict {

/// @brief Wrapper for SimpleCycpepPredictApplication that allows the app to use hierarchical MPI/pthreads based job distribution.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class SimpleCycpepPredictApplication_MPI : public protocols::cyclic_peptide_predict::HierarchicalHybridJDApplication {

public:

	/// @brief Default constructor.
	SimpleCycpepPredictApplication_MPI();

	/// @brief Constructor with options
	SimpleCycpepPredictApplication_MPI(
		int const MPI_rank,
		int const MPI_n_procs,
		core::scoring::ScoreFunctionCOP sfxn_in,
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
	SimpleCycpepPredictApplication_MPI(SimpleCycpepPredictApplication_MPI const & src);

	/// @brief Destructor.
	virtual ~SimpleCycpepPredictApplication_MPI();

	/// @brief Clone function: make a copy of this object and return an owning pointer to the copy.
	HierarchicalHybridJDApplicationOP clone() const override;

protected: //Methods

	/// @brief Get the protocol-specific settings.
	/// @details The director reads these from disk and broadcasts them to all other nodes.  This function should be called from all nodes;
	/// it figures out which behaviour it should be performing.
	/// @note Pure virtual in base class; implemented here.
	void get_protocol_specific_settings() override;

	/// @brief Create an instance of the SimpleCycpepPredictApplication app class, and carry out N jobs on a single process.
	/// @details This code is called in a single thread in multi-threaded mode, and is used in the single-threaded version too.
	/// @note This implements a function that is pure virutal function in the base class.
	void derived_worker_carry_out_n_jobs(
		core::Size const njobs_from_above,
		utility::vector1 < HierarchicalHybridJD_JobResultsSummaryOP > &jobsummaries,
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

private: // Variables

	/// @brief Allowed canonical residues at each position.
	/// @details Map key 0 stores default settings applied to all positions not specified.	
	std::map< core::Size, utility::vector1 < std::string > > allowed_canonicals_;

	/// @brief Allowed noncanonical residues at each position.
	/// @details Map key 0 stores default settings applied to all positions not specified.	
	std::map< core::Size, utility::vector1 < std::string > > allowed_noncanonicals_;

	/// @brief Has an aa_composition setup file ben provided for residues in the L-alpha helix region of
	/// Ramachadran space?
	bool L_alpha_comp_file_exists_;

	/// @brief Has an aa_composition setup file ben provided for residues in the D-alpha helix region of
	/// Ramachadran space?
	bool D_alpha_comp_file_exists_;

	/// @brief Has an aa_composition setup file ben provided for residues in the L-beta strand region of
	/// Ramachadran space?
	bool L_beta_comp_file_exists_;

	/// @brief Has an aa_composition setup file ben provided for residues in the D-beta strand region of
	/// Ramachadran space?
	bool D_beta_comp_file_exists_;

	/// @brief Storage for the composition constraint setup for the L-alpha helix region of Ramachandran space.
	/// @details Cached to prevent repeated read from disk.
	std::string comp_file_contents_L_alpha_;

	/// @brief Storage for the composition constraint setup for the D-alpha helix region of Ramachandran space.
	/// @details Cached to prevent repeated read from disk.
	std::string comp_file_contents_D_alpha_;

	/// @brief Storage for the composition constraint setup for the L-beta strand region of Ramachandran space.
	/// @details Cached to prevent repeated read from disk.
	std::string comp_file_contents_L_beta_;

	/// @brief Storage for the composition constraint setup for the D-beta strand region of Ramachandran space.
	/// @details Cached to prevent repeated read from disk.
	std::string comp_file_contents_D_beta_;

	/// @brief Storage for a bin definition file.
	/// @details Cached to prevent repeated read from disk.
	std::string abba_bins_;

};

} //protocols
} //cyclic_peptide_predict

#endif //INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_MPI_hh

#endif //USEMPI
