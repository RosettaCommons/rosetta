// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.hh
/// @brief Application-level code for the simple_cycpep_predict app.
/// @details  This application predicts structures of simple backbone-cyclized peptides made of alpha-, beta-, or gamma-amino acids (of any chirality)
/// using generalized kinematic closure (GenKIC) for cyclization, and enforcing user-defined requiresments for numbers of mainchain hydrogen bonds.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_hh
#define INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_hh

// Unit Headers
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Project Headers
#include <protocols/cyclic_peptide/DeclareBond.fwd.hh>
#include <protocols/filters/BasicFilters.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/Ramachandran.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <stdio.h>

#define SimpleCycpepPredictApplication_PEPBOND_LENGTH 1.328685
#define SimpleCycpepPredictApplication_PEPBOND_C_ANGLE 2.02807246864
#define SimpleCycpepPredictApplication_PEPBOND_N_ANGLE 2.12406564732

namespace protocols {
namespace cyclic_peptide_predict {

/// @brief Application-level code for simple_cycpep_predict application.
/// @details Also called by the BOINC minirosetta app.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)
class SimpleCycpepPredictApplication : public utility::pointer::ReferenceCount
{
public:
	/// @brief Constructor
	///
	SimpleCycpepPredictApplication();

	/// @brief Explicit virtual destructor.
	///
	virtual ~SimpleCycpepPredictApplication();

	/// @brief Explicit copy constructor.
	///
	SimpleCycpepPredictApplication( SimpleCycpepPredictApplication const &src );

	/// @brief Register the set of options that this application uses (for the help menu).
	///
	static void register_options();

	/// @brief Initialize the application.
	/// @details Initializes using the option system.
	void initialize_from_options();

	/// @brief Actually run the application.
	/// @details The initialize_from_options() function must be called before calling this.  (Called by default constructor.)
	void run() const;

private:
	/// ------------- Methods ----------------------------

	/// @brief Actually build the geometry that we'll be working with.
	///
	void build_polymer(
		core::pose::PoseOP pose,
		utility::vector1<std::string> const &restypes
	) const;

	/// @brief Given the name of a Rama_Table_Type, set the default Rama_Table_Type.
	/// @details Error if unknown type.
	void set_default_rama_table_type( std::string const &type_name);

	/// @brief Given a string vector that we need to parse, populate the rama_table_type_by_res_ map.
	/// @details The string vector must be of the format: [integer] [rama type name] [integer] [rama type name] etc.
	/// Throws error if could not parse.
	void set_rama_table_type_by_res( utility::vector1 <std::string> const &type_name_vector);

	/// @brief Given a Rama_Table_Type name, return the Rama_Table_Type, or an informative error message on failure.
	///
	core::scoring::Rama_Table_Type get_rama_table_type_from_name( std::string const &type_name ) const;

	/// @brief Read a sequence (as a series of full names, separated by whitespace) and store
	/// it in a string vector.
	void
	read_sequence (
		std::string const &seqfile,
		utility::vector1 < std::string > &resnames
	) const;

	/// @brief Set up the DeclareBond mover used to connect the termini.
	///
	void
	set_up_termini_mover (
		protocols::cyclic_peptide::DeclareBondOP termini,
		core::pose::PoseCOP pose,
		bool const native=false
	) const;

	/// @brief Takes a vector of residue names, chooses a random number for cyclic offset, and
	/// does a cyclic permutation.
	/// @details Returns the offset and stores the new string vector in resnames_copy.
	core::Size
	do_cyclic_permutation (
		utility::vector1 <std::string> const &resnames,
		utility::vector1 <std::string> &resnames_copy
	) const;

	/// @brief Imports the native pose and sets up a terminial peptide bond.
	///
	void
	import_and_set_up_native (
		std::string const &native_file,
		core::pose::PoseOP native_pose,
		core::Size const expected_residue_count
	) const;

	/// @brief Function to add cyclic constraints to a pose.
	///
	void add_cyclic_constraints ( core::pose::PoseOP pose ) const;

	/// @brief Sets all omega values to 180, and randomizes mainchain torsions.
	/// @details For alpha-amino acids, mainchain torsions are randomized by the Ramachandran plot.
	/// For other residue types, just randomizes mainchain torsions other than peptide bonds.
	void set_mainchain_torsions (
		core::pose::PoseOP pose,
		core::Size const cyclic_offset
	) const;

	/// @brief Set up the filters for the mainchain hydrogen bonds that will
	/// be used to discard solutions with too little mainchain hydrogen bonding.
	void
	set_up_hbond_filter(
		protocols::filters::CombinedFilterOP total_hbond,
		core::Size const nres,
		core::scoring::ScoreFunctionOP sfxn,
		core::Real const &min_hbonds
	) const;

	/// @brief Use GeneralizedKIC to close the pose.
	///
	bool
	genkic_close(
		core::pose::PoseOP pose,
		core::scoring::ScoreFunctionOP sfxn_highhbond,
		protocols::filters::CombinedFilterOP total_hbond,
		core::Size const cyclic_offset
	) const;

	/// @brief Given a pose, store a list of the disulfides in the pose.
	/// @details Clears the old_disulfides list and repopulates it.
	void
	store_disulfides (
		core::pose::PoseCOP pose,
		utility::vector1 < std::pair < core::Size, core::Size > > &old_disulfides
	) const;

	/// @brief Given a pose and a list of the disulfides in the pose, break the disulfides.
	///
	void
	break_disulfides (
		core::pose::PoseOP pose,
		utility::vector1 < std::pair < core::Size, core::Size > > const &disulfides
	) const;

	/// @brief Given a pose and a list of the disulfides that should be in the pose, form the disulfides.
	///
	void
	rebuild_disulfides (
		core::pose::PoseOP pose,
		utility::vector1 < std::pair < core::Size, core::Size > > const &disulfides
	) const;

	/// @brief Given a list of old disulfide positions, generate a list of new disulfide positions based on the offset.
	/// @details Replaces the new_disulfides list.
	void
	depermute_disulfide_list(
		utility::vector1 < std::pair < core::Size, core::Size > > const &old_disulfides,
		utility::vector1 < std::pair < core::Size, core::Size > > &new_disulfides,
		core::Size const offset,
		core::Size const nres
	) const;

	/// @brief Given a pose that has undergone an N-residue cyclic permutation, restore
	/// the original pose, without the permutation.
	void
	depermute (
		core::pose::PoseOP pose,
		core::Size const offset
	) const;

	/// @brief Align pose to native_pose, and return the RMSD between the two poses.
	/// @details Assumes that the pose has already been de-permuted (i.e. the native and the pose line up).
	core::Real
	align_and_calculate_rmsd(
		core::pose::PoseOP pose,
		core::pose::PoseCOP native_pose
	) const;

	/// @brief Create a new checkpoint file.
	///
	void new_checkpoint_file() const;

	/// @brief Initialize checkpointing for this run.
	/// @details  This function does several things.  First, it checks for an existing checkpoint
	/// file.  If one exists, it checks whether the unique job name in the file matches the current
	/// job.  If it does, then this job has already been attempted, and we're somewhere in the middle
	/// of it.  The function reads the last attempt number and success count from the checkpoint
	/// file, and returns these values.  Otherwise, it creates a new checkpoint file with the current
	/// job name and returns (0,0).  If checkpointing is disabled, this function does nothing, and
	/// returns (0,0).
	/// @param[out] lastjob The index of the last job run.  Set to zero if checkpointing is disabled
	/// or if we're creating a new checkpoint file (first job run).
	/// @param[out] successes The number of successes so far.  Set to zero if checkpointing is
	/// disabled or if we're creating a new checkpoint file (first job run).
	void initialize_checkpointing( core::Size &lastjob, core::Size &successes ) const;

	/// @brief Add a checkpoint to the checkpoint file.
	/// @details  The checkpoint file must already exist.  Does nothing if checkpointing is disabled.
	/// @param[in] curjob The index of the current job just run, for writing to the checkpoint file.
	/// @param[in] successes The number of successes so far, for writing to the checkpoint file.
	void checkpoint( core::Size const curjob, core::Size const successes ) const;

	/// @brief End checkpointing and delete the checkpoint file.
	/// @details Does nothing if checkpointing is disabled.
	void end_checkpointing() const;

	/// @brief Restore the state of the random generator from a previous run.
	///
	void get_random_seed_info() const;

	/// @brief Store the state of the random generator from a previous run.
	///
	void store_random_seed_info() const;

	/// @brief Erase the stored state of the random generator from a previous run.
	///
	void erase_random_seed_info() const;

	/// @brief Given an absolute position in the sequence and the current cyclic permuation offset,
	/// return the position in the current pose corresponding to that absolute postion.
	inline core::Size current_position(
		core::Size const absolute_position,
		core::Size const permutation_offset,
		core::Size const nresidue
	) const {
		signed long int returnval( absolute_position - permutation_offset );
		if ( returnval < 1 ) returnval += nresidue;
		return static_cast<core::Size>(returnval);
	}

	/// @brief Does a position have a custom Rama table defined?
	/// @details Does not include a default custom Rama table -- only position-specific Rama tables are checked.
	inline bool custom_rama_table_defined( core::Size const absolute_position ) const {
		return ( rama_table_type_by_res_.count( absolute_position )!=0 );
	}

	/// @brief Custom Rama table for a position.
	/// @brief Returns unknown_ramatable_type if none defined for the position.
	inline core::scoring::Rama_Table_Type rama_table_type_by_res( core::Size const absolute_position ) const {
		if ( !custom_rama_table_defined( absolute_position ) ) return core::scoring::unknown_ramatable_type;
		return rama_table_type_by_res_.at(absolute_position);
	}

	/// @brief Get the default custom Rama table type.
	/// @details Defaults to unknown_ramatable_type if not set.
	inline core::scoring::Rama_Table_Type default_rama_table_type() const {
		return default_rama_table_type_;
	}

	/// @brief Are we using a rama filter?  Defaults to true.
	///
	inline bool use_rama_filter() const { return use_rama_filter_; }

private:
	/// ------------- Data -------------------------------
	/// -------- When you add new data to this class, ----
	/// -------- you must update the copy constructor ----

	/// @brief Should this application produce silent file output?
	///
	bool silent_out_;

	/// @brief The prefix for the output filename.
	/// @details Defaults to "S_".
	std::string out_filename_;

	/// @brief The output score file name.
	/// @details Defaults to "default.sc".
	std::string out_scorefilename_;

	/// @brief Filename for the text file containing the sequence of the peptide.
	/// @details Must be provided with the -cyclic_peptide:sequence_file flag.
	std::string sequence_file_;

	/// @brief The number of attempts that will be made by the generalized kinematic closure machinery.
	/// @details Defaults to 1.
	core::Size genkic_closure_attempts_;

	/// @brief The minimum number of solutions that must be found by the generalized kinematic closure machinery
	/// before a single solution is chosen.
	/// @details Defaults to 1.
	core::Size genkic_min_solution_count_;

	/// @brief Should we consider cyclic permutations?
	/// @details Defaults to true.
	bool cyclic_permutations_;

	/// @brief Should solutions be filtered by Ramachandran energy (in the case of alpha-amino acids)?
	/// @details Defaults to true.
	bool use_rama_filter_;

	/// @brief The energy cutoff for the rama filter.
	/// @details Defaults to 0.3.
	core::Real rama_cutoff_;

	/// @brief For parts of the protocol that upweight the hydrogen bond terms, this is the factor by which these terms
	/// are upweighted.
	/// @details Defaults to 10.0.
	core::Real high_hbond_weight_multiplier_;

	/// @brief The minimum number of mainchain hydrogen bonds that a generalized kinematic closure solution must have.
	/// @details Defaults to 3.
	core::Real min_genkic_hbonds_;

	/// @brief The minimum number of mainchain hydrogen bonds that a final solution must have.
	/// @details Defaults to 0 (report only).
	core::Real min_final_hbonds_;

	/// @brief The hbond energy cutoff above which a hbond is no longer counted.
	/// @details Defaults to -0.25.
	core::Real hbond_energy_cutoff_;

	/// @brief The number of FastReleax rounds to be applied whenever FastRelax is used.
	/// @details Defaults to 3.
	core::Size fast_relax_rounds_;

	/// @brief Should sidechain-mainchain hydrogen bonds be counted when counting hydrogen bonds?
	/// @details Defaults to false.
	bool count_sc_hbonds_;

	/// @brief Has the user specified a native structure?
	///
	bool native_exists_;

	/// @brief Filename for the native, if one is provided.
	///
	std::string native_filename_;

	/// @brief Max number of structures to generate.
	///
	core::Size nstruct_;

	/// @brief A unique job name for checkpointing.
	/// @details If none is provided, job is not checkpointed.
	std::string checkpoint_job_identifier_;

	/// @brief The name of the checkpoint file
	/// @details Defaults to "checkpoint.txt".  Read from options.
	std::string checkpoint_filename_;

	/// @brief A default custom Ramachandran table to be used for sampling.
	/// @details Defaults to unknown_ramatable_type (which means that no custom Rama table
	/// is applied by default).
	core::scoring::Rama_Table_Type default_rama_table_type_;

	/// @brief Custom Ramachandran tables to be used for sampling, listed by residue.
	/// @details Defaults to an empty map.  Only mapped indicies have custom rama tables
	/// applied.  Overrides default_rama_table_type_ (if used) at the residue indices in
	/// the map.
	std::map <core::Size, core::scoring::Rama_Table_Type> rama_table_type_by_res_;

	/// @brief The name of the checkpoint file for the random number generator.
	/// @details Defaults to "rng.state.gz".  Read from options.
	std::string rand_checkpoint_file_;

	/// @brief Should we try all disulfide combinations during structure prediction?
	/// @details Default false.  Read from options.
	bool try_all_disulfides_;

	/// @brief The cutoff dslf_fa13 energy, pre-relaxation, above which closure solutions are rejected.
	/// @details Read from options.
	core::Real disulf_energy_cutoff_prerelax_;

	/// @brief The cutoff dslf_fa13 energy, post-relaxation, above which closure solutions are rejected.
	/// @details Read from options.
	core::Real disulf_energy_cutoff_postrelax_;


};

} //cyclic_peptide
} //protocols

#endif //INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_hh
