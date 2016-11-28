// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.hh
/// @brief Application-level code for the simple_cycpep_predict app.
/// @details  This application predicts structures of simple backbone-cyclized peptides made of alpha-, beta-, or gamma-amino acids (of any chirality)
/// using generalized kinematic closure (GenKIC) for cyclization, and enforcing user-defined requiresments for numbers of mainchain hydrogen bonds.
/// @author Vikram K. Mulligan, Baker laboratory (vmullig@uw.edu)

#ifndef INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_hh
#define INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_hh

// Unit Headers
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.fwd.hh>
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication_MPI_JobResultsSummary.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Project Headers
#include <protocols/cyclic_peptide/DeclareBond.fwd.hh>
#include <protocols/denovo_design/movers/FastDesign.fwd.hh>
#include <protocols/filters/BasicFilters.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/io/silent/SilentStruct.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <cstdio>
#include <map>

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
	/// @details If allow_file_read is true, initialization triggers reads
	/// from the filesystem.
	SimpleCycpepPredictApplication( bool const allow_file_read = true );

	/// @brief Explicit virtual destructor.
	///
	~SimpleCycpepPredictApplication() override;

	/// @brief Explicit copy constructor.
	///
	SimpleCycpepPredictApplication( SimpleCycpepPredictApplication const &src );

	/// @brief Register the set of options that this application uses (for the help menu).
	///
	static void register_options();

	/// @brief Initialize the application.
	/// @details Initializes using the option system.
	void initialize_from_options();

	/// @brief Sets the default scorefunction to use.
	/// @details The scorefunction is cloned.  The high-hbond version is constructed
	/// from this one.  If necessary, the aa_composition score term will be turned on
	/// in that one; it needn't be turned on in this one.
	void set_scorefxn( core::scoring::ScoreFunctionCOP sfxn_in );

	/// @brief Allows external code to provide a native, so that the SimpleCycpepPredictApplication doesn't have to read
	/// directly from disk.
	void set_native( core::pose::PoseCOP native );

	/// @brief Allows external code to provide a sequence, so that the SimpleCycpepPredictApplication doesn't have to read
	/// directly from disk.
	void set_sequence( std::string const &seq );

	/// @brief Allows external code to set the allowed residues by position, so that this needn't be read directly
	/// from disk.
	void set_allowed_residues_by_position( std::map< core::Size, utility::vector1< std::string > > const &allowed_canonicals, std::map< core::Size, utility::vector1< std::string > > const &allowed_noncanonicals );

	/// @brief Allows external code to specify that output should be appended to a list of SilentStructureOPs, so that the
	/// SimpleCycpepPredictApplication doesn't have to write directly to disk.
	void set_silentstructure_outputlist( utility::vector1 < core::io::silent::SilentStructOP > * silentlist, utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > * summarylist );

	/// @brief Allows external code to suppress checkpointing, to prevent direct file I/O from disk.
	/// @details Useful on Blue Gene.
	void set_suppress_checkpoints( bool const suppress_checkpoints);

	/// @brief If called by MPI code, the rank of the current process can be stored here.
	/// @details Used for output of job summaries.
	void set_my_rank( int const rank_in );

	/// @brief Allows external code to override the number of structures that this should generate (otherwise
	/// set by options system.
	void set_nstruct( core::Size const nstruct_in );

	/// @brief Allows external code to set the file contents for the L-alpha aa_composition file.
	///
	void set_L_alpha_compfile_contents( std::string const &contents_in );

	/// @brief Allows external code to set the file contents for the D-alpha aa_composition file.
	///
	void set_D_alpha_compfile_contents( std::string const &contents_in );

	/// @brief Allows external code to set the file contents for the L-beta aa_composition file.
	///
	void set_L_beta_compfile_contents( std::string const &contents_in );

	/// @brief Allows external code to set the file contents for the D-beta aa_composition file.
	///
	void set_D_beta_compfile_contents( std::string const &contents_in );

	/// @brief Allows external code to set the ABBA bin parameters without having to read the
	/// bin params file directly from disk.
	void set_abba_bins_binfile_contents( std::string const &contents_in );

	/// @brief Set the frequency with which we sample cis proline.
	/// @details Implicitly sets sample_cis_pro_ to "true" if freq_in is not 0.0, "false" if it is.
	void set_sample_cis_pro_frequency( core::Real const &freq_in );

	/// @brief Set cis proline sampling OFF.
	///
	void disable_cis_pro_sampling();

	/// @brief Set the total energy cutoff.
	/// @details Also sets use_total_energy_cutoff_ to 'true'.
	void set_total_energy_cutoff( core::Real const &value_in );

	/// @brief Sets use_total_energy_cutoff_ to 'false'.
	///
	void disable_total_energy_cutoff();

	/// @brief Set the number of rounds of relaxation with flexible
	/// bond angles.
	void set_angle_relax_rounds( core::Size const rounds_in );

	/// @brief Set the number of rounds of relaxation with flexible
	/// bond angles and bond lengths.
	void set_angle_length_relax_rounds( core::Size const rounds_in );

	/// @brief Set the number of rounds of Cartesian relaxation.
	///
	void set_cartesian_relax_rounds( core::Size const rounds_in );

	/// @brief Set whether we're using RamaPrePro tables for sampling.
	/// @details Setting this to "false" lets us use classic rama tables.  True by default.
	void set_use_rama_prepro_for_sampling( bool const setting );

	/// @brief Actually run the application.
	/// @details The initialize_from_options() function must be called before calling this.  (Called by default constructor.)
	void run() const;

private:
	/// ------------- Methods ----------------------------

	/// @brief Count the number of cis-peptide bonds in the pose.
	///
	core::Size count_cis_peptide_bonds(
		core::pose::PoseCOP pose
	) const;

	/// @brief Carry out the final FastRelax.
	///
	void
	do_final_fastrelax(
		core::pose::PoseOP pose,
		core::scoring::ScoreFunctionOP sfxn,
		core::Size const relax_rounds,
		bool const angle_min,
		bool const length_min,
		bool const cartesian_min
	) const;

	/// @brief Actually build the geometry that we'll be working with.
	///
	void build_polymer(
		core::pose::PoseOP pose,
		utility::vector1<std::string> const &restypes
	) const;

	/// @brief Add N-methylation.
	///
	void add_n_methylation(
		core::pose::PoseOP pose,
		core::Size const cyclic_offset
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
		core::scoring::ScoreFunctionOP sfxn_highhbond_cart,
		core::scoring::ScoreFunctionCOP sfxn_default,
		protocols::filters::CombinedFilterOP total_hbond,
		core::Size const cyclic_offset
	) const;

	/// @brief Set up the TaskOperations that conrol the design process, given user inputs.
	/// @details Default behaviour is designing all positions with L-canonicals and their
	/// D-equivalents EXCEPT cys and met (and gly), unless the user overrides this.
	void set_up_design_taskoperations( protocols::denovo_design::movers::FastDesignOP fdes, core::Size const cyclic_offset, core::Size const nres ) const;

	/// @brief Given a vector of full residue names of canonical residues, give me a concatenated list of one-letter codes.
	/// @details Does no checking for duplicates.
	std::string get_oneletter_codes( utility::vector1< std::string > const &fullnames ) const;

	/// @brief Given a vector of full residue names, give me a string of the form "NC <3-letter code> NC <3-letter code> NC <3-letter code> ..."
	/// @details Does no checking for duplicates.  Will fail gracelessly with invalid names.
	std::string get_nc_threeletter_codes( utility::vector1< std::string> const &fullnames ) const;

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

	/// @brief Given an position in the current (perturbed) pose, return the position in the original (unperturbed) pose.
	///
	inline core::Size original_position(
		core::Size const curr_position,
		core::Size const permutation_offset,
		core::Size const nresidue
	) const {
		core::Size returnval( curr_position + permutation_offset );
		if ( returnval > nresidue ) returnval -= nresidue;
		return returnval;
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

	/// @brief Are we sampling cis-prolines?
	/// @details Defaults to false.
	inline bool sample_cis_pro() const { return sample_cis_pro_; }

	/// @brief Frequency for sampling cis prolines.
	///
	inline core::Real const & sample_cis_pro_frequency() const { return sample_cis_pro_frequency_; }

	/// @brief Bond angle relax rounds.
	///
	inline core::Size angle_relax_rounds() const { return angle_relax_rounds_; }

	/// @brief Bond angle / bond length relax rounds.
	///
	inline core::Size angle_length_relax_rounds() const { return angle_length_relax_rounds_; }

	/// @brief Cartesian relax rounds.
	///
	inline core::Size cartesian_relax_rounds() const { return cartesian_relax_rounds_; }

	/// @brief Are we using RamaPrePro tables for sampling?
	/// @details True if we are, false if we're using classic rama tables.
	inline bool use_rama_prepro_for_sampling() const { return use_rama_prepro_for_sampling_; }

	/// @brief The length of the sequence, excluding crosslinkers.
	/// @details This is also the index of the last sequence residue.  Crosslinker residues follow this position.
	inline core::Size sequence_length() const { return sequence_length_; }

	/// @brief Given a pose with TBMB in it and another pose without TBMB, copy the TBMB residues from the first to the second,
	/// and add back covalent bonds.
	/// @details This function is called at the end of the protocol, and therefore doesn't bother to add back constraints.
	void re_append_tbmb_residues( core::pose::PoseCOP pose, core::pose::PoseOP newpose, core::Size const offset ) const;

private:
	/// ------------- Data -------------------------------
	/// -------- When you add new data to this class, ----
	/// -------- you must update the copy constructor ----

	/// @brief If this is called by MPI code, this can store the rank of the current process.  Zero otherwise.
	///
	int my_rank_;

	/// @brief The default ScoreFunction to use.  The high h-bond version is constructed from this.
	///
	core::scoring::ScoreFunctionOP scorefxn_;

	/// @brief Allows external code to suppress checkpointing, so that the SimpleCycpepPredictApplication doesn't write directly to disk.
	///
	bool suppress_checkpoints_;

	/// @brief Should this application produce silent file output?
	///
	bool silent_out_;

	/// @brief Should this application produce silent structure OP output, appending to a list?
	///
	bool silentlist_out_;

	/// @brief If silentlist_out_ is used, this is the list to append SilentStructurOPs to.
	/// @details Note that this is a non-const pointer to a vector.
	utility::vector1 < core::io::silent::SilentStructOP > * silentlist_;

	/// @brief If silentlist_out_ is used, this is the list of job summaries to which summaries should be appended.
	/// @details Note that this is a non-const pointer to a vector.
	utility::vector1 < SimpleCycpepPredictApplication_MPI_JobResultsSummaryOP > * summarylist_;

	/// @brief A native pose, provided by external code.
	/// @details If provided, this prevents the app from reading from the filesystem in the run() function.
	core::pose::PoseCOP native_pose_;

	/// @brief The prefix for the output filename.
	/// @details Defaults to "S_".
	std::string out_filename_;

	/// @brief The output score file name.
	/// @details Defaults to "default.sc".
	std::string out_scorefilename_;

	/// @brief Filename for the text file containing the sequence of the peptide.
	/// @details Must be provided with the -cyclic_peptide:sequence_file flag.
	std::string sequence_file_;

	/// @brief The string that would be read from a sequence file.
	/// @details If provided by external code, prevents filesystem read in run() function.
	std::string sequence_string_;

	/// @brief The length of the sequence, excluding any cross-linkers.
	/// @details Computed internally and stored for reference.  This is also the index of the last sequence residue (with cross-linkers
	/// following in linear sequence).
	mutable core::Size sequence_length_;

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

	/// @brief The total energy cutoff above which solutions are discarded.
	/// @details Defaults to 0.0.
	core::Real total_energy_cutoff_;

	/// @brief Determines whether the total energy cutoff should be used.
	/// @details Defaults to 'false'.
	bool use_total_energy_cutoff_;

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

	/// @brief Map of (seqpos -> phi/psi/omega triplet).
	/// @details The user can set certain alpha-amino acid mainchain dihedral values, if he or she so wishes.
	std::map< core::Size, utility::vector1 < core::Real > > user_set_alpha_dihedrals_;

	/// @brief A small random value added to all user-set dihedral values.  Defaults to 0.
	///
	core::Real user_set_dihedral_perturbation_;

	/// @brief Should we filter out solutions that have more than the allowed number of hydrogen bonds
	/// to an acceptor?
	/// @details Default true.
	bool filter_oversaturated_hbond_acceptors_;

	/// @brief The energy above which we no longer count a hydrogen bond to an acceptor for the filter.
	/// @details Default -0.1.
	core::Real oversaturated_hbond_cutoff_energy_;

	/// @brief At residues that are more likely to have cis omega values, sample cis some fraction of the time.
	/// @details Currently applies only to L-proline and D-proline.  Could be extended to beta-3-proline and to
	/// N-methylated amino acids, once those are added.  Read from database.
	bool sample_cis_pro_;

	/// @brief At residues that are more likely to have cis omega values, sample cis this fraction of the time.
	/// @details Read from database.
	core::Real sample_cis_pro_frequency_;

	/// @brief Should we design the peptide after sampling the conformation, or just relax it?
	/// @details Default "false" (relax only; no design).
	bool design_peptide_;

	/// @brief If we're designing, this is the filename for the file that tells us what residues are
	/// allowed at what postion.
	/// @details If left blank, no file is read.
	std::string design_filename_;

	/// @brief Used to prevent this process from directly reading the input file if it has been read by
	/// another process and the information has been transmitted to this process.
	/// @details Default false; set to true by set_allowed_residues_by_position().
	bool prevent_design_file_read_;

	/// @brief Allowed canonical residues at each position.
	/// @details Only used for design.  Map key "0" is used for default settings applied to
	/// any position not explicitly specified.
	std::map < core::Size, utility::vector1 < std::string > > allowed_canonicals_by_position_;

	/// @brief Allowed noncanonical residues at each position.
	/// @details Only used for design.  Map key "0" is used for default settings applied to
	/// any position not explicitly specified.
	std::map < core::Size, utility::vector1 < std::string > > allowed_noncanonicals_by_position_;

	/// @brief Should D-residues be prohibited at positions with negative phi values?
	/// @details Default true.
	bool prohibit_D_at_negative_phi_;

	/// @brief Should L-residues be prohibited at positions with positive phi values?
	/// @details Default true.
	bool prohibit_L_at_positive_phi_;

	/// @brief Should we use the aa_composition score term in design?
	/// @details Default false.
	bool use_aa_comp_;

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

	/// @brief Prevent counting of hydrogen bonds to adjacent residues?
	/// @details Default true.
	bool do_not_count_adjacent_res_hbonds_;

	/// @brief Number of rounds of relaxation with flexible bond angles.
	/// @details Default 0.
	core::Size angle_relax_rounds_;

	/// @brief Number of rounds of relaxation with flexible bond angles and bond lengths.
	/// @details Default 0.
	core::Size angle_length_relax_rounds_;

	/// @brief Number of rounds of Cartesian-space relaxation.
	/// @details Default 0.
	core::Size cartesian_relax_rounds_;

	/// @brief Should we use rama_prepro tables for sampling?
	/// @details Default true.
	bool use_rama_prepro_for_sampling_;

	/// @brief List of positions (in original sequence indexing -- not permuted) that are N-methylated.
	/// @details Defaults to empty list.
	utility::vector1< core::Size > n_methyl_positions_;

	/// @brief List of positions linked by 1,3,5-tris(bromomethyl)benzene.
	/// @details This is a vector of lists of three residues.
	utility::vector1< utility::vector1 < core::Size >  > tbmb_positions_;

	/// @brief If true, filters are applied based on distance between TBMB cysteines and on constraints to discard
	/// GenKIC solutions that can't be crosslinked easily.
	/// @details True by default.
	bool use_tbmb_filters_;

	/// @brief Multiplier to make the TBMB distance filter more permissive.
	/// @details Default 1.0.
	core::Real tbmb_sidechain_distance_filter_multiplier_;

	/// @brief Multiplier to make the TBMB constraints energy filter more permissive.
	/// @details Default 1.0.
	core::Real tbmb_constraints_energy_filter_multiplier_;

};

} //cyclic_peptide
} //protocols

#endif //INCLUDED_protocols_cyclic_peptide_predict_SimpleCycpepPredictApplication_hh
