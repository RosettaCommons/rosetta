// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/helical_bundle_predict/HelicalBundlePredictApplication.hh
/// @brief The meat-and-potatoes for the helical_bundle_predict application, used to predict structures of helical bundles
/// made from canonical or noncanonical building-blocks.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_helical_bundle_predict_HelicalBundlePredictApplication_hh
#define INCLUDED_protocols_helical_bundle_predict_HelicalBundlePredictApplication_hh

#include <protocols/helical_bundle_predict/HelicalBundlePredictApplication.fwd.hh>

// Protocols headers
#include <protocols/helical_bundle_predict/HBPHelixAssignments.hh>
#include <protocols/helical_bundle_predict/HBP_TemperatureScheduleGenerator.fwd.hh>
#include <protocols/helical_bundle_predict/HBP_MoveGenerator.fwd.hh>
#include <protocols/cyclic_peptide_predict/HierarchicalHybridJD_JobResultsSummary.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/import_pose/import_pose.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace protocols {
namespace helical_bundle_predict {

/// @brief Options for the application.
/// @details Prevents repeated calls to the global options system.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org).
class HelicalBundlePredictApplicationOptions : public utility::VirtualBase {

public:

	/// @brief Constructor
	/// @details Triggers read from options system!
	HelicalBundlePredictApplicationOptions();

	/// @brief Destructor
	~HelicalBundlePredictApplicationOptions() override;

	/// @brief Clone operator: copy this object and return a smart pointer to the copy.
	HelicalBundlePredictApplicationOptionsOP clone() const;

public: //Member functions

	/// @brief Indicate which commandline flags are relevant (i.e. which should be listed with the --help flag).
	/// @details This is a static function that must be called BEFORE devel_init().
	static void register_options();

	/// @brief Read from the options system to initialize this object.
	void initialize_from_options();

	/// @brief Set the file containing the FASTA sequence.
	void set_fasta_file( std::string const & file_in );

	/// @brief Set the file containing the full-basename sequence.
	void set_sequence_file( std::string const & file_in );

	/// @brief Set the file containing the helix assignments.
	void set_helix_assignment_file( std::string const & file_in );

	/// @brief Set the contents of the FASTA file.
	void set_fasta_file_contents( std::string const & contents_in );

	/// @brief Set the contents of the seqeunce file.
	void set_sequence_file_contents( std::string const & contents_in );

	/// @brief Set the contents of the helix assignment file.
	void set_helix_assignment_file_contents( std::string const & contents_in );

	/// @brief Get the residues to ignore in the native pose when setting up the alignment for RMSD.
	/// @details Throws errors if any are zero or negative.
	void set_rmsd_residues_to_ignore_native( utility::vector1< signed long > const &input );

	/// @brief Get the residues to ignore in the generated poses when setting up the alignment for RMSD.
	/// @details Throws errors if any are zero or negative.
	void set_rmsd_residues_to_ignore_prediction( utility::vector1< signed long > const &input );

	/// @brief Get the file containing the sequence.
	inline std::string const & fasta_file() const { return fasta_file_; }

	/// @brief Get the file containing the helix assignments.
	inline std::string const & helix_assignment_file() const { return helix_assignment_file_; }

	/// @brief Get the contents of the FASTA file.
	inline std::string const & fasta_file_contents() const { return fasta_file_contents_; }

	/// @brief Get the contents of the sequence file.
	inline std::string const & sequence_file_contents() const { return sequence_file_contents_; }

	/// @brief Get the contents of the helix assignment file.
	inline std::string const & helix_assignment_file_contents() const { return helix_assignment_file_contents_; }

	/// @brief Given input filenames, read the files.
	/// @details INVOLVES READS FROM DISK!  WARNING!
	void read_inputs();

	/// @brief Get the number of simulated annealing rounds in centroid mode.
	inline core::Size num_simulated_annealing_rounds_centroid() const { return num_simulated_annealing_rounds_centroid_; }

	/// @brief Get the number of steps per simulated annealing round in centroid mode.
	inline core::Size num_steps_per_simulated_annealing_round_centroid() const { return num_steps_per_simulated_annealing_round_centroid_; }

	/// @brief Get the max temperature during simulated annealing steps in centroid mode.
	inline core::Real const & centroid_max_temperature() const { return centroid_max_temperature_; }

	/// @brief Get the min temperature during simulated annealing steps in centroid mode.
	inline core::Real const & centroid_min_temperature() const { return centroid_min_temperature_; }

	/// @brief Get the number of repeated applications of this protocol.
	inline core::Size nstruct() const { return nstruct_; }

	/// @brief Will be doing fullatom refinement?
	inline bool do_fullatom_refinement() const { return do_fullatom_refinement_; }

	/// @brief If we are doing fullatom refinement, how many rounds of FastRelax should we apply?
	inline core::Size fullatom_fast_relax_rounds() const { return fullatom_fast_relax_rounds_; }

	/// @brief If we are doing fullatom refinement, should we try disulfide permutations?
	inline bool fullatom_find_disulfides() const { return fullatom_find_disulfides_; }

	/// @brief Get the native filename.
	inline std::string const & native_file() const { return native_file_; }

	/// @brief Get the residues to ignore in the native pose when setting up the alignment for RMSD.
	inline utility::vector1< core::Size > const & rmsd_residues_to_ignore_native() const { return rmsd_residues_to_ignore_native_; }

	/// @brief Get the residues to ignore in the generated poses when setting up the alignment for RMSD.
	inline utility::vector1< core::Size > const & rmsd_residues_to_ignore_prediction() const { return rmsd_residues_to_ignore_prediction_; }

private: //Functions

	/// @brief Read a FASTA file from disk.
	void read_fasta();

	/// @brief Read a sequence file from disk.
	void read_sequence_file();

	/// @brief Given a set of characters, find the first instance of any of them in a string
	/// and return the (zero-based) index of that character.
	core::Size findchar( std::string const & curstring, utility::vector1< char > const & chars ) const;

	/// @brief Given FASTA file contents, remove comment lines.
	void clean_fasta_file_contents();

	/// @brief Read a helix assignemnt file from disk.
	void read_helix_assignments();

private: //Variables

	/// @brief The native structure.
	std::string native_file_;

	/// @brief The file containing the sequence (in FASTA format).
	std::string fasta_file_;

	/// @brief The file containing the sequence (in full name format).
	std::string sequence_file_;

	/// @brief The file containing the helix assignments.
	std::string helix_assignment_file_;

	/// @brief The contents of the FASTA file.
	std::string fasta_file_contents_;

	/// @brief The contents of the sequence file.
	std::string sequence_file_contents_;

	/// @brief The contents of the helix assignment file.
	std::string helix_assignment_file_contents_;

	/// @brief Number of simulated annealing rounds in centroid mode.
	core::Size num_simulated_annealing_rounds_centroid_=3;

	/// @brief Number of steps in each simulated annealing round in centroid mode.
	core::Size num_steps_per_simulated_annealing_round_centroid_=1000;

	/// @brief The max temperature during simulated annealing steps in centroid mode.
	core::Real centroid_max_temperature_=50;

	/// @brief The min temperature during simulated annealing steps in centroid mode.
	core::Real centroid_min_temperature_=0.62;

	/// @brief The number of attempts to make.
	core::Size nstruct_=1;

	/// @brief Will be doing fullatom refinement?
	bool do_fullatom_refinement_=true;

	/// @brief If we are doing fullatom refinement, how many rounds of FastRelax should we apply?
	core::Size fullatom_fast_relax_rounds_=3;

	/// @brief If we are doing fullatom refinement, should we try disulfide permutations?
	bool fullatom_find_disulfides_=true;

	/// @brief Residues to ignore in the native pose when setting up the alignment for RMSD.
	utility::vector1< core::Size > rmsd_residues_to_ignore_native_;

	/// @brief Residues to ignore in the generated poses when setting up the alignment for RMSD.
	utility::vector1< core::Size > rmsd_residues_to_ignore_prediction_;

};

/// @brief The meat-and-potatoes for the helical_bundle_predict application, used to predict structures of helical bundles made from canonical or noncanonical building-blocks.
class HelicalBundlePredictApplication : public utility::VirtualBase {

public:

	/// @brief Delete the default constructor, and require an options constructor.
	HelicalBundlePredictApplication() = delete;

	/// @brief Options constructor.
	/// @details Stores input options directly; doesn't clone.
	/// @note Triggers read from disk to set up move generator!
	HelicalBundlePredictApplication( HelicalBundlePredictApplicationOptionsCOP options_in );

	/// @brief Options + move generator constructor.
	/// @details Stores input options directly; doesn't clone.
	/// @note Avoids read from disk by using a move generator that was already set up.  The input move
	/// generator and the input scorefunctions are used directly (not cloned).
	HelicalBundlePredictApplication( HelicalBundlePredictApplicationOptionsCOP options_in, HBP_MoveGeneratorOP centroid_move_generator_in, core::scoring::ScoreFunctionOP centroid_sfxn_in, core::scoring::ScoreFunctionOP fullatom_sfxn_in );

	/// @brief Copy constructor.
	HelicalBundlePredictApplication(HelicalBundlePredictApplication const & src);

	/// @brief Destructor.
	~HelicalBundlePredictApplication() override;

	HelicalBundlePredictApplicationOP
	clone() const;

public: //Member functions:

	/// @brief Actually run the application and produce output.
	void run();

	/// @brief Set the native pose.
	/// @details Does not clone the input; sets owning pointer directly.
	void set_native( core::pose::PoseCOP native );

	/// @brief Given a pose, align it to the native pose.
	/// @details Throws an error if there's a mismatch between the pose lengths or mainchain atom counts.
	/// @returns RMSD to native.
	core::Real align_to_native_pose( core::pose::Pose & pose ) const;

	/// @brief Create a new move generator.  Static, so this can be called from other classes (e.g.
	/// HelicalBundlePredictApplication_MPI.)
	/// @details Triggers read from disk!
	static HBP_MoveGeneratorOP create_centroid_move_generator();

	/// @brief Create the move generator used for the final fullatom refinement step.
	/// @details Not static, since it depends on the options_ object.
	/// @note In its current form, this should not trigger a read from disk.
	HBP_MoveGeneratorOP create_final_fullatom_refinement_move_generator() const;

	/// @brief Create the scorefunction used during centroid mode.
	/// @details Reads from disk!  Do not use repeatedly!  Store the result rather than regenerating it!
	static core::scoring::ScoreFunctionOP create_centroid_scorefunction();

	/// @brief Create the scorefunction used during fullatom mode.
	/// @details Reads from disk!  Do not use repeatedly!  Store the result rather than regenerating it!
	/// @note This returns whatever the current default fullatom scorefunction is, currently.  This function is here to make it easy
	/// to hard-code a specialized scorefunction in the future if necessary.
	static core::scoring::ScoreFunctionOP create_fullatom_scorefunction();

	/// @brief Set the number of repeats to try.
	inline void set_nstruct( core::Size const nstruct_in ) { nstruct_ = nstruct_in; }

#ifdef USEMPI
	/// @brief Set the job summary list and full structure list to which this protocol should write output.
	/// @details This is a rare instance in which using raw pointers is appropriate.  DO NOT EMULATE THIS UNLESS YOU KNOW WHAT YOU'RE DOING.
	void
	set_output(
		utility::vector1 < protocols::cyclic_peptide_predict::HierarchicalHybridJD_JobResultsSummaryOP > * jobsummaries,
		utility::vector1 < core::io::silent::SilentStructOP > * all_output
	);

	/// @brief Write output, not to disk, but to a list of binary silent structures in memory.
	/// @details This also writes a summary of the job (basically, energy + RMSD + jobindex) to a job summary list.
	void
	output_to_silent_list(
		core::pose::Pose const &pose,
		utility::vector1 < protocols::cyclic_peptide_predict::HierarchicalHybridJD_JobResultsSummaryOP > & jobsummaries,
		utility::vector1 < core::io::silent::SilentStructOP > & all_output,
		core::Size const repeat_index,
		bool const include_rmsd,
		core::Real const & rmsd
	) const;

	/// @brief Set the rank of this MPI process.
	inline void set_my_rank( core::Size const rank_in ) { my_rank_ = rank_in; }

	/// @brief Set the number of jobs already done on this process or thread.
	inline void set_already_completed_job_count( core::Size const count_in ) { already_completed_job_count_ = count_in; }

#else //Output options for non-MPI build:

	/// @brief Set the output format.  (Automatically sets extension.)
	/// @details This one sets silent output to false.
	void set_output_format( core::import_pose::FileType const type );

	/// @brief Indicate that we're using silent output.  This overrides set_ouput_format(), and
	/// sets the extension automatically.
	void set_silent_output();

	/// @brief Set the prefix and suffix for output.
	void set_output_prefix_and_suffix( std::string const & prefix, std::string const & suffix );

#endif //USEMPI


private: //Member functions:

	/// @brief Read the native pose in from disk, based on the native_file_ set in the options_ object.
	/// @details READS FROM DISK!  This can be avoided by passing in a pointer to the native pose.
	void load_native_pose_from_disk();

	/// @brief Given the helix assignment file's contents, set up the HBP_HelixCoilMoveGenerator used in centroid mode..
	/// @details The helix_assignment_file_content_ variable must be populated first!
	void set_up_centroid_move_generator();

	/// @brief Construct a pose from the contents of a FASTA file.  The conformation is set to linear at this point.
	core::pose::PoseOP make_pose_from_fasta_contents() const;

	/// @brief Construct a pose from the contents of a sequence file.  The conformation is set to linear at this point.
	core::pose::PoseOP make_pose_from_sequence_file_contents() const;

	/// @brief Do a single simulated annealing round (a Monte Carlo trajectory with temperature ramping from high to low, once).
	void do_a_simulated_annealing_round(
		core::pose::PoseOP & pose,
		core::Size num_steps,
		HBP_MoveGeneratorCOP move_generator,
		HBP_TemperatureScheduleGeneratorCOP temperature_generator,
		core::scoring::ScoreFunction const &sfxn
	) const;

	/// @brief Entry point into the simulated annealing.
	void do_simulated_annealing(
		core::pose::PoseOP & pose,
		core::Size const num_rounds,
		core::Size const num_steps_per_round,
		HBP_MoveGeneratorOP move_generator,
		HBP_TemperatureScheduleGeneratorOP temperature_generator,
		core::scoring::ScoreFunction const &sfxn
	) const;

	/// @brief Carry out the final full-atom refinement steps.
	void do_final_fullatom_refinement( core::pose::PoseOP & pose ) const;

	/// @brief Given the old energy, the new energy, and the temperature, apply the Metropolis criterion.
	bool apply_metropolis_criterion( core::Real const &old_energy, core::Real const &new_energy, core::Real const &temperature ) const;

	/// @brief Check that the list of residues to ignore in calculating RMSD is reasonable.
	/// @details Throws if residues are outside the size of the posem or if residues are
	/// provided but there's no pose.
	void
	check_ignore_residues_reasonable(
		utility::vector1< core::Size > const & ignore_residues,
		core::pose::PoseCOP const & pose
	) const;

private: //Variables:

	/// @brief Options for the application.
	HelicalBundlePredictApplicationOptionsCOP options_;

	/// @brief The native pose.  Loaded from disk when run() is called unless provided separately.
	core::pose::PoseCOP native_pose_;

	/// @brief The pose that we're working on.
	core::pose::PoseOP pose_;

	/// @brief Move generator for centroid-mode Monte Carlo.
	HBP_MoveGeneratorOP centroid_move_generator_;

	/// @brief Temperature schedule generator for centroid-mode Monte Carlo.
	HBP_TemperatureScheduleGeneratorOP centroid_temperature_generator_;

	/// @brief The scorefunction to use in centroid mode.
	core::scoring::ScoreFunctionOP centroid_sfxn_;

	/// @brief The scoring function to use for fullatom refinement steps.
	core::scoring::ScoreFunctionOP fullatom_sfxn_;

	/// @brief Move generator for the full-atom refinement move at the end.
	/// @note This generates a move that is applied ONCE.
	HBP_MoveGeneratorOP final_fullatom_refinement_move_generator_;

	/// @brief Number of repeats to attempt.
	core::Size nstruct_;



#ifdef USEMPI
	//Variables used only in MPI mode:

	/// @brief Rank of this MPI process
	core::Size my_rank_;

	/// @brief The number of jobs already completed on this thread or process.
	core::Size already_completed_job_count_;

	/// @brief Pointer to the vector where I should store job summaries.
	/// @details This is a rare instance in which using raw pointers is appropriate.  DO NOT EMULATE THIS UNLESS YOU KNOW WHAT YOU'RE DOING.
	utility::vector1 < protocols::cyclic_peptide_predict::HierarchicalHybridJD_JobResultsSummaryOP > * jobsummaries_;

	/// @brief Pointer to the vector where I should store full silent structures of all output.
	/// @details This is a rare instance in which using raw pointers is appropriate.  DO NOT EMULATE THIS UNLESS YOU KNOW WHAT YOU'RE DOING.
	utility::vector1 < core::io::silent::SilentStructOP > * all_output_;

#else
	//Variables for output in non-MPI mode:

	/// @brief The output prefix.
	/// @details Not used in MPI build.
	std::string outfile_prefix_ = "result_";

	/// @brief The output suffix.
	/// @details Not used in MPI build.
	std::string outfile_suffix_ = "";

	/// @brief The output extension.
	/// @details Not used in MPI build.
	std::string outfile_extension_ = "pdb";

	/// @brief The output format.
	/// @details Not used in MPI build.
	core::import_pose::FileType output_filetype_ = core::import_pose::PDB_file;

	/// @brief Are we doing silent output?
	/// @details Not used in MPI build.
	bool silent_output_ = false;
#endif

};


} //protocols
} //helical_bundle_predict



#endif //INCLUDED_protocols_helical_bundle_predict_HelicalBundlePredictApplication_hh





