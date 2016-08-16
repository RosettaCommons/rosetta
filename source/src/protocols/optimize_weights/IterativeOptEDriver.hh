// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/optimize_weights/IterativeOptEDriver.hh
/// @brief  A class for optimizing weights in an iterative weight fitting, sequence recovery test protocol
/// @author Andrew Leaver-Fay -- emulating a protocol by Jim Havranek and Brian Kuhlman.


#ifndef INCLUDED_protocols_optimize_weights_IterativeOptEDriver_hh
#define INCLUDED_protocols_optimize_weights_IterativeOptEDriver_hh


// MPI Headers have to be #included first
#ifdef USEMPI
#include <mpi.h>
#endif

// Unit headers
#include <protocols/optimize_weights/IterativeOptEDriver.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/optimization/types.hh>
#ifdef WIN32
#include <core/pose/Pose.hh>
#endif
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

#include <protocols/optimize_weights/OptEData.fwd.hh>

// Utility headers

/// STL Headers
//#include <map>

//Auto Headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/optimize_weights/OptEMultifunc.fwd.hh>


namespace protocols {
namespace optimize_weights {

enum AddStatus
{
	ADDED_STRUCTURE_AS_INTENDED,
	ADDED_STRUCTURE_OPPOSITE_AS_INTENDED,
	DID_NOT_ADD_STRUCTURE
};

class IterativeOptEDriver {
public:
	typedef core::Real                          Real;
	typedef core::Size                          Size;
	typedef core::scoring::ScoreTypes           ScoreTypes;
	typedef core::scoring::EnergyMap            EnergyMap;

public:
	IterativeOptEDriver();
	~IterativeOptEDriver();

	/// @brief take care when using a custom TaskFactory:
	/// TaskOperations must not 'accumulate state' as they will be reused repeatedly
	void task_factory( core::pack::task::TaskFactoryCOP );

	void go();

private:
	Size num_outer_iterations() const;
	Size num_inner_iterations() const;

	void read_tagfile_to_taskfactory(std::string tagfile_name, core::pack::task::TaskFactoryOP task_factory);
	void load_pose( core::pose::Pose & pose, std::string const & filename, bool ignore_centroid_input_flag );

	void divide_up_pdbs();
	void intialize_free_and_fixed_energy_terms();
	void setup_derived_free_and_fixed_data();

	void collect_rotamer_energies();
	void compute_rotamer_energies_for_assigned_pdbs();
	void collect_rotamer_energies_from_slave_cpus();
	void collect_rotamer_energies_from_slave_cpu( Size const which_cpu );
	void send_rotamer_energies_to_master_cpu();

	void optimize_weights();
	void score_position_data();
	utility::vector1< Real > read_reference_energies_from_file( std::string const & fname ) const;
	void write_new_scorefile();

	void test_sequence_recovery();
	void collect_sequence_recovery_data_from_slave_cpus();
	void collect_sequence_recovery_data_from_slave_cpu( Size const which_cpu );
	void run_design_on_assigned_pdbs();
	void send_sequence_recovery_data_to_master_cpu();

	void repack_assigned_pdbs();
	void collect_rotamer_recovery_data_from_slave_cpus();
	void collect_rotamer_recovery_data_from_slave_cpu( Size const which_cpu );
	void send_rotamer_recovery_data_to_master_cpu();

	bool decide_if_sequence_recovery_improved();

	void test_weight_sensitivity(
		OptEMultifunc const & func,
		core::optimization::Multivec const & dofs
	) const;

	std::string get_scorefile_name();
	void barrier();
	void exit_gracefully();

	utility::vector1< std::string >
	get_native_pdb_names();

	void zero_aa_counts();

	Real
	measure_sequence_recovery(
		utility::vector1< std::string > const & native_pdb_names,
		utility::vector1< std::string > const & names_for_output_pdbs,
		core::scoring::ScoreFunctionOP sfxn,
		Size & npos,
		Size & nrecovered
	);

	Real
	measure_rotamer_recovery(
		utility::vector1< std::string > const & native_pdb_names,
		utility::vector1< std::string > const & names_for_output_pdbs,
		core::scoring::ScoreFunctionOP sfxn,
		Size & npos,
		Size & nrecovered
	);

	Real
	opte_weight_mixing_factor(
		Size outer_loop_counter,
		Size inner_loop_counter
	);


	void
	initialize_free_and_fixed(
		core::scoring::EnergyMap & free_parameters,
		core::scoring::EnergyMap & fixed_parameters
	);

	bool
	converged(
		core::scoring::EnergyMap & free_parameters_prev,
		core::scoring::EnergyMap & free_parameters_curr,
		utility::vector1< Real > const & reference_energies_prev,
		utility::vector1< Real > const & reference_energies_curr
	);


	void
	setup_pdbnames_next_round(
		Size const outer_loop_counter,
		utility::vector1< std::string  > & pdbs_next_round,
		utility::vector1< std::string > const & native_pdb_names
	);


	void
	write_parameters_to_std_out(
		core::scoring::EnergyMap & free_parameters,
		utility::vector1< Real > const & reference_energies
	);

	/// @brief refactor this into an optE common...
	void
	get_nat_aa_opte_data(
		std::string const & pdb_name,
		core::pose::Pose & pose,
		core::pose::Pose & native_pose,
		core::scoring::ScoreFunction const & scorefxn,
		ScoreTypes & score_list,
		ScoreTypes & fixed_score_vec,
		OptEData & opte_data
	);

	void
	get_nat_rot_opte_data(
		std::string const & pdb_name,
		core::pose::Pose & pose,
		core::pose::Pose & native_pose,
		utility::vector1<bool> include_rsd,
		core::scoring::ScoreFunction const & scorefxn,
		ScoreTypes & score_list,
		ScoreTypes & fixed_score_vec,
		OptEData & opte_data
	);

	void
	set_aa_periodicity(
		PNatRotOptEPositionDataOP pos_data,
		core::chemical::AA aa
	) const;

	/// @brief True/False: a particular residue be skipped for inclusion
	/// into the native-rotamer recovery test
	/// because its dunbrack energy is really high
	bool
	residue_has_unacceptably_bad_dunbrack_energy(
		core::pose::Pose const & pose,
		Size const resid
	) const;

	/// @brief True/False: a particular residue be skipped for inclusion
	/// into the native-rotamer recovery test
	/// because it contains an atom with a high B-factor
	bool
	residue_has_bad_bfactor(
		core::pose::Pose const & pose,
		Size const resid
	) const;


	core::scoring::EnergyMap
	free_terms_energy_map_from_dofs(
		core::optimization::Multivec const & dofs
	) const;

	void
	free_weights_and_refEs_from_vars(
		utility::vector1< Real > const & vars,
		core::scoring::EnergyMap & weights,
		utility::vector1< Real > & reference_energies
	) const;

	core::scoring::ScoreFunctionOP
	configure_new_scorefunction() const;

	core::scoring::ScoreFunctionOP
	create_unweighted_scorefunction() const;

	core::scoring::ScoreFunctionOP
	create_weighted_scorefunction() const;

	void
	collect_decoy_discrimination_data();

	SingleStructureDataOP
	single_structure_data_for_pose(
		core::scoring::ScoreFunctionOP scorefxn,
		core::pose::Pose & pose,
		core::pose::Pose const & crystal_native,
		utility::vector1< Real > & free_data, // scratch space; avoids new
		utility::vector1< Real > & fixed_data, // scratch space; avoids new
		std::string const & structure_tag = ""
	) const;

	AddStatus
	add_structure_based_on_rms(
		SingleStructureDataOP ssd,
		PNatStructureOptEDataOP structure_data,
		bool intended_native
	) const;

	void
	compute_rotamers_around_ligands();

	void
	collect_ligand_discrimination_data();

	core::scoring::EnergyMap
	score_ligand_interface(
		core::scoring::ScoreFunction const & scorefxn,
		core::pose::Pose & pose // scoring is a non-const operation
	);

	SingleStructureDataOP
	make_simple_ssd_from_pdb(
		std::string const & pdb_filename,
		core::scoring::ScoreFunctionOP sfxn,
		bool pretend_no_fa_rep = false
	) const;

	void
	collect_dG_of_binding_data();

	void
	collect_ddG_of_mutation_data();

	void
	collect_ddG_of_binding_data();

	void
	load_pssm_data( std::string const & native_filename, Size const which_protein );

	core::pack::task::PackerTaskOP
	copy_native_packertask_logic(core::pose::Pose native_pose,
		core::pose::Pose context_pose,
		core::pack::task::TaskFactoryOP native_taskfactory);

#ifdef USEMPI
public:
	static
	void
	send_string_to_node( int destination, std::string const & string_to_send );

	static
	std::string
	receive_string_from_node( int source );
private:
#endif

	void
	repack_and_minimize_pose( core::pose::Pose & pose, core::scoring::ScoreFunctionOP sfxn ) const;

	std::string
	node_name( int rank );

	void print_energies(
		core::pose::Pose &,
		core::scoring::ScoreFunctionOP,
		std::ostream & os = std::cout
	);

	void output_weighted_unfolded_energies();

private:
	core::pack::task::TaskFactoryOP task_factory_;
	utility::vector1< std::string      > native_pdbs_;
	utility::vector1< core::pose::Pose > native_poses_; // the pose of the xtal native used in the design step
	utility::vector1< core::pose::Pose > context_poses_; // After the first round, the poses from the previous round of design
	utility::vector1< std::string      > pdbs_this_round_;
	utility::vector1< std::string      > next_iteration_pdbs_;

	utility::vector1< core::pose::Pose > rotamer_recovery_context_poses_;

	/// For rotamer recovery around ligands
	utility::vector1< std::string      > ligand_repack_pdbs_;
	utility::vector1< core::pose::Pose > ligand_repack_native_poses_;
	//utility::vector1< core::pose::Pose > ligand_repack_context_poses_;

	/// A list of file-list pairs.
	/// first = list of relaxed-native pdbs; second = list of decoy pdbs.
	utility::vector1< std::pair< std::string, std::string > >  decdisc_native_decoy_pairs_;
	utility::vector1< std::string > decdisc_crystal_natives_;
	utility::vector1< utility::vector1< core::pose::Pose > > decdisc_native_poses_;
	utility::vector1< utility::vector1< core::pose::Pose > > decdisc_decoy_poses_;
	utility::vector1< core::pose::Pose >decdisc_xtal_natives_;

	/// A list of file-list pairs.
	/// first = list of relaxed-native pdbs; second = list of decoy pdbs.
	utility::vector1< std::pair< std::string, std::string > >  ligand_native_decoy_pairs_;
	utility::vector1< std::string > ligand_crystal_natives_;

	/// A list of file-list pairs:
	/// first = list of relaxed mutant pdbs; second = list of native pdbs
	utility::vector1< std::pair< std::string, std::string > > ddg_mut_wt_pairs_;
	utility::vector1< Real > ddGs_;

	/// A list of lists for use in optimizing ddG binding predictions
	utility::vector1< utility::vector1< std::string > > ddG_bind_files_;
	utility::vector1< Real > ddGs_binding_;

	/// A list of PDB file pairs:
	/// first = "bound" structure; second = "unbound" structure with partners well separated
	utility::vector1< std::pair< std::string, std::string > > dG_bound_unbound_pairs_;
	utility::vector1< Real > dG_binding_;

	/// preserve the decoy scores by component between iterations to minimize I/O.
	OptEDataOP decoy_discrim_data_;
	/// preserve the decoy scores by component between iterations to minimize I/O.
	OptEDataOP ligand_discrim_data_;
	/// preserve the dGbinding data between iterations to minimize I/O
	OptEDataOP dG_binding_data_;
	/// preserve the ddMutation data between iterations to minimize I/O
	OptEDataOP ddG_mutation_data_;
	/// preserve the ddG bind mutation data between iterations to minimize I/O
	OptEDataOP ddG_bind_optE_data_;

	/// These two get set by the user
	EnergyMap free_parameters_;
	EnergyMap fixed_parameters_;

	/// These get updated from the above two
	EnergyMap include_terms_;
	Size include_count_;
	Size fixed_count_;
	Size free_count_;
	ScoreTypes free_score_list_;
	ScoreTypes fixed_score_list_;
	utility::vector1< Real > component_weights_;

	/// store the PSSM data so it's only read once.
	utility::vector1< utility::vector1< std::pair< core::chemical::AA, utility::vector1< Real > > > > all_pssm_data_;
	utility::vector1< std::pair< core::chemical::AA, utility::vector1< Real > > > pssm_data_;

#ifdef USEMPI
	MPI_Status stat_;
	int tag_;
#endif

	Size MPI_rank_;
	Size MPI_nprocs_;

	Size outer_loop_counter_;
	Size inner_loop_counter_;
	OptEDataOP optE_data_;

	Size total_positions_;
	Size count_recovered_;
	utility::vector1< Size > aa_obs_; /// the counts for each amino acid from the previous round of design (observed)
	utility::vector1< Size > aa_exp_; /// the counts for each amino acid in the input data set (expected)
	utility::vector1< Real > aa_freq_obs_; /// the frequency for each amino acid from the previous round of design
	utility::vector1< Real > aa_freq_exp_; /// the frequency for each amino acid in the input data set

	Size total_rotamer_positions_;
	Size count_rotamers_recovered_;

	/// Hold the result of minimization -- interpolate between the
	/// before and after weight sets during iterative sequence recovery.
	utility::vector1< Real > before_minimization_reference_energies_;
	utility::vector1< Real > after_minimization_reference_energies_;
	utility::vector1< Real > reference_energies_inner_loop_;

	/// When using the optEMultifunc wrapper, hold on to the dofs that
	/// the minimizer has access to and interpolate in the space of
	/// those dofs; then convert the interpolated dofs into weights.
	/// Interpolating the weights directly might result in weights that
	/// violate the mixing rules for the wrapper.
	WrapperOptEMultifuncOP   wrapped_opt_min_;
	utility::vector1< Real > minimizer_dofs_before_minimization_;
	utility::vector1< Real > minimizer_dofs_after_minimization_;
	utility::vector1< Real > minimizer_dofs_mixed_;

	EnergyMap free_weights_before_minimization_;
	EnergyMap free_weights_after_minimization_;
	EnergyMap free_weights_inner_loop_;

	Real mixing_factor_;

	Real outer_loop_last_sequence_recovery_rate_;
	Real outer_loop_seq_profile_cross_entropy_;
	Real inner_loop_sequence_recovery_rate_;
	Real inner_loop_rotamer_recovery_rate_; // output for curriousity's sake.

	bool using_unfolded_energy_term_;

};

/// @brief Read options[ optE::component_weights ] file from input file.
/// Not a member of the above driver class since its independent of the driver;
/// possibly belongs in a separate source file.
/// Any component specified in the weights file is set to the corresponding weight.
/// Any component not specified in the weights file is set to 1.
void
load_component_weights(
	utility::vector1< core::Real > & component_weights
);

} // namespace optimize_weights
} // namespace protocols

#endif
