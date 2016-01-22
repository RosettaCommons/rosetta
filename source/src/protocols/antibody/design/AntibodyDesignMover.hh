// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyDesignMover.hh
/// @brief Class that initially designs antibodies through grafting using an AntibodyDatabase + North_AHO numbering scheme
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_design_AntibodyDesignMover_hh
#define INCLUDED_protocols_antibody_design_AntibodyDesignMover_hh

// Project Includes
#include <protocols/antibody/design/AntibodyDesignMover.fwd.hh>
#include <protocols/antibody/design/AntibodyDesignEnum.hh>
#include <protocols/antibody/design/AntibodyDesignEnumManager.fwd.hh>
#include <protocols/antibody/design/GeneralAntibodyModeler.hh>
#include <protocols/antibody/database/CDRSetOptions.hh>
#include <protocols/antibody/design/CDRGraftDesignOptions.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>

#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/design/MutateFrameworkForCluster.fwd.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>

#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/constraints/ParatopeSiteConstraintMover.fwd.hh>
#include <protocols/antibody/constraints/ParatopeEpitopeSiteConstraintMover.fwd.hh>
#include <protocols/antibody/constraints/CDRDihedralConstraintMover.fwd.hh>

// Protocol Includes
#include <protocols/grafting/CCDEndsGraftMover.fwd.hh>
#include <protocols/grafting/AnchoredGraftMover.fwd.hh>
#include <protocols/simple_moves/MinMover.hh>

// Core Includes
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/Constraint.hh>

// Protocol Includes
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/types.hh>

namespace protocols {
namespace antibody {
namespace design {

//Sampling weights not currently implemented
//struct SamplingWeights {
// std::map< Size,  Real > types; //Used to sample more of one type then others.  Sample most in type1, but use a bit of type2/3.
// Real native_cluster; //May be used to weight the sampling of native clusters vs others.
// Real cdr; //Used to sample more of one or particular CDRs rather then others (H3)
// Real center; //Used to sample more of cluster centers then other CDRs.
//};


/// @brief This class designs antibodies by grafting, using cluster information and constraints to help.
/// It represents the first step in The Rosetta Antibody Designer, but can be used outside of the application.
/// A 2011 Antibody Database is included in the rosetta datase.  Up-to-date versions can be downloaded from
///
///
/// @details To use this class:
/// 1) Use default instruction path, or set one.  Loads CDRs from AntibodyDatabase
///     *See rosetta_database/sampling/antibodies for file format.
/// 2) Use class interface settings to control sampling.  Loads CDRs from AntibodyDatabase.
/// 3) Use set_cdr_sets to set your own CDRs to graft.
///
class AntibodyDesignMover: public protocols::moves::Mover {

public:

	AntibodyDesignMover();

	AntibodyDesignMover( AntibodyInfoCOP ab_info );

	virtual ~AntibodyDesignMover();

	/// @brief Reads default CDRGraftInstruction file and creates structs.
	void
	set_defaults();

	/// @brief Parse my tag for RosettaScripts.  Main RS interface to full antibody design protocol.
	virtual void
	parse_my_tag(
		TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & filters,
		moves::Movers_map const & movers,
		Pose const & pose
	);

	virtual void
	apply(core::pose::Pose & pose);

	protocols::moves::MoverOP
	clone() const;

	protocols::moves::MoverOP
	fresh_instance() const;


public:

	/// @brief Set the design protocol.  Default is classic_monte_carlo.
	/// Available options are: generalized_monte_carlo, deterministc_graft.
	///
	/// @details. Deterministic Protocol only available if GraftDesigning a single CDR.
	void
	set_design_protocol( AntibodyDesignProtocolEnum design_protocol);

	/// @brief Set the options which will be used for querying the database
	void
	set_cdr_set_options(AntibodyCDRSetOptions cdr_set_options);

	/// @brief Set CDR-specific GraftDesign options
	void
	set_graft_design_options(AntibodyCDRGraftDesignOptions graft_design_options);

	/// @brief Set CDR-specific SeqDesign options
	void
	set_seq_design_options(AntibodyCDRSeqDesignOptions seq_design_options);

	/// @brief Will not initialize CDRs from the AntibodyDatabase.  Use if you have your own CDR's you are interested in grafting.
	/// @details Overhang residues will be used for superposition.  To identify your CDRs, use the functions in AntibodyInfo or the pilot app identify_cdr_clusters
	/// in pilot/jadolfbr/cluster_utils
	//void
	//set_cdr_set(CDRDBPoseSet & cdr_set, core::Size overhang);

	/// @brief Graft and Sequence design these CDRs.  Overrides all options settings.
	void
	set_cdr_override(utility::vector1<CDRNameEnum> cdrs_to_design);

	/// @brief Set the instruction file path to load options from instead of setting these via associated options functions
	/// Or loading via cmd-line option.  Gotta love supporting 3 different Rosetta interfaces.
	void
	set_instruction_file(std::string instruction_file);

public:
	////////////////////////////////////////////////////////////////////////////
	// Modeling Settings
	//
	//

	void
	set_scorefunction(core::scoring::ScoreFunctionOP scorefxn);

	void
	set_scorefunction_min(core::scoring::ScoreFunctionOP min_scorefxn);

	/// @brief Set the number of overall rounds.  Each round will choose a CDR and design it and others based on options
	void
	set_outer_cycles(core::Size outer_cycles);

	/// @brief Set the number of inner cycles, which optimizes the structure and sequence of the antibody after any graft.
	/// Essentially min cycles.  Recommend 3 cycles if not doing relax as the mintype.
	void
	set_inner_cycles(core::Size inner_cycles);


	/// @brief Set the algorithm to run a low-resolution docking step after each graft.  Default false.
	/// @details Uses command-line options for docking.  Lower inner or outer cycles if you are grafting many CDRs.  Still quicker then relax.
	void
	set_dock_post_graft(bool dock_post_graft);

	void
	set_dock_rounds(core::Size dock_rounds);

	/// @brief Set the algorithm to run a final rigid-body minimization of antigen/antibody interface post graft.  Useful if not docking - quicker, but less reliable then full dock.
	void
	set_rb_min_post_graft(bool rb_min_post_graft);


	/// @brief Sets the protocol to keep a specific number of top designs.  Default is 10
	void
	set_keep_top_designs(core::Size top_designs);

	/// @brief Get the top designs found.  You can then use them in other protocols, dump them, etc. They are in order.
	/// @details - This should be refactored to get_additional_output.
	utility::vector1< core::pose::PoseOP>
	get_top_designs(){
		return top_designs_;
	};


public:
	////////////////////////////////////////////////////////////////////////////
	// Paratope and Epitope constraints
	//
	//


	/// @brief Set any paratope CDRs.  If not set, will use all CDRs as the paratope where needed.
	/// Used mainly for constraints
	void
	set_paratope_cdrs(utility::vector1<bool> const & cdrs);

	/// @brief Set any epitope residues in PDB numbering.  If not set, they will be detected automatically via the interface.
	void
	set_epitope_residues(utility::vector1<PDBNumbering > epitope_residues);

	/// @brief Setting to use epitope constraints.  Without this false, will not use any set epitope residues.
	void
	set_use_epitope_constraints(bool use_epitope_csts);

public:
	void set_interface_dis(core::Real interface_dis);
	void set_neighbor_dis(core::Real neighbor_dis);

public:

	virtual std::string get_name() const;

	virtual void
	show(std::ostream & output=std::cout) const;

private:

	void
	read_command_line_options();

	/// @brief Setup ALL options classes
	void
	setup_options_classes();

	//Determine each CDRs cluster if not already calculated and set in AntibodyInfo instance.
	void
	setup_native_clusters(core::pose::Pose & pose);

	void
	setup_native_sequence(core::pose::Pose & pose);

	void
	setup_epitope_residues(core::pose::Pose const & pose);

	void
	setup_paratope_epitope_constraints(core::pose::Pose & pose);

	void
	setup_cart_minimizer();

	void
	setup_modeler();

	void
	finalize_setup(core::pose::Pose & pose);

	void
	setup_scorefxn();

	/// @brief Extend the native CDRs to be designed for benchmarking.
	void
	setup_random_start_pose(core::pose::Pose & pose, vector1<CDRNameEnum> & cdrs_to_design);

	void
	setup_default_graft_settings();

	/// @brief Setup alternative cdr pose indexes and sampling strategies.
	void
	setup_cdr_pose_sampling_strategies();


	/// @brief Uses instructions to Query the AntibodyDatabase and load poses.
	void
	initialize_cdr_set(core::pose::Pose const & pose);

	void
	check_for_top_designs(core::pose::Pose & pose);

	/// @brief Gets a list of vectors whose indexes correspond to CDRNameEnum, and whose values correspond to the cdr_set index.  If the value is 0, it means no cdr in set.
	utility::vector1< utility::vector1< core::Size > >
	get_cdr_set_index_list();


	/// @details Applies graft, modeling, etc. to a single CDR.
	bool
	apply_to_cdr(core::pose::Pose & pose, CDRNameEnum cdr, core::Size index, bool min_post_graft = true);

	/// @brief Grafts a single CDR into framework. Index is the vector index of CDRSet/CDRClusterMap.
	/// @details Return success or failure
	bool
	graft_in_cdr(core::pose::Pose & pose, CDRNameEnum const cdr, CDRDBPose & cdr_pose);

	//Runs graft then cartmin on AnchoredGraftMover or CCDEndsGraftMover.  Returns boolean of closure.
	std::pair<bool, core::Size>
	run_graft(core::pose::Pose & pose, CDRNameEnum const cdr, CDRDBPose & cdr_pose, grafting::AnchoredGraftMoverOP grafter);

	/// @brief Run the inner structure/sequence optimization cycle
	void
	run_optimization_cycle(core::pose::Pose & pose, protocols::moves::MonteCarlo & mc, CDRNameEnum cdr);

	/// @brief Mutates framework residues needed to stay within a particular cluster.  Only one  (L1-11-1) is known to have framework dependencies.  For now.
	/// Will be replaced by AntibodyDesignOptimizer
	//void
	//mutate_framework_residues(pose::Pose pose, CDRClusterEnum cluster);

private:
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Main Algorithms
	//
	//

	/// @brief.  If rounds >= number of possible combinations - Try them all.
	/// @details Grafts CDRs on the input structure, in a random order.
	void
	run_deterministic_graft_algorithm(core::pose::Pose & pose, utility::vector1<CDRNameEnum>& cdrs_to_design);

	/// @brief Basic mc algorithm that randomly samples from the cdr set.
	void
	run_basic_mc_algorithm(core::pose::Pose & pose, utility::vector1<CDRNameEnum>& cdrs_to_design, AntibodyDesignProtocolEnum mc_algorithm);



	void
	print_str_vec(std::string const name, utility::vector1<std::string> const & vec, std::ostream & output=std::cout) const;


private:

	AntibodyInfoOP ab_info_;
	AntibodyDesignEnumManagerOP design_enum_manager_;

	AntibodyCDRSetOptions cdr_set_options_;
	AntibodyCDRGraftDesignOptions cdr_graft_design_options_;
	AntibodyCDRSeqDesignOptions cdr_seq_design_options_;
	AntibodySeqDesignTFCreatorOP seq_design_creator_;

	std::map< CDRNameEnum, utility::vector1< CDRDBPose > > cdr_set_;

	//CDRSet Indexing for CDR sampling strategies:

	/// @brief Even cluster sampling (even_cluster_mc protocol)
	std::map< CDRNameEnum,
		std::map< clusters::CDRClusterEnum,
		utility::vector1< core::Size > > > cluster_based_CDRDBPose_indexes_;

	/// @brief Even Length and Cluster sampling (even_length_cluster_mc protocol)
	std::map< CDRNameEnum,
		std::map< core::Size,
		std::map< clusters::CDRClusterEnum,
		utility::vector1< core::Size > > > > length_based_CDRDBPose_indexes_;


	protocols::grafting::CCDEndsGraftMoverOP graft_mover_;
	protocols::grafting::AnchoredGraftMoverOP anchored_graft_mover_;

	MutateFrameworkForClusterOP framework_mutator_;

	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreFunctionOP scorefxn_min_;
	core::scoring::ScoreFunctionOP scorefxn_cart_graft_;

	GeneralAntibodyModelerOP modeler_;
	protocols::simple_moves::MinMoverOP cart_min_graft_;

	protocols::moves::MonteCarloOP mc_;

	constraints::ParatopeEpitopeSiteConstraintMoverOP paratope_epitope_cst_mover_;
	constraints::ParatopeSiteConstraintMoverOP paratope_cst_mover_;
	constraints::CDRDihedralConstraintMoverOP cdr_dihedral_cst_mover_;

	core::Size overhang_;
	core::Size outer_cycles_;
	core::Size inner_cycles_;
	core::Size dock_cycles_;
	core::Size num_top_designs_; //Number of top designs to keep.

	core::Real interface_dis_;
	core::Real neighbor_dis_;

	//Can be a struct.  Just not now.
	utility::vector1< core::pose::PoseOP > top_designs_;
	utility::vector1< core::Real> top_scores_;

	core::Size total_permutations_; //Total number of possible combinations

	//Overall Booleans
	bool dock_post_graft_; //Run a low and high resolution docking step after the graft?
	bool rb_min_post_graft_; //Run an rb_min step post graft.  Useful if not using docking.
	//bool design_neighbor_cdrs_;

	//Any cmd-line set Paratope and Epitope residues
	utility::vector1<bool> paratope_cdrs_;
	utility::vector1<PDBNumbering > epitope_residues_; //Vector of resnum, chain pairs.

	bool adapt_graft_;
	bool enable_adapt_graft_cartesian_;
	bool benchmark_;
	bool use_light_chain_type_;
	bool use_epitope_constraints_;
	bool print_tracer_info_; //Quick change to not print tracer info on  random starting cdr



	bool idealize_graft_cdrs_; //Idealize CDRs during graft?

	bool add_log_to_pose_;
	utility::vector1<std::string> graft_log_;
	utility::vector1<std::string> accept_log_;

	//Monte Carlo KT
	core::Real outer_kt_;
	core::Real inner_kt_;

	bool enable_full_protocol_atom_pair_cst_;
	AntibodyDesignProtocolEnum design_protocol_;

	utility::vector1<CDRNameEnum> design_override_;
	utility::vector1<CDRNameEnum> cdrs_to_design_;
	std::string instruction_file_;
	bool dock_min_dock_;

	core::Size stats_cutoff_;
	bool mutate_framework_for_cluster_;

};
}
}
}



#endif //INCLUDED_protocols_antibody_design_AntibodyDesignMover.hh
