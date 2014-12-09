// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyGraftDesignMover.hh
/// @brief Class that initially designs antibodies through grafting using an AntibodyDatabase + North_AHO numbering scheme
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_design_AntibodyGraftDesignMover_hh
#define INCLUDED_protocols_antibody_design_AntibodyGraftDesignMover_hh

// Project Includes
#include <protocols/antibody/design/AntibodyGraftDesignMover.fwd.hh>
#include <protocols/antibody/design/AntibodyDesignEnum.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/design/AntibodyDesignMoverGenerator.fwd.hh>
#include <protocols/antibody/database/CDRSetOptions.hh>
#include <protocols/antibody/design/CDRGraftDesignOptions.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>

#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>

#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/constraints/ParatopeSiteConstraintMover.fwd.hh>
#include <protocols/antibody/constraints/ParatopeEpitopeSiteConstraintMover.fwd.hh>

// Protocol Includes
#include <protocols/grafting/CCDEndsGraftMover.fwd.hh>
#include <protocols/grafting/AnchoredGraftMover.fwd.hh>

// Core Includes
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/Constraint.hh>

// Protocol Includes
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/types.hh>

namespace protocols{
namespace antibody{
namespace design{

	//Sampling weights not currently implemented
	//struct SamplingWeights {
	//	std::map< Size,  Real > types; //Used to sample more of one type then others.  Sample most in type1, but use a bit of type2/3.
	//	Real native_cluster; //May be used to weight the sampling of native clusters vs others.
	//	Real cdr; //Used to sample more of one or particular CDRs rather then others (H3)
	//	Real center; //Used to sample more of cluster centers then other CDRs.
	//};	
	

///@brief This class designs antibodies by grafting, using cluster information and constraints to help.
/// It represents the first step in The Rosetta Antibody Designer, but can be used outside of the application.
/// A 2011 Antibody Database is included in the rosetta datase.  Up-to-date versions can be downloaded from 
/// 
///
///@details To use this class:
///	1) Use default instruction path, or set one.  Loads CDRs from AntibodyDatabase
///	    *See rosetta_database/sampling/antibodies for file format.
///	2) Use class interface settings to control sampling.  Loads CDRs from AntibodyDatabase.
///	3) Use set_cdr_sets to set your own CDRs to graft.
///  
class AntibodyGraftDesignMover: public protocols::moves::Mover {

public:
	
	AntibodyGraftDesignMover(AntibodyInfoOP ab_info);
			
	virtual ~AntibodyGraftDesignMover();
	
	///@brief Reads default CDRGraftInstruction file and creates structs.  
	void
	set_defaults();
	
	//Parse my tag
	
	virtual void
	apply(pose::Pose & pose);
	
	
public:
	
	///@brief Set the options which will be used for querying the database
	void
	set_cdr_set_options(AntibodyCDRSetOptions cdr_set_options);
	
	///@brief Set CDR-specific GraftDesign options
	void
	set_graft_design_options(AntibodyCDRGraftDesignOptions graft_design_options);
	
	///@brief Set CDR-specific SeqDesign options
	void
	set_seq_design_options(AntibodyCDRSeqDesignOptions seq_design_options);
	
	///@brief Will not initialize CDRs from the AntibodyDatabase.  Use if you have your own CDR's you are interested in grafting.
	///@details Overhang residues will be used for superposition.  To identify your CDRs, use the functions in AntibodyInfo or the pilot app identify_cdr_clusters
	/// in pilot/jadolfbr/cluster_utils
	void
	set_cdr_set(CDRSet & cdr_set, core::Size overhang); 
	
public:
	////////////////////////////////////////////////////////////////////////////
	// Modeling Settings
	//
	//
	
	void
	set_scorefunction(core::scoring::ScoreFunctionOP scorefxn);
	
	void
	set_graft_rounds(core::Size graft_rounds);
	
	void
	set_post_graft_modeling_cycles(core::Size cycles);
	
	
	///@brief Set the algorithm to run a low-resolution docking step after each graft.  Default false.
	///@details Uses command-line options for docking.  Lower inner or outer cycles if you are grafting many CDRs.  Still quicker then relax.
	void
	set_dock_post_graft(bool dock_post_graft);
	
	void
	set_dock_rounds(core::Size dock_rounds);
	
	///@brief Set the algorithm to run a final rigid-body minimization of antigen/antibody interface post graft.  Useful if not docking - quicker, but less reliable then full dock.
	void
	set_rb_min_post_graft(bool rb_min_post_graft);
	
	
	///@brief Sets the protocol to keep a specific number of top designs.  Default is 10
	void
	set_keep_top_designs(core::Size top_designs);
	
	///@brief Get the top designs found.  You can then use them in other protocols, dump them, etc. They are in order. 
	///@details - This should be refactored to get_additional_output. 
	vector1< pose::PoseOP>
	get_top_designs(){
		return top_designs_;
	};
	

public:
	////////////////////////////////////////////////////////////////////////////
	// Paratope and Epitope constraints
	//
	//
	
	
	///@brief Set any paratope CDRs.  If not set, will use all CDRs as the paratope where needed.
	/// Used mainly for constraints
	void
	set_paratope_cdrs(vector1<bool> const & cdrs);
	
	///@brief Set any epitope residues in PDB numbering.  If not set, they will be detected automatically.
	void
	set_epitope_residues(vector1<PDBNumbering > epitope_residues);
	
public:
	
	virtual string get_name() const;

	virtual void
	show(std::ostream & output=std::cout) const;
	
private:
	
	void
	read_command_line_options();
	
	///@brief Setup ALL options classes
	void
	setup_options_classes();
	
	//Determine each CDRs cluster if not already calculated and set in AntibodyInfo instance.
	void
	setup_native_clusters(pose::Pose & pose);
	
	void
	setup_paratope_epitope_constraints(core::pose::Pose & pose);
	
	void
	setup_scorefxn();
	
	///@brief Extend the native CDRs to be designed for benchmarking.
	void
	setup_random_start_pose(pose::Pose & pose, vector1<CDRNameEnum> & cdrs_to_design);
	
	void
	setup_default_graft_settings();
	
	
	///@brief Uses instructions to Query the AntibodyDatabase and load poses.  
	void
	initialize_cdr_set(pose::Pose const & pose);
	
	void
	check_for_top_designs(pose::Pose & pose);
	
	///@brief Gets a list of vectors whose indexes correspond to CDRNameEnum, and whose values correspond to the cdr_set index.  If the value is 0, it means no cdr in set.
	vector1< vector1< Size > >
	get_cdr_set_index_list();
	

	///@details Applies graft, modeling, etc. to a single CDR.
	bool
	apply_to_cdr(pose::Pose & pose, CDRNameEnum cdr, core::Size index, bool min_post_graft = true);
	
	///@brief Grafts a single CDR into framework. Index is the vector index of CDRSet/CDRClusterMap. 
	///@details Return success or failure
	bool
	graft_in_cdr(pose::Pose & pose, CDRNameEnum const cdr, CDRPose & cdr_pose);
	
	//Runs graft then cartmin on AnchoredGraftMover or CCDEndsGraftMover.  Returns boolean of closure.
	std::pair<bool, core::Size>
	run_graft(pose::Pose & pose, CDRNameEnum const cdr, CDRPose & cdr_pose, grafting::AnchoredGraftMoverOP grafter);
	
	void
	run_post_graft_min(pose::Pose & pose, protocols::moves::MonteCarlo & mc, CDRNameEnum cdr);
	
	///@brief Mutates framework residues needed to stay within a particular cluster.  Only one  (L1-11-1) is known to have framework dependencies.  For now.
	/// Will be replaced by AntibodyDesignOptimizer
	//void
	//mutate_framework_residues(pose::Pose pose, CDRClusterEnum cluster);
	
private:
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Main Algorithms
	//
	//
	
	///@brief.  If rounds >= number of possible combinations - Try them all. 
	///@details Grafts CDRs on the input structure, in a random order.  
	void
	run_deterministic_graft_algorithm(pose::Pose & pose, vector1<CDRNameEnum>& cdrs_to_design);
	
	///@brief Basic mc algorithm that randomly samples from the cdr set.
	void
	run_basic_mc_algorithm(pose::Pose & pose, vector1<CDRNameEnum>& cdrs_to_design);
	
	
	
	void
	print_str_vec(std::string const name, utility::vector1<std::string> const & vec, std::ostream & output=std::cout) const;
	
	
private:
	
	AntibodyInfoOP ab_info_;
	AntibodyCDRSetOptions cdr_set_options_;
	AntibodyCDRGraftDesignOptions cdr_graft_design_options_;
	AntibodyCDRSeqDesignOptions cdr_seq_design_options_;
	AntibodySeqDesignTFCreatorOP seq_design_creator_;
	
	std::map< CDRNameEnum, vector1< CDRPose > > cdr_set_;
	
	protocols::grafting::CCDEndsGraftMoverOP graft_mover_;
	protocols::grafting::AnchoredGraftMoverOP anchored_graft_mover_;
	
	core::scoring::ScoreFunctionOP scorefxn_;
	AntibodyDesignModelerOP modeler_;
	AntibodyDesignMoverGeneratorOP mover_generator_;
	
	protocols::moves::MonteCarloOP mc_;
	
	constraints::ParatopeEpitopeSiteConstraintMoverOP paratope_epitope_cst_mover_;
	constraints::ParatopeSiteConstraintMoverOP paratope_cst_mover_;
	
	core::Size overhang_;
	core::Size graft_rounds_;
	core::Size dock_rounds_;
	core::Size num_top_designs_; //Number of top designs to keep.
	
	//Can be a struct.  Just not now.
	utility::vector1< pose::PoseOP > top_designs_; 
	utility::vector1< core::Real> top_scores_;
	
	core::Size total_permutations_; //Total number of possible combinations
	
	//Overall Booleans
	bool dock_post_graft_; //Run a low-resolution docking step after the graft?
	bool rb_min_post_graft_; //Run an rb_min step post graft.  Useful if not using docking.
	bool initial_perturb_; //Run DockingInitialPerturber post graft
	bool use_deterministic_algorithm_;
	bool design_post_graft_;
	//bool design_neighbor_cdrs_;
	
	//Any cmd-line set Paratope and Epitope residues
	vector1<bool> paratope_cdrs_;
	utility::vector1<PDBNumbering > epitope_residues_; //Vector of resnum, chain pairs.
	
	bool adapt_graft_;
	bool benchmark_;
	bool use_light_chain_type_;
	bool use_epitope_constraints_;
	bool print_tracer_info_; //Quick change to not print tracer info on  random starting cdr
	
	core::Size post_graft_modeling_cycles_;
	
	
};
}
}
}



#endif	//INCLUDED_protocols_antibody_design_AntibodyGraftDesignMover.hh
