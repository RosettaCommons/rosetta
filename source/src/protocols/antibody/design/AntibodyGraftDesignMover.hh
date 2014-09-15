// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyGraftDesignMover.hh
/// @brief Class that initially designs antibodies through grafting using an AntibodyDatabase + North_AHO numbering scheme
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_design_AntibodyGraftDesignMover_hh
#define INCLUDED_protocols_antibody_design_AntibodyGraftDesignMover_hh

// Project Includes
#include <protocols/antibody/design/AntibodyGraftDesignMover.fwd.hh>
#include <protocols/antibody/design/AntibodyDesignEnum.hh>
#include <protocols/antibody/design/AntibodyDesignModeler.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>

// Protocol Includes
#include <protocols/grafting/CCDEndsGraftMover.fwd.hh>

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
	using namespace protocols::antibody;
	using namespace protocols::grafting;
	using namespace protocols::antibody::clusters;
	
	using namespace utility;
	using namespace core;
	using namespace core::scoring;
	using std::map;

	///@brief This struct gives antibody designers complete control either through this classes interface or an instruction file.
	struct CDRGraftInstructions {
		bool graft; //Use grafting algorithm on this CDR?
		bool stay_native_cluster; //Only sample within the native cluster?
		bool cluster_centers_only;//Only use the center clusters for sampling?
		bool min_neighbor_sc; //Minimize neighbor sidechains during chosen minimization if possible
		bool min_rb; // Minimize rigid body antigen-antibody interface during chosen minimization if possible
		
		MinTypeEnum mintype; //What to do after graft?
		std::map< Size,  bool > cluster_types; //Only sample within cluster types (Cluster types 1,2,3)
		vector1< std::string > leave_out_pdb_ids;
		vector1< std::string > include_only_pdb_ids;
		vector1< CDRClusterEnum > leave_out_clusters;//Leave out these clusters.  
		vector1< CDRClusterEnum > include_only_clusters; //Include only these clusters.
		core::Size min_length;
		core::Size max_length;
	};

	//Sampling weights not currently implemented
	struct SamplingWeights {
		std::map< Size,  Real > types; //Used to sample more of one type then others.  Sample most in type1, but use a bit of type2/3.
		Real native_cluster; //May be used to weight the sampling of native clusters vs others.
		Real cdr; //Used to sample more of one or particular CDRs rather then others (H3)
		Real center; //Used to sample more of cluster centers then other CDRs.
	};	
	
	
	typedef map< CDRNameEnum, vector1< pose::PoseOP > > CDRSet;
	typedef map< CDRNameEnum, vector1< CDRClusterEnum > > CDRClusterMap; 
	typedef std::map< CDRNameEnum, vector1< std::string > > PDBMap;
	
	
	
	
	
	
	
///@brief This class designs antibodies by grafting, using cluster information and constraints to help.
/// It represents the first step in The Rosetta Antibody Designer, but can be used outside of the application.
/// The Antibody Database that this class uses will be available online (Too large for the rosetta_database)
/// Till then, email me for access to the github repository.
///
///@details To use this class:
///			1) Use default instruction path, or set one.  Loads CDRs from AntibodyDatabase
///			    *See rosetta_database/sampling/antibodies for file format.
///			2) Use class interface settings to control sampling.  Loads CDRs from AntibodyDatabase.
///			3) Use set_cdr_sets to set your own CDRs to graft.
///  
class AntibodyGraftDesignMover: public protocols::moves::Mover {

public:
	
	AntibodyGraftDesignMover(AntibodyInfoOP ab_info);
	
	AntibodyGraftDesignMover(AntibodyInfoOP ab_info, std::string instruction_path);
			
	virtual ~AntibodyGraftDesignMover();
	
	///@brief Reads default CDRGraftInstruction file and creates structs.  
	void
	set_defaults();
	
	//Parse my tag
	
	virtual void
	apply(pose::Pose & pose);
	
	
public:
	////////////////////////////////////////////////////////////////////////////
	// Full Settings
	//
	//
	
	///@brief Will not initialize CDRs from the AntibodyDatabase.  Use if you have your own CDR's you are interested in grafting.
	///@details Overhang residues will be used for superposition.  To identify your CDRs, use the functions in AntibodyInfo or the pilot app identify_cdr_clusters
	/// in pilot/jadolfbr/cluster_utils
	void
	set_cdr_set(CDRSet & cdr_set, CDRClusterMap & cdr_cluster_map, core::Size overhang); 
	
public:
	////////////////////////////////////////////////////////////////////////////
	// Modeling Settings
	//
	//
	
	void
	set_scorefunction(ScoreFunctionOP scorefxn);
	
	void
	set_graft_rounds(core::Size graft_rounds);
	
	//void
	//set_filters();
	
	///@brief Options are: relax, minimize, repack. Default is minimize
	void
	set_mintype(CDRNameEnum const cdr_name, MinTypeEnum mintype);
	
	void
	set_mintype_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, MinTypeEnum mintype);
	
	
	///@brief  Use rigid body optimization during minimization of the CDR (Minimize jump bt antigen and antibody)?
	///@details dock_chains is LH_antigen.  Does nothing if repacking the cdr.  Default False.
	void
	set_min_rb(CDRNameEnum const cdr_name, bool const setting);
	
	void
	set_min_rb_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool const setting);
	
	
	///@brief Use Neighbor sidechains during minimization of this CDR post-graft?  
	void
	set_min_neighbor_sc(CDRNameEnum const cdr_name, bool const setting);
	
	void
	set_min_neighbor_sc_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool const setting);
	
	
	///@brief Set the algorithm to run a low-resolution docking step after each graft.  Default false.
	///@details Uses command-line options for docking.  Lower inner or outer cycles if you are grafting many CDRs.  Still quicker then relax.
	void
	set_dock_post_graft(bool dock_post_graft);
	
	void
	set_dock_rounds(core::Size dock_rounds);
	
	///@brief Set the algorithm to run a packing step with neighbor detection post graft.  Before any minimization or docking.  
	void
	set_pack_post_graft(bool pack_post_graft);
	
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
	// Cluster Settings
	//
	//
	
	///@brief Leave out or include only clusters added here
	void
	set_cluster(CDRNameEnum const cdr, CDRClusterEnum const cluster, bool const setting);
	
	///@brief Leave out or include only clusters added here
	void
	set_cluster_range(
		CDRNameEnum const cdr, 
		CDRClusterEnum const cluster_start,
		CDRClusterEnum const cluster_end,
		bool const setting);

	
public:
	////////////////////////////////////////////////////////////////////////////
	// CDR Settings
	//
	//
	
	///@brief Use this CDR during design.
	void
	set_cdr(CDRNameEnum const cdr, bool const setting);
	
	///@brief Use this range of CDRs during design.  Default is all of them.
	void
	set_cdr_range(CDRNameEnum const cdr_start, CDRNameEnum cdr_end, bool const setting);
	

	///@brief Set to only sample within the native cluster of this CDR.  
	void
	set_cdr_stay_in_native_cluster(CDRNameEnum const cdr, bool const setting);
	
	///@brief Set a range of CDRs for which to only sample their native clusters.
	void
	set_cdr_range_stay_in_native_cluster(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool const setting);
	
	
	///@brief Set to only sample with clusters of the given type for this CDR. Default is to only use type 1.
	void
	set_cdr_stay_in_type(CDRNameEnum const cdr, core::Size const type, bool const setting);
	
	///@brief set a range of CDRs for which to stay within a particular type. Default is to only use type 1.
	void
	set_cdr_range_stay_in_type(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, core::Size const type, bool const setting);
	
	
	///@brief set the minimum cdr length to sample.  Nothing shorter then this will be used during graft.
	void
	set_cdr_min_length(CDRNameEnum const cdr, core::Size length);
	
	void
	set_cdr_min_length_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, core::Size length);
	
	///@brief set the maximum cdr length to sample.  Nothing longer then this will be used.  Useful for H3.
	void
	set_cdr_max_length(CDRNameEnum const cdr, core::Size length);
	
	void
	set_cdr_max_length_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, core::Size length);
	
	
	///@brief Only sample cluster centers for this CDR
	void
	set_cdr_cluster_centers_only(CDRNameEnum const cdr, bool setting);
	
	void
	set_cdr_cluster_centers_only_range(CDRNameEnum const cdr_start, CDRNameEnum const cdr_end, bool setting);
	
	///@brief Leave out these PDBs from sampling.  Mainly for testing, but who knows.
	//void
	//leave_out_pdb_ids(CDRNameEnum const cdr, vector1< std::string > pdbids);
	
public:
	////////////////////////////////////////////////////////////////////////////
	// Sampling Settings
	//
	//
	
	//void
	//set_type_weight(CDRNameEnum const cdr, Size type, Real weight);
	
	//void
	//set_native_cluster_weight(CDRNameEnum const cdr, Real weight);
	
	//void
	//set_cdr_weight(CDRNameEnum const cdr, Real weight);
	
	//void
	//set_center_weight(CDRNameEnum const cdr, Real weight);

public:
	////////////////////////////////////////////////////////////////////////////////
	// Boiler Plate
	//
	//
	
	virtual string get_name() const;

	
public:
	////////////////////////////////////////////////////////////////////////////////
	// Rosetta Scripts
	//
	//
	
private:
	
	void
	read_command_line_options();
	
	///@brief Call DesignInstructionParser to determine settings for the design run.
	void
	read_instructions(std::string instruction_path);
	
	//Determine each CDRs cluster if not already calculated and set in AntibodyInfo instance.
	void
	setup_native_clusters(pose::Pose & pose);
	
	///@brief Uses instructions to Query the AntibodyDatabase and load poses.  
	void
	initialize_cdr_set();
	
	void
	check_for_top_designs(pose::Pose & pose);
	
	///@brief Grafts a single CDR. Index is the vector index of CDRSet/CDRClusterMap
	///@details Return success or failure
	bool
	graft_cdr(pose::Pose & pose, CDRNameEnum cdr, core::Size index);
	
	
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
	
	//void
	//run_non_repeating_mc_algorithm(pose::Pose & pose, vector1< CDRNameEnum> & cdrs_to_design);
	
	///@brief Fix PDB info as somehow the code for pdbinfo changed and my other method no longer works.  No idea why.
	void
	fix_pdb_info(pose::Pose & pose, CDRNameEnum cdr, CDRClusterEnum cluster, core::Size original_start, core::Size original_pdb_end);
	
	///@brief Extend the native CDRs to be designed for benchmarking.
	void
	extend_native_cdrs(pose::Pose & pose, vector1<CDRNameEnum> & cdrs_to_design);
	
	void
	set_default_graft_settings();
	
	///@brief Gets a list of vectors whose indexes correspond to CDRNameEnum, and whose values correspond to the cdr_set index.  If the value is 0, it means no cdr in set.
	vector1< vector1< Size > >
	get_cdr_set_index_list();
	
	//void
	//run_stochastic_graft_algorithm(pose::Pose & pose, vector1< CDRNameEnum > & cdrs_to_design);
	
	
	map< CDRNameEnum, CDRClusterOP> native_clusters_; //Native North cluster
	map< CDRNameEnum,  CDRGraftInstructions > cdr_instructions_;
	map< CDRNameEnum, SamplingWeights > sampling_instructions_;
	
	map< CDRNameEnum, utility::vector1< core::scoring::constraints::ConstraintCOP > >constraint_map_; //Currently set constraints.  Used for removal when nessessary.
	
	std::string instruction_path_;
	CDRSet cdr_set_; //Actual poses that will be grafted.
	CDRClusterMap cdr_cluster_map_; //Maps CDRSet to CDRClusterEnums to constrained relaxes
	PDBMap pdbmap_; //Maps CDRSet to PDB tags
	
	AntibodyInfoOP ab_info_;
	CCDEndsGraftMoverOP graft_mover_; 
	ScoreFunctionOP scorefxn_;
	AntibodyDesignModelerOP modeler_;
	protocols::moves::MonteCarloOP mc_;
	
	core::Size overhang_;
	core::Size graft_rounds_;
	core::Size dock_rounds_;
	core::Size num_top_designs_; //Number of top designs to keep.
	
	//Can be a struct.  Just not now.
	vector1< pose::PoseOP > top_designs_; 
	vector1< core::Real> top_scores_;
	
	core::Size total_permutations_; //Total number of possible combinations
	
	//Overall Booleans
	bool dock_post_graft_; //Run a low-resolution docking step after the graft?
	bool pack_post_graft_; //Run a packing step with neighbor detection after the graft.  
	bool rb_min_post_graft_; //Run an rb_min step post graft.  Useful if not using docking.
	bool initial_perturb_; //Run DockingInitialPerturber post graft
	bool use_deterministic_algorithm_;
	core::Real max_linear_chainbreak_;  //Sometimes the graft completely fails unfortunately.  May edit the graft code itself soon.

	//Benchmarking
	bool extend_native_cdrs_;
	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// bool benchmark_;
};
}
}
}



#endif	//INCLUDED_protocols_antibody_design_AntibodyGraftDesignMover.hh
