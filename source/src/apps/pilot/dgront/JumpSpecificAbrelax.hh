// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file AbrelaxApplication
/// @brief Application-level code for Abrelax, fold_cst and JumpingFoldCst protocols
/// @detailed
///	   use -help to see options
///    usage of class:
///    call AbrelaxApplication::register_options() before core::init::init
///    in main program make instance and call run() method.
///
/// @author Oliver Lange

#ifndef INCLUDED_apps_pilot_dgront_JumpSpecificAbrelax_hh
#define INCLUDED_apps_pilot_dgront_JumpSpecificAbrelax_hh

// Unit Headers

// Package Headers

// Project Headers
#include <protocols/jobdist/JobDistributors.hh> // keep first

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <protocols/jobdist/Jobs.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileData.fwd.hh>

#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintForest.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/fragment/FragSet.hh>

#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/PCA.hh>

#include <protocols/abinitio/ClassicAbinitio.fwd.hh>
#include <protocols/Protocol.fwd.hh>
#include <protocols/abinitio/Templates.hh>
#include <protocols/loops/LoopClass.hh>

#include <protocols/checkpoint/CheckPointer.hh>

#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/SecondaryStructure.fwd.hh>


// ObjexxFCL Headers

// Utility headers
// #include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <vector>

namespace protocols {
namespace abinitio {

/// @brief application level code  for Abrelax, Foldconstraints and JumpingFoldconstraints

class JumpSpecificAbrelax  {
public:
	JumpSpecificAbrelax();
	static void register_options();

	/// @brief diagnostic stuff, i.e., computing numbers like an RMSD for each decoy and storing in silent-score-file
	/// put everything in here. --- actually it mainly calls the evaluator_->apply method.
	/// add diagnostic stuff either here as explicit code or
	/// in form of a PoseEvaluator to evaluator_ ( see setup () )
	/// the latter has the advantage that the specific evaluation can be carried out during the run of the protocol
	/// e.g., for abinitio:debug ( output to stage1_outfile stage2_outfile... )
	void process_decoy(
		core::pose::Pose &pose,
		core::scoring::ScoreFunction const&,
		std::string tag,
		core::io::silent::ProteinSilentStruct&
	) const;

	/// @brief read constraint set (self-initializing) and connect it to pose
	void add_constraints( core::pose::Pose &pose );

	/// @brief initialization of application: read some pdb files,  set evaluator_
	void setup();

	/// @brief run application code: calls do_rerun() or fold() depending on cmd-options
	void run();

	/// @brief run process_decoy on all poses in silent-in file
	void do_rerun();

	/// @brief run process_decoy on all poses in silent-in file -- use of JobDistributor
	void do_distributed_rerun();

	/// @brief run abrelax-type protocols
	void fold();

	/// @brief return pose with simple fold-tree that has small <0.1 RMSD to input pose
	bool close_loops( core::pose::Pose &pose, core::scoring::ScoreFunctionOP scorefxn, std::string const& tag );

	/// @brief relax structure ( fast / classic as controlled by cmd-options )
	void relax( core::pose::Pose &pose, core::scoring::ScoreFunctionOP, std::string const& tag );

	/// @brief relax multiple structures that are stored in abinitio_protocol.structure_store
	bool multi_fast_relax(
	 Protocol& abinitio_protocol,
	 core::scoring::ScoreFunctionOP,
	 jobdist::PlainSilentFileJobDistributor jobdist,
	 int& curr_nstruct,
	 jobdist::BasicJobOP& curr_job
	);

	/// @brief little helper: minimize structure to have lower chainbreak score
	/// ( seems particularly necessary after reading from silent-file )
	//obsolet  void fix_chainbreaks( core::pose::Pose &pose );

	/// @brief add a PoseEvaluator derived instance for decoy-processing
	void add_evaluation( evaluation::PoseEvaluatorOP );

	/// @brief check if the given pose passes the set of abinitio filters.
	bool check_filters( core::pose::Pose & pose );

private:
	/// @brief create score-functions for centroid and fullatom level
	core::scoring::ScoreFunctionOP generate_scorefxn( bool fullatom = false );

	/// @brief setup everything needed for fold() --- calls helper functions below
	void setup_fold( core::pose::Pose &extended_pose, ProtocolOP& prot_ptr );

	/// ---- Helper functions for setup_fold

	/// @brief steal native torsions from native_pose_ and apply to the extended_pose.
	void copy_native_structure( core::pose::Pose &extended_pose ) const;

	/// @brief copy torsions from the desired_pose, copy them into the extended_pose.
	void copy_structure( core::pose::Pose & extended_pose, core::pose::Pose & desired_pose ) const;

	/// @brief steal native torsions from native_pose_ and apply to the "extended_pose"
	void generate_extended_pose( core::pose::Pose &extended_pose, std::string const& sequence ) const;

	/// @brief read fragment data
	void setup_fragments(); // core::fragment::FragSetOP& fragset_large, core::fragment::FragSetOP& fragset_small );

	/// @brief read jump definitions and set jump_def_
	void setup_jumps( core::pose::Pose const& extended_pose );

	/// @brief initialize template_
	void setup_templates();

	/// @brief insert fragments from aligned regions
	void insert_template_frags( core::pose::Pose&, core::kinematics::MoveMapOP movemap, std::string tag/*for logs*/ ) const;

	void initialize_constraint_forest( core::pose::Pose & pose );

	/// ------------- Data ------------------------------

	// a score file ( written to in process_decoy )
	core::io::silent::SilentFileDataOP silent_score_file_;

	// native_pose: steal fragments, compute rmsd, start-structure, sequence
	core::pose::PoseOP native_pose_;

	// start-pose non-extended pose for start of runs
	core::pose::PoseOP init_pose_obj_;

	// loops ( yes, we can do loop-modelling now )
	loops::Loops loops_in_;

	// if specified the structures are projected to PCA-eigenvector
	evaluation::PCA_OP pca_;

	// are we doing fa-relax ?
	bool bRelax_;

	// the sequence of the target protein
	std::string sequence_;

	// the constraint set --- if available
	core::scoring::constraints::ConstraintSetOP cstset_;

	core::scoring::constraints::ConstraintForestOP constraint_forest_;

	// jump definitions --- if available
	jumping::BaseJumpSetupOP jump_def_;

	// ss_def ( bascially used for loop_fraction )
	jumping::SecondaryStructureOP ss_def_;

	// info about homologues structures --- if available
	TemplatesOP templates_;

	//probably 9mer fragments
	core::fragment::FragSetOP fragset_large_;

	//probably 3mer fragments
	core::fragment::FragSetOP fragset_small_;

	// pure template fragments... or merged with fragset_large_ only used if vary_frag_size == true
	core::fragment::FragSetOP fragset_templates_;

// a bunch of PoseEvaluators for process_decoy() --- if available
	evaluation::MetaPoseEvaluatorOP evaluator_;

	// checkpoints for close_loop
	checkpoint::CheckPointer abrelax_checkpoints_;
};

} //abinitio
} //protocols

#endif
