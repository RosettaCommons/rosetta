// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AbrelaxApplication
/// @brief Application-level code for Abrelax, fold_cst and JumpingFoldCst protocols
/// @details
///    use -help to see options
///    usage of class:
///    call AbrelaxApplication::register_options() before core::init::init
///    in main program make instance and call run() method.
///
/// @author Oliver Lange

#ifndef INCLUDED_protocols_abinitio_AbrelaxApplication_hh
#define INCLUDED_protocols_abinitio_AbrelaxApplication_hh

// Unit Headers

// Package Headers

// Project Headers
#include <protocols/jobdist/JobDistributors.hh> // keep first

#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#include <protocols/jobdist/Jobs.fwd.hh>

#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileOptions.fwd.hh>

#include <core/scoring/constraints/ConstraintSet.fwd.hh>
//#include <core/scoring/constraints/ConstraintForest.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>

#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PCA.fwd.hh>

#include <protocols/abinitio/Protocol.fwd.hh>
#include <protocols/abinitio/Templates.fwd.hh>

/// The instance of Loops contained by AbrelaxApplication should be replaced by a LoopsOP
#include <protocols/loops/Loops.hh>

/// The instance of CheckPointer contained by AbrelaxApplication should be replaced by a CheckPointerOP
#include <protocols/checkpoint/CheckPointer.hh>

#include <protocols/jumping/JumpSetup.fwd.hh>
#include <protocols/jumping/MembraneJump.fwd.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>

// ObjexxFCL Headers

// Utility headers
// #include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <string>
#include <vector>

#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

/// @brief application level code  for Abrelax, Foldconstraints and JumpingFoldconstraints
/// WARNING WARNING WARNING. THREAD UNSAFE. INVOKES ConstraintFactory::replace_creator.
/// CODE THAT ABUSES SINGLETONS LIKE THIS OUGHT TO BE SHOT.
class AbrelaxApplication  {
public:
	AbrelaxApplication();

	/// @brief Explicit virtual destructor since AbrelaxApplication contains OPs
	/// NOTE: any time you define a class that is derived from by other classes
	/// and which contains polymorphic functions, it needs to have a virtual destructor.
	/// If your class derives from ReferenceCount, then it will inherit a virtual destructor
	/// If it does not, as AbrelaxApplication does not, then you must declare the
	/// destructor virtual.
	virtual ~AbrelaxApplication();

	/// @brief Explicit copy constructor since AbrelaxApplication contains OPs
	AbrelaxApplication( AbrelaxApplication const & );

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
		core::io::silent::SilentStruct&
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

	/// @brief setup everything needed for fold() --- calls helper functions below
	void setup_fold( core::pose::Pose &extended_pose, ProtocolOP& prot_ptr );

	/// @brief run abrelax-type protocols
	void fold( core::pose::Pose& extended_pose, ProtocolOP prot_ptr );

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


	/// @brief read in membrane topology
	void setup_membrane_topology( core::pose::Pose & pose, std::string spanfile ) const;

	/// @brief initialize template_
	void setup_templates();

	/// @brief insert fragments from aligned regions
	void insert_template_frags( core::pose::Pose&, core::kinematics::MoveMapOP movemap, std::string tag/*for logs*/ ) const;

	void initialize_constraint_forest( core::pose::Pose & pose );

private:
	/// ------------- Data -------------------------------
	/// -------- When you add new data to this class, ----
	/// -------- you must update the copy constructor ----

	// a score file ( written to in process_decoy )
	core::io::silent::SilentFileOptionsOP silent_options_;
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

	//core::scoring::constraints::ConstraintForestOP constraint_forest_;
	//membrane jumping? do we need to fiddle around with the fold tree...

	jumping::MembraneJumpOP membrane_jumps_;

	// jump definitions --- if available
	jumping::BaseJumpSetupOP jump_def_;

	// ss_def ( bascially used for loop_fraction )
	core::fragment::SecondaryStructureOP ss_def_;

	// info about homologues structures --- if available
	TemplatesOP templates_;

	//probably 9mer fragments
	core::fragment::FragSetOP fragset_large_;

	//probably 3mer fragments top25=25
	core::fragment::FragSetOP fragset_small_top25_;

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
