// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/MPIWorkPoolJobDistributor.hh
/// @brief  header for MPIWorkPoolJobDistributor - intended for continuous resamplig jobs  that spawn new jobs based on a pool/archive of
///         structures
/// @author Oliver Lange olange@u.washington.edu

#ifndef INCLUDED_protocols_jd2_archive_EvaluatedArchive_hh
#define INCLUDED_protocols_jd2_archive_EvaluatedArchive_hh

// Unit headers
#include <protocols/jd2/archive/ArchiveManager.fwd.hh>
#include <protocols/jd2/archive/ArchiveBase.hh>
#include <protocols/jd2/archive/VarianceStatisticsArchive.fwd.hh>
//#include <protocols/jd2/archive/EvaluatedArchive.fwd.hh>

// Package headers
#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.fwd.hh>


// Utility headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>

// C++ headers
#include <string>
#include <list>

#include <core/scoring/ResidualDipolarCoupling.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {
namespace archive {
//class ArchiveManager;


/// @brief Tags used to tag messeges sent by MPI functions used to decide whether a slave is requesting a new job id or
///flagging as job as being a bad input

/// @details This job distributor is meant for running jobs where the machine you are using has a large number of
///processors, the number of jobs is much greater than the number of processors, or the runtimes of the individual jobs
///could vary greatly. It dedicates the head node (whichever processor gets processor rank #0) to handling job requests
///from the slave nodes (all nonzero ranks). Unlike the MPIWorkPartitionJobDistributor, this JD will not work at all
///without MPI and the implementations of all but the interface functions have been put inside of ifdef directives.
///Generally each function has a master and slave version, and the interface functions call one or the other depending
///on processor rank.

class EvaluatedArchive : public ArchiveBase {
	typedef ArchiveBase Parent;
public:

	/// @brief Constructor  and Destructor
	EvaluatedArchive( ArchiveManagerAP ptr );
	EvaluatedArchive();
	~EvaluatedArchive();

	/// @brief Archive specific options
	static void register_options();

	/// @brief add decoy to Archive
	/// @detail evaluate decoy and call add_evaluated_structure
	virtual bool add_structure(
		core::io::silent::SilentStructOP new_decoy,
		core::io::silent::SilentStructOP alternative_decoy,
		Batch const& batch
	);


	/// @brief  compute score according to select_weights --- this can contain any evaluator columns
	core::Real select_score( core::io::silent::SilentStructOP evaluated_decoy );

	/// @brief set common evaluators: eg. ConstraintEvaluator if -cst_file is present
	void setup_default_evaluators();

	/// @brief overloaded that we can sort the pool after reading
	virtual bool restore_from_file();

	/// @brief only overloaded this to add some verbosity each time we read structures
	virtual void read_structures(
		core::io::silent::SilentFileData& sfd,
		core::io::silent::SilentFileData& alternative_decoys,
		Batch const& batch
	);

	/// @brief overloaded to make input decoys appear the same as decoys coming from batches
	virtual void init_from_decoy_set( core::io::silent::SilentFileData const& sfd ) = 0;

	/// @brief typedefs for Evaluators and Weights
	typedef std::map< std::string, core::Real > WeightMap;
	typedef std::map< std::string, evaluation::PoseEvaluatorCOP > EvaluatorMap;

	void start_evaluation_timer() const;
	/// @brief yields an "evaluated" silent-struct which can be queried with select_score
	/// @detail will run scoring-process if evaluate_local() otherwise just returns the intpu-silent-struct
	core::io::silent::SilentStructOP evaluate_silent_struct( core::io::silent::SilentStructOP from_batch ) const;

	/// @brief add an evaluated decoy to Archive (i.e, evaluated_decoy = evaluate( some_decoy ) );
	virtual bool add_evaluated_structure(
		core::io::silent::SilentStructOP evaluated_decoy,
		core::io::silent::SilentStructOP alternative_decoy,
		Batch const&
	);

	/// @brief specify if decoys are evaluated on the master or (non-local i.e., on the individual slave nodes)
	bool evaluate_local() const {
		return b_evaluate_incoming_decoys_;
	}

	void set_evaluate_local( bool setting ) {
		b_evaluate_incoming_decoys_ = setting;
	}

	/// @brief recompute all score-values of all decoys and re-order the archive by (new) select_score
	virtual void rescore();

	/// @brief add new PoseEvaluation to set of evaluators, specify weight for contribution to select_score()
	void add_evaluation( evaluation::PoseEvaluatorCOP, core::Real weight = 0.0 );

	/// @brief remove Evaluator
	void remove_evaluation( std::string const& column );

	/// @brief is a certain elvaluator present ?
	bool has_evaluator( std::string const& column );

	/// @brief set weight of an evaluator or a column otherwise present in silent-structs
	/// (i.e, score, chainbreak, external evaluation like score_final )
	void set_weight( std::string const& column, core::Real weight );

	core::Real get_weight( std::string const& column ) const;

	/// @brief set scorefxn used for evaluation
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn_ );

	core::scoring::ScoreFunction const & scorefxn() const;

	virtual WeightMap const& score_variations() const;

	virtual core::Real score_variation( std::string const& col ) const;

	WeightMap const& weights() const {
		return select_weights_;
	}

	EvaluatorMap const& evaluators() const {
		return evaluators_;
	}

	void set_weights( WeightMap const& setting );
	void set_evaluators( EvaluatorMap const&, WeightMap const& );

	///overloaded to save / restore the variance_archive_
	virtual void save_to_file( std::string suffix = "" );

protected:
	core::scoring::ScoreFunctionOP scorefxn_non_const();

	/// @brief score a pose
	virtual void score( core::pose::Pose& pose ) const;

	virtual void invalidate_score_variations() {}

private:
	/// @brief call score( pose ) and collect energies into result
	/// this is low-level function: it expects that result already contains the coordinates of the pose
	/// for convenience the pointer result is also return value
	virtual core::io::silent::SilentStructOP evaluate_pose( core::io::silent::SilentStructOP result, core::pose::Pose& input_pose ) const;

	/// @brief re-sort decoys based on select_score
	void sort();

	/// @brief scorefxn_ for evaluate( SilentStruct, Pose const& )
	core::scoring::ScoreFunctionOP scorefxn_;

	/// @brief Evaluators and weights for select_score and evaluate
	WeightMap select_weights_;
	EvaluatorMap evaluators_;

	WeightMap dummy_score_variations_;

	/// @brief keep track wether cached scores in _archive_select_score_ are up-to-date
	mutable bool scores_are_clean_; //false after add_evaluation or change of scorefxn_

	// KAB - below line commented out by warnings removal script (-Wunused-private-field) on 2014-09-11
	// mutable bool score_variations_are_clean_;

	/// @brief local evaluation or is evaluation outsourced to slave nodes?
	bool b_evaluate_incoming_decoys_;

	/// @brief keep track whether our options have been registered at start up
	static bool options_registered_;

	mutable time_t start_eval_time_;

	VarianceStatisticsArchiveOP variance_archive_;
};


}//archive
}//jd2
}//protocols


#endif //INCLUDED_protocols_jd2_Archive_HH
