// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/ConsensusDesignMover.hh
/// @brief header file for ConsensusDesignMover
/// @author Florian Richter (floric@u.washington.edu), april 2011

#ifndef INCLUDED_protocols_simple_moves_ConsensusDesignMover_hh
#define INCLUDED_protocols_simple_moves_ConsensusDesignMover_hh

// Unit header
#include <protocols/simple_moves/ConsensusDesignMover.fwd.hh>

// Project headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/sequence/SequenceProfile.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

/// @brief This mover will modify a given task according to a sequence profile
/// and then call the PackRotamersMover.
/// At every position that is designable in the task, AAs that have a probability > min_aa_probability_
/// and higher than the native in the sequence profile  will be allowed
class ConsensusDesignMover : public moves::Mover {

public:

	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef protocols::moves::MoverOP MoverOP;


public:

	// default constructor
	ConsensusDesignMover();

	ConsensusDesignMover(
		core::pack::task::PackerTaskCOP ptask,
		core::scoring::ScoreFunctionCOP sfxn
	);

	~ConsensusDesignMover() override;

	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	// @brief main operation
	void apply( core::pose::Pose & pose ) override;

	void
	set_sasa_cutoff(
		core::Real cutoff ) {
		sasa_cutoff_ = cutoff; }

	void
	set_invert_task( bool setting ){
		invert_task_ = setting; }

	void
	set_use_seqprof_constraints( bool setting ){
		use_seqprof_constraints_ = setting; }

	core::pack::task::PackerTaskCOP
	create_consensus_design_task(
		core::pose::Pose const & pose
	);

	core::scoring::constraints::ConstraintCOPs
	create_sequence_profile_constraints(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask const & task
	) const;



	void parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & ) override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	core::pack::task::PackerTaskCOP ptask_;
	core::pack::task::TaskFactoryOP task_factory_;
	core::scoring::ScoreFunctionCOP sfxn_;

	bool invert_task_; //if a task has been specified externally and this variable is true, all non-packable positions in the task will become design positions
	bool use_seqprof_constraints_; // add a sequence profile constraint during the packing step
	core::Real sasa_cutoff_; //only residues that have sasa above the cutoff will be touched
	core::sequence::SequenceProfileCOP seqprof_;
	bool ignore_pose_profile_length_mismatch_;
};


} // moves
} // protocols


#endif
