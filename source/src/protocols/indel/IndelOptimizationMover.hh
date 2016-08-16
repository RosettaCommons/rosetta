// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_indel_IndelOptimizationMover_hh
#define INCLUDED_protocols_indel_IndelOptimizationMover_hh

#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/Mover.fwd.hh>

#include <protocols/indel/IndelOptimizationMover.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace indel {

class IndelOptimizationMover : public moves::Mover {

public:

	//default ctor
	IndelOptimizationMover( Size start_res, Size end_res, Size loop_length,
		std::string remodel, std::string intermedrelax, std::string refine, std::string relax, bool frag_files, Size num_to_dock, bool dump_initial_results ):
		Mover("IndelOptimizationMover"),
		start_res_( start_res ),
		end_res_( end_res ),
		loop_length_( loop_length ),
		remodel_( remodel ),
		intermedrelax_( intermedrelax ),
		refine_( refine ),
		relax_( relax ),
		frag_files_( frag_files ),
		num_to_dock_( num_to_dock ),
		dump_initial_results_( dump_initial_results )
	{}

	//default dtor
	virtual ~IndelOptimizationMover(){}

	//methods
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const { return "IndelOptimizationMover"; }

	// getters
	Size start_res() {
		return start_res_;
	}
	Size end_res() {
		return end_res_;
	}
	Size loop_length() {
		return loop_length_;
	}
	std::string remodel() {
		return remodel_;
	}
	std::string intermedrelax() {
		return intermedrelax_;
	}
	std::string refine() {
		return refine_;
	}
	std::string relax() {
		return relax_;
	}
	bool frag_files() {
		return frag_files_;
	}
	Size num_to_dock() {
		return num_to_dock_;
	}
	bool dump_initial_results() {
		return dump_initial_results_;
	}

	// setters
	void start_res( Size i ) {
		start_res_ = i;
	}
	void end_res( Size i ) {
		end_res_ = i;
	}
	void loop_length( Size i ) {
		loop_length_ = i;
	}
	void remodel( std::string i ) {
		remodel_ = i;
	}
	void intermedrelax( std::string i ) {
		intermedrelax_ = i;
	}
	void refine( std::string i ) {
		refine_ = i;
	}
	void relax( std::string i ) {
		relax_ = i;
	}
	void frag_files( bool i ) {
		frag_files_ = i;
	}
	void num_to_dock( Size i ) {
		num_to_dock_ = i;
	}
	void dump_initial_results( bool i ) {
		dump_initial_results_ = i;
	}

private:

	Size start_res_;
	Size end_res_;
	Size loop_length_;
	std::string remodel_;
	std::string intermedrelax_;
	std::string refine_;
	std::string relax_;
	bool frag_files_;
	Size num_to_dock_;
	bool dump_initial_results_;

};

}
}

#endif
