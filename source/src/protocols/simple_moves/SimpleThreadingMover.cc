// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/SimpleThreadingMover.cc
/// @brief Very Simple class for threading a regional sequence onto a structure
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/simple_moves/SimpleThreadingMoverCreator.hh>
#include <protocols/simple_moves/SimpleThreadingMover.hh>

#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/chemical/AA.hh>
#include <core/select/util.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/loops/loops_main.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

static THREAD_LOCAL basic::Tracer TR("SimpleThreadingMover");

namespace protocols {
namespace simple_moves {
using namespace core::select;

SimpleThreadingMover::SimpleThreadingMover():
	protocols::moves::Mover("SimpleThreadingMover")
{
	set_defaults();
}

SimpleThreadingMover::SimpleThreadingMover(std::string thread_sequence, core::Size start_position):
	protocols::moves::Mover("SimpleThreadingMover"),
	start_position_(start_position)


{
	set_defaults();
	thread_sequence_ = thread_sequence;
}

SimpleThreadingMover::~SimpleThreadingMover()= default;

void
SimpleThreadingMover::set_defaults(){
	pack_neighbors_ = false;
	parsed_position_ = "NA";
	thread_sequence_ = "NA";
	scorefxn_ = core::scoring::get_score_function();
	skip_unknown_mutant_ = false;
	neighbor_dis_ = 6.0;
	pack_rounds_ = 5.0;

}


void
SimpleThreadingMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	pack_neighbors_ = tag->getOption< bool >("pack_neighbors", pack_neighbors_);
	neighbor_dis_ = tag->getOption< core::Real >("neighbor_dis", neighbor_dis_);

	parsed_position_ = tag->getOption< std::string >("start_position");

	thread_sequence_ = tag->getOption< std::string >("thread_sequence", thread_sequence_);
	if ( tag->hasOption("scorefxn") ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data );
	}

	skip_unknown_mutant_ = tag->getOption< bool >("skip_unknown_mutant", skip_unknown_mutant_);
	pack_rounds_ = tag->getOption< core::Size >("pack_rounds", pack_rounds_);

}



SimpleThreadingMover::SimpleThreadingMover(SimpleThreadingMover const & )= default;

protocols::moves::MoverOP
SimpleThreadingMover::clone() const{
	return protocols::moves::MoverOP( new SimpleThreadingMover(*this) );

}

//SimpleThreadingMover & operator=( SimpleThreadingMover const & src){
// return SimpleThreadingMover(src);

//}

moves::MoverOP
SimpleThreadingMover::fresh_instance() const {
	return protocols::moves::MoverOP( new SimpleThreadingMover);

}

std::string
SimpleThreadingMover::get_name() const {
	return "SimpleThreadingMover";

}

void
SimpleThreadingMover::set_sequence(std::string thread_sequence, core::Size start_position){
	start_position_ = start_position;
	thread_sequence_ = thread_sequence;

}

void
SimpleThreadingMover::set_pack_neighbors(bool pack_neighbors){
	pack_neighbors_ = pack_neighbors;
}

void
SimpleThreadingMover::set_pack_rounds(core::Size pack_rounds){
	pack_rounds_ = pack_rounds;
}

void
SimpleThreadingMover::set_neighbor_distance(core::Real neighbor_dis){
	neighbor_dis_ = neighbor_dis;
}

void
SimpleThreadingMover::set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn){
	scorefxn_ = scorefxn;
}

bool
SimpleThreadingMover::get_pack_neighbors() const{
	return pack_neighbors_;

}

core::Real
SimpleThreadingMover::get_neighbor_distance() const {
	return neighbor_dis_;

}

void
SimpleThreadingMover::apply(core::pose::Pose& pose){
	using namespace core::pack::task;
	using namespace core::scoring;
	using namespace core::pack::task::operation;



	//This could have just as easily have been a task op.

	if ( thread_sequence_ == "NA" ) {
		utility_exit_with_message("Sequence not set for threading.  Cannot continue.");
	}

	if ( parsed_position_ !=  "NA" ) {
		start_position_ = core::pose::parse_resnum(parsed_position_, pose);
	}

	TR << "Threading Sequence :"<<thread_sequence_<<":"<<std::endl;
	TaskFactoryOP tf = TaskFactoryOP( new TaskFactory());
	tf->push_back(TaskOperationCOP(new InitializeFromCommandline()) );
	PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);


	utility::vector1< bool > mutant_resnums( pose.size(), false);
	std::map< core::Size, core::chemical::AA > mutations;
	std::map< core::Size, core::chemical::AA >::iterator it;

	core::Size seq_position = 0;
	for ( core::Size resnum = start_position_; resnum <= start_position_+thread_sequence_.size() -1; ++resnum ) {

		mutant_resnums[ resnum ] = true;

		if ( core::chemical::oneletter_code_specifies_aa(thread_sequence_[ seq_position ]) ) {
			mutations[ resnum ] = core::chemical::aa_from_oneletter_code(thread_sequence_[ seq_position ]);
		} else if ( thread_sequence_[ seq_position ] == '-' ) {
			seq_position += 1;
			continue;
		} else {
			TR << "Amino Acid not understood: "<< thread_sequence_[ seq_position ] << std::endl;
			if ( skip_unknown_mutant_ ) continue;
			else utility_exit_with_message("SimpleThreadingMover: Unknown Amino Acid at position "+utility::to_string( seq_position + 1)+" Pass skip_unknown_mutant to skip this instead of fail.");
		}
		seq_position += 1;
	}

	//Enable positions to design - but only into the residues in list.
	for ( it = mutations.begin(); it != mutations.end(); ++it ) {
		utility::vector1< bool > allowed_aminos(20, false);
		allowed_aminos[ it->second ] = true;
		task->nonconst_residue_task(it->first).restrict_absent_canonical_aas(allowed_aminos);
	}

	core::pack::task::operation::PreventRepacking turn_off_packing;
	core::pack::task::operation::RestrictResidueToRepacking turn_off_design;

	//Select Set the packer to turn off design everywhere but our residues.
	for ( core::Size resnum=1; resnum <= pose.size(); ++resnum ) {
		if ( ! mutant_resnums[ resnum ] ) {
			turn_off_design.include_residue( resnum ); //Turn all design off except residues we are forcing.
			if ( ! pack_neighbors_ ) {
				turn_off_packing.include_residue( resnum );
			}
		}
	}

	scorefxn_->score(pose); //Segfault Protection.

	//If we pack neighbors
	if ( pack_neighbors_ ) {
		utility::vector1< bool > mutant_resnums_and_neighbors = mutant_resnums;

		core::select::fill_neighbor_residues(pose, mutant_resnums_and_neighbors, neighbor_dis_);
		for ( core::Size resnum = 1; resnum <= pose.size(); ++resnum ) {
			if ( ! mutant_resnums_and_neighbors[ resnum ] ) {
				turn_off_packing.include_residue( resnum );
			}
		}
	}

	turn_off_design.apply(pose, *task);
	turn_off_packing.apply(pose, *task);

	protocols::simple_moves::PackRotamersMover packer = PackRotamersMover(scorefxn_, task, pack_rounds_);
	packer.apply(pose);
	TR << "Complete" <<std::endl;

}


////////////// Creator /////////

protocols::moves::MoverOP
SimpleThreadingMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SimpleThreadingMover );
}

std::string
SimpleThreadingMoverCreator::keyname() const {
	return SimpleThreadingMoverCreator::mover_name();
}

std::string
SimpleThreadingMoverCreator::mover_name(){
	return "SimpleThreadingMover";
}


}//simple_moves
}//protocols

