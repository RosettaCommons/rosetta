// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MoverContainer.cc
/// @brief apply functions for classes of type MoverContainer
/// @detailed
/// @author Monica Berrondo

#include <protocols/moves/MoverContainer.hh>

// Rosetta Headers
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <basic/basic.hh>
#include <basic/Tracer.hh>

// AUTO-REMOVED #include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverStatus.hh>
// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/format.hh>

// Random number generator
#include <numeric/random/random.hh>

//
#include <string>

static basic::Tracer tr("protocols.moves.MoverContainer");


namespace protocols {
namespace moves {

static numeric::random::RandomGenerator RG(114);  // <- Mike's Magic number, do not change it (and dont try and use it anywhere else) !!! what a retarded system ... %-|

using namespace core;
using basic::T;
using basic::Error;
using basic::Warning;

void MoverContainer::add_mover( MoverOP mover_in , Real weight_in ) // do we need weights?
{
	movers_.push_back(mover_in);
	weight_.push_back(weight_in);
}

// Sets the input pose for both the container and the contained movers
// overriding this method fixes an annoying bug (barak)
	// TODO: does it make sense to cal this also from add_mover? I think not (barak)
void MoverContainer::set_input_pose( PoseCOP pose ){
	this->Mover::set_input_pose( pose );
	for ( Size i=0; i<movers_.size(); ++i ) {
		movers_[i]->set_input_pose( pose );
	}
}

// Sets the native pose for both the container and the contained movers
// overriding this method fixes an annoying bug (barak)
// TODO: does it make sense to call this also from add_mover? I think not (barak)
void MoverContainer::set_native_pose( PoseCOP pose ){
	this->Mover::set_native_pose( pose );
	for ( Size i=0; i<movers_.size(); ++i ) {
		movers_[i]->set_native_pose( pose );
	}
}

void SequenceMover::apply( core::pose::Pose & pose )
{

	using protocols::moves::MS_SUCCESS;
	using protocols::moves::FAIL_DO_NOT_RETRY;
	using protocols::moves::FAIL_BAD_INPUT;
	using protocols::moves::FAIL_RETRY;

	type("");
	if( use_mover_status_ ){

		Size i = 0;
		while( i < movers_.size() ){

			core::pose::Pose storepose( pose );
			movers_[i]->apply( pose );

			MoverStatus ms = movers_[i]->get_last_move_status();
			set_last_move_status( ms );
			if( ms == MS_SUCCESS ){
 				type( type()+movers_[i]->type() );
				i++;
			}else if( ms == FAIL_RETRY ){
				tr << "Mover failed. Same mover is performed again. " << std::endl;
				pose = storepose;
			}else if( ms == FAIL_DO_NOT_RETRY || ms == FAIL_BAD_INPUT ){
				tr << "Mover failed. Exit from SequenceMover." << std::endl;
				break;
			}
		}

	}else{

		for ( Size i=0; i<movers_.size(); ++i ) {
			movers_[i]->apply( pose );
			type( type()+movers_[i]->type() );
		}

	}
}

std::string
SequenceMover::get_name() const {
	return "SequenceMover";
}

void RandomMover::apply( core::pose::Pose & pose )
{
	Real weight_sum(0.0);
	size_t m;
	for(m=0;m< movers_.size(); m++){
		weight_sum += weight_[m];
	}
	type("");
	for(Size i=0;i< nmoves_; i++)
	{
		// choose a move
		Real movechoice = RG.uniform()*weight_sum;
		Real sum=0.0;
		for(m=0;m< movers_.size(); m++){
			if(movechoice < sum) break;
			sum += weight_[m];
		}
		m--;
		//		tr.Trace << "choose move " << m+1 << " of " << nr_moves() << std::endl;
		// apply the chosen move
		movers_[m]->apply( pose );
		type( type() + movers_[m]->type());
		last_proposal_density_ratio_ = movers_[m]->last_proposal_density_ratio();//ek
	}

}

std::string
RandomMover::get_name() const {
	return "RandomMover";
}

//ek added this function must be called AFTER apply
core::Real RandomMover::last_proposal_density_ratio(){
		return last_proposal_density_ratio_;
}

void CycleMover::apply( core::pose::Pose& pose )
{
	next_move_ %= movers_.size();
	movers_[ next_move_ ]->apply(pose);
	++next_move_;
}
	
	
//JQX added this to reset the next_move_
void CycleMover::reset_cycle_index()   
{
	next_move_ = 0;
}

std::string
CycleMover::get_name() const {
	return "CycleMover";
}

}  // namespace moves
}  // namespace protocols
