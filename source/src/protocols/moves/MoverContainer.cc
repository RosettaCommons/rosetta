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
/// @details
/// @author Monica Berrondo

#include <protocols/moves/MoverContainerCreator.hh>
#include <protocols/moves/MoverContainer.hh>

// Rosetta Headers
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <protocols/moves/MoverStatus.hh>

// Random number generator
#include <numeric/random/random.hh>

#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/util.hh>


#include <string>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

static thread_local basic::Tracer tr( "protocols.moves.MoverContainer" );


namespace protocols {
namespace moves {


using namespace core;
using basic::T;
using basic::Error;
using basic::Warning;

/// @details Invoke clone() on each of the movers that are contained by this MoverContainer
/// to create deep copies.
MoverContainer::MoverContainer( MoverContainer const & source ) : protocols::moves::Mover()
{
	//remember that movers_ is a vector0
	movers_.resize( source.movers_.size() );
	for ( Size ii = 0; ii < source.movers_.size(); ++ii ) {
		movers_[ ii ] = source.movers_[ ii ]->clone();
	}
	weight_ = source.weight_;
}

void MoverContainer::add_mover( MoverOP mover_in , Real weight_in ) // do we need weights?
{
	movers_.push_back(mover_in);
	weight_.push_back(weight_in);
}

std::string MoverContainer::get_mover( core::Size index ) const {
	return movers_[index]->get_name();
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

// Sets the current_tag_ for both the container and the contained movers
void MoverContainer::set_current_tag( std::string const & new_tag ){
  Mover::set_current_tag( new_tag );
  for ( Size i=0; i<movers_.size(); ++i ) {
    movers_[i]->set_current_tag( new_tag );
  }

}

SequenceMover::SequenceMover( SequenceMover const & source ) :
	MoverContainer( source ),
	use_mover_status_( source.use_mover_status_ )
{}

MoverOP
SequenceMover::clone() const {
	return MoverOP( new SequenceMover( *this ) );
}

MoverOP
SequenceMover::fresh_instance() const {
	return MoverOP( new SequenceMover() );
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
		// set the mover status from the last mover applied, but do not act on the mover status
		for ( Size i=0; i<movers_.size(); ++i ) {
			movers_[i]->apply( pose );
			type( type()+movers_[i]->type() );
		}
		if ( movers_.size() > 0 ) {
			set_last_move_status( movers_[movers_.size()-1]->get_last_move_status() );
		}

	}
}

std::string
SequenceMover::get_name() const {
	return "SequenceMover";
}


std::string RandomMoverCreator::mover_name() {
  return "RandomMover";
}

std::string RandomMoverCreator::keyname() const {
  return mover_name();
}

protocols::moves::MoverOP RandomMoverCreator::create_mover() const {
  return protocols::moves::MoverOP( new RandomMover() );
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
		Real movechoice = numeric::random::rg().uniform()*weight_sum;
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

std::string RandomMover::get_name_individual_mover(core::Size index)
{
	std::string name = movers_[index]->get_name();
	return name;
}

RandomMover::RandomMover( RandomMover const & source ) :
	MoverContainer( source ),
	nmoves_( source.nmoves_ ),
	last_proposal_density_ratio_( source.last_proposal_density_ratio_ )
{}

MoverOP
RandomMover::clone() const {
	return MoverOP( new RandomMover( *this ) );
}

MoverOP
RandomMover::fresh_instance() const {
	return MoverOP( new RandomMover() );
}

std::string
RandomMover::get_name() const {
	return "RandomMover";
}

void
RandomMover::parse_my_tag( utility::tag::TagCOP tag,
                           basic::datacache::DataMap &,
                           protocols::filters::Filters_map const &,
                           protocols::moves::Movers_map const &movers,
                           core::pose::Pose const & ) {
	using namespace protocols::filters;

	utility::vector1<std::string> mover_names( utility::string_split( tag->getOption< std::string >("movers"), ',') );
	utility::vector1<std::string> mover_weights;
	if (tag->hasOption ("weights"))
		mover_weights = utility::string_split( tag->getOption< std::string >( "weights" ), ',');

	// make sure # movers matches # weights
	runtime_assert( mover_weights.size() == 0 || mover_weights.size() == mover_names.size() );

	nmoves_ = tag->getOption< core::Size >( "repeats",1 );

	for (core::Size i=1; i<=mover_names.size(); ++i) {
		protocols::moves::MoverOP mover = find_mover_or_die(mover_names[i], tag, movers);
		core::Real weight = (mover_weights.size() == 0)? 1.0 : std::atof( mover_weights[i].c_str() );
		add_mover( mover, weight );
	}
}


//ek added this function must be called AFTER apply
core::Real RandomMover::last_proposal_density_ratio(){
		return last_proposal_density_ratio_;
}

CycleMover::CycleMover( CycleMover const & source ) :
	MoverContainer( source ),
	next_move_( source.next_move_ )
{}

MoverOP
CycleMover::clone() const {
	return MoverOP( new CycleMover( *this ) );
}

MoverOP
CycleMover::fresh_instance() const {
	return MoverOP( new CycleMover() );
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

std::ostream &operator<< (std::ostream &os, MoverContainer const &mover)
{
	//moves::operator<<(os, mover);
	os << "Mover name: " << mover.get_name() << ", Mover type: " << mover.get_type() << ", Mover current tag: " << mover.get_current_tag() << std::endl;
	os << mover.size() << " movers are contained in the following order:" << std::endl;
	for ( Size i=0; i<mover.size(); ++i ) {
		os << "   Mover[" << i+1 << "]: " << mover.get_mover(i) << std::endl;
	}

	return os;
}

}  // namespace moves
}  // namespace protocols
