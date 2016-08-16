// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file FilterMover
/// @brief apply filter after mover
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

// Unit headers
#include <protocols/moves/FilterMover.hh>

// Project Headers
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.FilterMover" );


#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

FilterMover::FilterMover():
	max_tries_( 1 ),
	ms_whenfail_( FAIL_DO_NOT_RETRY )
{}

FilterMover::FilterMover(
	MoverOP const & my_mover,
	FilterOP const & my_filter,
	Size const max_tries,
	MoverStatus const mover_status
):Mover( my_mover->type() ),
	my_mover_( my_mover ),
	my_filter_( my_filter ),
	max_tries_( max_tries ),
	ms_whenfail_( mover_status )
{}

FilterMover::~FilterMover(){}

void FilterMover::add_filter( FilterOP const & my_filter ){
	my_filter_ = my_filter;
}

void FilterMover::set_mover( MoverOP const & my_mover ){
	my_mover_ = my_mover;
}

void FilterMover::max_tries( Size const mt ){
	max_tries_ = mt;
}

/// @details  Keep trying to make a move with my_mover until my_filter return true or max_tries is exceeded
/// @note  At the expense of a few additional pose copies we could save the pose with the best and use
/// that if we exceed max_tries

void FilterMover::apply( Pose & pose )
{

	Pose const start_pose( pose );

	Size ntries( 0 );

	while ( true ) { // keep looping until we succeed

		++ntries;
		TR << " ntries/max_tries = " << ntries << '/' << max_tries_ << std::endl;
		my_mover_->apply( pose );

		//pose.dump_pdb( "hoge.pdb" );
		if ( my_mover_->get_last_move_status() == MS_SUCCESS ) {

			bool filter_status = my_filter_->apply( pose );

			if ( filter_status ) {
				set_last_move_status( MS_SUCCESS );
				break;
			} else {
				if ( ntries >= max_tries_ ) {
					set_last_move_status( ms_whenfail_ );
					break;
				}
			} // filter_status

		} else {
			set_last_move_status( my_mover_->get_last_move_status() );
			break;
		}// my_mover_->get_last_move_status() == MS_SUCCESS

		pose = start_pose;

	}

} //apply

std::string
FilterMover::get_name() const {
	return "FilterMover";
}

}  // namespace moves
}  // namespace protocols
