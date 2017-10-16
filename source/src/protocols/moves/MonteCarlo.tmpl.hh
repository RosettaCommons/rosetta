// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/MonteCarlo.tmpl.hh
/// @brief  implentation for MonteCarlo template functions
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_moves_MonteCarlo_tmpl_hh
#define INCLUDED_protocols_moves_MonteCarlo_tmpl_hh

// unit headers
#include <protocols/moves/MonteCarlo.hh>

// project headers
#include <core/pose/Pose.hh>

namespace protocols {
namespace moves {


/// @brief attach observer to last accepted conformation
/// @tparam ConformationObserver any class implementing <tt> void attach_to( Conformation & ) </tt>
template< typename ConformationObserver >
void
MonteCarlo::attach_observer_to_last_accepted_conformation( ConformationObserver & obs ) {
	obs.attach_to( last_accepted_pose_->conformation() );
}


/// @brief attach observer to lowest score conformation
/// @tparam ConformationObserver any class implementing <tt> void attach_to( Conformation & ) </tt>
template< typename ConformationObserver >
void
MonteCarlo::attach_observer_to_lowest_score_conformation( ConformationObserver & obs ) {
	obs.attach_to( lowest_score_pose_->conformation() );
}


/// @brief attach observer to last accepted pose
/// @tparam PoseObserver any class implementing <tt> void attach_to( Pose & ) </tt>
template< typename PoseObserver >
void
MonteCarlo::attach_observer_to_last_accepted_pose( PoseObserver & obs ) {
	obs.attach_to( *last_accepted_pose_ );
}


/// @brief attach observer to lowest score pose
/// @tparam PoseObserver any class implementing <tt> void attach_to( Pose & ) </tt>
template< typename PoseObserver >
void
MonteCarlo::attach_observer_to_lowest_score_pose( PoseObserver & obs ) {
	obs.attach_to( *lowest_score_pose_ );
}


} // namespace moves
} // protocols


#endif /* INCLUDED_protocols_moves_MonteCarlo_TMPL_HH */
