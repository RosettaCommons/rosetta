// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FilterMover
/// @brief performs MoverOP and then apply FilterOP
/// @author Nobuyasu Koga ( nobuyasau@uw.edu )

#ifndef INCLUDED_protocols_moves_FilterMover_hh
#define INCLUDED_protocols_moves_FilterMover_hh

// Unit header
#include <protocols/moves/FilterMover.fwd.hh>

// Project headers
#include <core/types.hh>
#include <protocols/moves/MoverStatus.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {


class FilterMover : public Mover {

public:

	typedef core::Size Size;
	typedef core::pose::Pose Pose;
	typedef protocols::moves::MoverOP MoverOP;
	typedef protocols::filters::FilterOP FilterOP;

public:

	// default constructor
	FilterMover();

	// @brief constructor with arguments
	// mover_status when filter failed
	FilterMover( MoverOP const & my_mover, FilterOP const & my_filter, Size const max_tries,
		MoverStatus const mover_status = FAIL_DO_NOT_RETRY );

	~FilterMover();

	// @brief main operation
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

	// @brief add filter
	void add_filter( FilterOP const & my_filter );

	/// @brief set mover
	void set_mover( MoverOP const & my_mover );

	/// @brief set MoverStatus when the mover failed
	//void moverstatus_whenfail( MoverStatus const & mover_status );

	/// @brief set maximum tries of making a move with my_mover
	void max_tries( Size const mt );

private:

	MoverOP my_mover_;
	FilterOP my_filter_;
	Size max_tries_;
	MoverStatus ms_whenfail_;

};


} // moves
} // protocols


#endif
