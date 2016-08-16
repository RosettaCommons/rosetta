// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/RepeatMover.hh
/// @brief Declarations and simple accessor/mutator definitions for RepeatMover

#ifndef INCLUDED_protocols_moves_RepeatMover_hh
#define INCLUDED_protocols_moves_RepeatMover_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/RepeatMover.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>


namespace protocols {
namespace moves {

/// @brief A Mover that repeats an input Mover a user-specified number of times
///
/// @details Common Methods:
///     RepeatMover.apply
class RepeatMover : public Mover {
public:
	/// @brief Empty constructor (nmoves=1)
	RepeatMover();

	/// @brief Constructs a RepeatMover
	RepeatMover(MoverOP mover_in, int nmoves_in);

	/// @brief Copy constructor
	RepeatMover(RepeatMover const & object_to_copy);

	~RepeatMover();

	/// @brief Repeats the input Mover a specified number of times
	virtual void apply( core::pose::Pose & pose );

	virtual std::string get_name() const;
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;

	core::Size get_nmoves() const;
	std::string get_mover() const;

private:
	MoverOP mover_;
	int nmoves_;
};


std::ostream &operator<< (std::ostream &os, RepeatMover const &mover);

} // moves
} // protocols

#endif
