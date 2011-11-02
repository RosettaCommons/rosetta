// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/RepeatMover.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_moves_RepeatMover_hh
#define INCLUDED_protocols_moves_RepeatMover_hh

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/moves/RepeatMover.fwd.hh>

// Project headers
// AUTO-REMOVED #include <core/types.hh>

#include <core/pose/Pose.fwd.hh>

// AUTO-REMOVED #include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.fwd.hh>

// ObjexxFCL Headers

// C++ Headers
#include <string>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

/// @brief A Mover that repeats an input Mover a user-specified number of times
///
/// Common Methods:
///     RepeatMover.apply
class RepeatMover : public Mover {
public:
	// default constructor (nmoves=1)
	RepeatMover();
	/// @brief Constructs a RepeatMover
	/// repeatmover = RepeatMover( mover_in , nmoves_in )
	///
	/// Mover    mover_in    /object defining what move to make
	/// int      nmoves_in   /how many times to apply mover_in
	RepeatMover( MoverOP mover_in, int nmoves_in );
	~RepeatMover();

	/// @brief Repeats the input Mover a specified number of times
	///
	/// example(s):
	///     repeatmover.apply(pose)
	/// See Also:
	///     MinMover
	///     SequenceMover
	///     ShearMover
	///     SmallMover
	///     TrialMover
	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	MoverOP mover_;
	int nmoves_;
};

} // moves
} // protocols


#endif
