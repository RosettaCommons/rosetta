// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/CyclicMover.hh
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_PROTOCOLS_MOVES_CYCLICMOVER_HH
#define INCLUDED_PROTOCOLS_MOVES_CYCLICMOVER_HH

// Unit header
#include <protocols/moves/CyclicMover.fwd.hh>

// C/C++ headers
#include <string>
#include <vector>

// Project headers
#include <core/pose/Pose.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace moves {

/// @detail A simple class for cycling between movers in consecutive calls to apply()
class CyclicMover : public Mover {
	typedef std::vector<MoverOP> Movers;

public:
	/// @brief Creates a new instance with no enqueued movers
	CyclicMover();

	/// @brief Enqueue the specified mover for execution
	void enqueue(MoverOP mover);

	// -- mover -- //
	std::string get_name() const;
	void apply(core::pose::Pose& pose);

	// -- jd2 -- //
	MoverOP clone() const;
	MoverOP fresh_instance() const;

private:
	/// @brief Tracks the number of calls to apply()
	long iterations_;

	/// @brief List of movers, which are executed in order in consecutive calls to apply()
	Movers movers_;
};

}  // namespace moves
}  // namespace protocols

#endif  // PROTOCOLS_MOVES_CYCLIC_MOVER_HH_
