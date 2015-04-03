// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/simple_moves/RepulsiveOnlyMover.hh
/// @brief
/// @author

#ifndef INCLUDED_protocols_simple_moves_RepulsiveOnlyMover_HH
#define INCLUDED_protocols_simple_moves_RepulsiveOnlyMover_HH

// Unit Headers
#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/RepulsiveOnlyMover.fwd.hh>

// Project headers

#include <core/pose/Pose.fwd.hh>


// ObjexxFCL Headers

// C++ Headers
#include <string>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace simple_moves {

/// @brief this mover replace designated residues with glycines, and then assigns a "REPLONLY" residue type to those glycines. Those "REPLONLY" residues will make only repulsive contributions to the overall energy.
class RepulsiveOnlyMover : public protocols::moves::Mover {
public:
	// default constructor (nmoves=1)
	RepulsiveOnlyMover();
	~RepulsiveOnlyMover();//deconstructor

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
	void set_mutate_to_glycine( bool setting ) {
		mutate_to_glycine_ = setting;
	}

private:
	bool mutate_to_glycine_;
};

} // moves
} // protocols


#endif
