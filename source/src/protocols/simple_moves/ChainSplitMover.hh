// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/ChainSplitMover.hh
/// @brief Declaration of ChainSplitMover which splits a pose into two chains at a given cutpoint


#ifndef INCLUDED_protocols_simple_moves_ChainSplitMover_hh
#define INCLUDED_protocols_simple_moves_ChainSplitMover_hh

// Unit Headers
#include <protocols/simple_moves/ChainSplitMover.fwd.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>

namespace protocols {
namespace simple_moves {

/// @brief ChainSplitMover splits a pose into two chains given a cutpoint. This may be necessary for evaluating a domain interface
/// using classes such as the InterfaceAnalyzerMover
class ChainSplitMover : public protocols::moves::Mover
{
public:

	ChainSplitMover();
	ChainSplitMover( core::Size cutpoint );
	ChainSplitMover( ChainSplitMover const & src );
	~ChainSplitMover();
	virtual void apply( core::pose::Pose & pose );
	virtual protocols::moves::MoverOP clone() const;
	virtual protocols::moves::MoverOP fresh_instance() const;
	virtual std::string get_name() const;

	/// @brief Retrieve the cutpoint at which the pose will be split
	// Undefined, commenting out to fix PyRosetta build  core::Size cutpoint( ) const ;

	/// @brief Set the cutpoint at which the pose will be split. Pose coordinates!
	// Undefined, commenting out to fix PyRosetta build  void cutpoint( core::Size splitpoint );
private:
	core::Size cutpoint_;
};

} // simple_moves
} // protocols


#endif
