// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/moves/MoverApplyingMover.hh
/// @brief A MoverApplyingMover is nothing more than a mover that applies another (variable) mover.
///
/// @author Justin R. Porter

#ifndef INCLUDED_protocols_moves_MoverApplyingMover_hh
#define INCLUDED_protocols_moves_MoverApplyingMover_hh

// Unit Headers
#include <protocols/moves/MoverApplyingMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Package headers

// Project headers

// ObjexxFCL Headers

// Utility Headers

// C++ Headers

namespace protocols {
namespace moves {

class MoverApplyingMover : public Mover {
	typedef Mover Parent;
public:
	MoverApplyingMover( std::string const& name ) : Parent( name ) {}
	MoverApplyingMover( MoverApplyingMover const & ) = default;

	~MoverApplyingMover() override = default;

	virtual void set_mover( MoverOP ) = 0;

	virtual MoverOP mover() const = 0;

}; // end MoverApplyingMover base class

} // moves
} // protocols

#endif //INCLUDED_protocols_moves_MoverApplyingMover_HH
