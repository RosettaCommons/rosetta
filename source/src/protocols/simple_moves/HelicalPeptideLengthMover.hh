// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/HelicalPeptideLengthMover.hh
/// @brief Definition of classes for changing the length of a helical peptide that's connected to the rest of the pose by a jump
/// @author Jacob Corn

#ifndef INCLUDED_protocols_simple_moves_HelicalPeptideLengthMover_hh
#define INCLUDED_protocols_simple_moves_HelicalPeptideLengthMover_hh

// Unit Headers
#include <protocols/simple_moves/HelicalPeptideLengthMover.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Package headers
#include <protocols/moves/Mover.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>

// ObjexxFCL Headers

// C++ Headers

// Utility Headers

namespace protocols {
namespace simple_moves {

enum residue_addition_type {
	append=1,
	prepend,
	insert
};


class HelicalPeptideLengthMover : public protocols::moves::Mover {
public:
	HelicalPeptideLengthMover( core::Size const seqpos, std::string const resname, residue_addition_type const addition_type );

	virtual void apply( core::pose::Pose & pose );
	virtual protocols::moves::MoverOP clone() const {
		return( new protocols::simple_moves::HelicalPeptideLengthMover( seqpos_, resname_, addition_type_ ) );
	}

private:
	// private variables
	core::Size seqpos_;
	std::string resname_;
	residue_addition_type addition_type_;

	// private fxns
	core::conformation::ResidueCOP create_residue_from_resname_( std::string const resname );
	void add_residue_( core::pose::Pose & pose, core::Size const seqpos, core::conformation::ResidueCOP residue, residue_addition_type const append_prepend );

}; // HelicalPeptideLengthMover

} // moves
} // protocols


#endif
