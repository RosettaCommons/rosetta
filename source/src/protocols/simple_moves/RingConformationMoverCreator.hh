// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/simple_moves/RingConformationMoverCreator.hh
/// @brief   Method declarations for RingConformationMoverCreator.
/// @author  Labonte <JWLabonte@jhu.edu>

#ifndef INCLUDED_protocols_simple_moves_RingConformationMoverCreator_HH
#define INCLUDED_protocols_simple_moves_RingConformationMoverCreator_HH

// Project headers
#include <protocols/moves/MoverCreator.hh>


namespace protocols {
namespace simple_moves {

/// @brief  MoverCreator allowing the MoverFactory to create a RingConformationMover
class RingConformationMoverCreator: public protocols::moves::MoverCreator {
public:
	/// @brief  Return an up-casted owning pointer (MoverOP) to the mover.
	virtual protocols::moves::MoverOP create_mover() const;

	/// @brief  Return the string identifier for the associated Mover (RingConformationMover).
	virtual std::string keyname() const;

	/// @brief  Static method that returns the keyname for performance reasons.
	static std::string mover_name();
};

}  // namespace simple_moves
}  // namespace protocols

#endif  // INCLUDED_simple_moves_protocols_RingConformationMoverCreator_HH
