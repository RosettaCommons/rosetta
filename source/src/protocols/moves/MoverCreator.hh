// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/moves/MoverCreator.hh
/// @brief  Base class for MoverCreators for the Mover load-time factory registration scheme
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_protocols_moves_MoverCreator_hh
#define INCLUDED_protocols_moves_MoverCreator_hh

// Unit Headers
#include <protocols/moves/Mover.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

// c++ headers
#include <string>

namespace protocols {
namespace moves {

/// @brief Abstract base class for a Mover factory; the Creator class is responsible for
/// creating a particular mover class.
class MoverCreator : public utility::pointer::ReferenceCount
{
public:
	MoverCreator();
	virtual ~MoverCreator();

	/// @brief Return a new mover.
	virtual MoverOP create_mover() const = 0;

	/// @brief Return the tag name associated with this factory.
	virtual std::string keyname() const = 0;
};

typedef utility::pointer::shared_ptr< MoverCreator > MoverCreatorOP;
typedef utility::pointer::shared_ptr< MoverCreator const > MoverCreatorCOP;

} //namespace moves
} //namespace protocols

#endif
