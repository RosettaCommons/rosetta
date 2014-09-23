// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/rosetta_scripts/PoseSelectorCreator.hh
/// @brief  Base class for PoseSelectors for the load-time factory registration scheme
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_moves_PoseSelectorCreator_hh
#define INCLUDED_protocols_moves_PoseSelectorCreator_hh

// Unit Headers
#include <protocols/rosetta_scripts/PoseSelector.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

#include <string>

namespace protocols {
namespace rosetta_scripts {

/// @brief Abstract base class for a PoseSelector factory; the Creator class is responsible for
/// creating a particular PoseSelector class.
class PoseSelectorCreator : public utility::pointer::ReferenceCount
{
public:
	PoseSelectorCreator();
	virtual ~PoseSelectorCreator();

	virtual protocols::rosetta_scripts::PoseSelectorOP create_selector() const = 0;
	virtual std::string keyname() const = 0;
};

typedef utility::pointer::shared_ptr< PoseSelectorCreator > PoseSelectorCreatorOP;
typedef utility::pointer::shared_ptr< PoseSelectorCreator const > PoseSelectorCreatorCOP;

} //namespace rosetta_scripts
} //namespace protocols

#endif
