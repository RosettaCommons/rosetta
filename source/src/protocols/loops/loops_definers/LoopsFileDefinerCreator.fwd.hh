// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/loops/loops_definers/LoopsFileDefinerCreator.fwd.hh
/// @brief  Forward Header for base class for LoopsFileDefinerCreator
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_loops_loops_definers_LoopsFileDefinerCreator_FWD_HH
#define INCLUDED_protocols_loops_loops_definers_LoopsFileDefinerCreator_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace loops {
namespace loops_definers {

/// @brief Abstract base class for a LoopsDefiner factory; the Creator class is responsible for
/// creating a particular LoopsDefiner class.
class LoopsFileDefinerCreator;

typedef utility::pointer::shared_ptr< LoopsFileDefinerCreator > LoopsFileDefinerCreatorOP;
typedef utility::pointer::shared_ptr< LoopsFileDefinerCreator const > LoopsFileDefinerCreatorCOP;

} //namespace
} //namespace
} //namespace

#endif
