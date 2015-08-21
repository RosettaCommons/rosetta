// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/methods/FreeDOF_Options.fwd.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_methods_FreeDOF_Options_FWD_HH
#define INCLUDED_core_scoring_methods_FreeDOF_Options_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace methods {

class FreeDOF_Options;
typedef utility::pointer::shared_ptr< FreeDOF_Options > FreeDOF_OptionsOP;
typedef utility::pointer::shared_ptr< FreeDOF_Options const > FreeDOF_OptionsCOP;

} //methods
} //scoring
} //core

#endif
