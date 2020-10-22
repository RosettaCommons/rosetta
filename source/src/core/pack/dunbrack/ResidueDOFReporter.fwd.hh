// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/scoring/dunbrack/ResidueDOFReporter.fwd.hh
/// @brief   Class to measure the DOFs used by a RotamerLibrary
/// @author  Andrew Leaver-Fay

#ifndef INCLUDED_core_pack_dunbrack_ResidueDOFReporter_FWD_HH
#define INCLUDED_core_pack_dunbrack_ResidueDOFReporter_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace dunbrack {

class ResidueDOFReporter;

typedef utility::pointer::shared_ptr< ResidueDOFReporter > ResidueDOFReporterOP;
typedef utility::pointer::shared_ptr< ResidueDOFReporter const > ResidueDOFReporterCOP;

}
}
}


#endif
