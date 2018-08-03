// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/chemically_conjugated_docking/UBQ_GTPaseMover.fwd.hh
/// @brief A mover typically used for binding ubiquitin to a substrate protein.
/// @author Hope Anderson (hdanderson)

#ifndef INCLUDED_protocols_chemically_conjugated_docking_UBQ_GTPaseMover_fwd_hh
#define INCLUDED_protocols_chemically_conjugated_docking_UBQ_GTPaseMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace chemically_conjugated_docking {

class UBQ_GTPaseMover;

typedef utility::pointer::shared_ptr< UBQ_GTPaseMover > UBQ_GTPaseMoverOP;
typedef utility::pointer::shared_ptr< UBQ_GTPaseMover const > UBQ_GTPaseMoverCOP;

} //protocols
} //chemically_conjugated_docking

#endif //INCLUDED_protocols_chemically_conjugated_docking_UBQ_GTPaseMover_fwd_hh
