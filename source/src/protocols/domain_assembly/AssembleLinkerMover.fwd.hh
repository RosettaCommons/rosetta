/// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/domain_assembly/AssembleLinkerMover.fwd.hh
/// @brief  A mover for assembling domains
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_domain_assembly_AssembleLinkerMover_fwd_hh
#define INCLUDED_protocols_domain_assembly_AssembleLinkerMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace domain_assembly {

class AssembleLinkerMover;
typedef utility::pointer::shared_ptr< AssembleLinkerMover > AssembleLinkerMoverOP;
typedef utility::pointer::shared_ptr< AssembleLinkerMover const > AssembleLinkerMoverCOP;


}//namespace domain_assembly
}//namespace protocols

#endif

