// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/denovo/movers/RNA_DeNovoProtocolMover.fwd.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu
/// @author Andy Watkins, amw579@stanford.edu


#ifndef INCLUDED_protocols_rna_denovo_movers_RNA_DeNovoProtocolMover_FWD_HH
#define INCLUDED_protocols_rna_denovo_movers_RNA_DeNovoProtocolMover_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace rna {
namespace denovo {
namespace movers {

class RNA_DeNovoProtocolMover;
typedef utility::pointer::shared_ptr< RNA_DeNovoProtocolMover > RNA_DeNovoProtocolMoverOP;
typedef utility::pointer::shared_ptr< RNA_DeNovoProtocolMover const > RNA_DeNovoProtocolMoverCOP;

} //movers
} //denovo
} //rna
} //protocols

#endif
