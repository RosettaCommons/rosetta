// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/struct_fragment/StructFragmentMover.fwd.hh
/// @brief  Creating fragments from supplied structures.
/// @author andreas scheck (graziano.bud@gmail.com), Correia's LPDI/EPFL


#ifndef INCLUDED_protocols_struct_fragment_StructFragmentMover_FWD_hh
#define INCLUDED_protocols_struct_fragment_StructFragmentMover_FWD_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace struct_fragment {

class StructFragmentMover;
typedef utility::pointer::shared_ptr< StructFragmentMover > StructFragmentMoverOP;
typedef utility::pointer::shared_ptr< StructFragmentMover const > StructFragmentMoverCOP;

}
}

#endif
