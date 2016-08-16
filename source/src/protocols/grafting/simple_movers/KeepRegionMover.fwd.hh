// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/grafting/simple_movers/KeepRegionMover
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_grafting_simple_movers_KeepRegionMover_fwd_hh
#define INCLUDED_protocols_grafting_simple_movers_KeepRegionMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace grafting {
namespace simple_movers {

class KeepRegionMover;

typedef utility::pointer::shared_ptr< KeepRegionMover > KeepRegionMoverOP;
typedef utility::pointer::shared_ptr< KeepRegionMover const > KeepRegionMoverCOP;

}
}
}

#endif


