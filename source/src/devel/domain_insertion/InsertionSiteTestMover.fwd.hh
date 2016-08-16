// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/domain_insertion/InsertionSiteTestMover.fwd.hh
/// @brief  fwd hh file for InsertionSiteTestMover
/// @author Florian Richter, flosopher@gmail.com, november 2013

#ifndef INCLUDED_devel_domain_insertion_InsertionSiteTestMover_fwd_hh
#define INCLUDED_devel_domain_insertion_InsertionSiteTestMover_fwd_hh


// Utility headers
//#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>


namespace devel {
namespace domain_insertion {


class InsertionSiteTestMover;

typedef utility::pointer::shared_ptr< InsertionSiteTestMover > InsertionSiteTestMoverOP;
typedef utility::pointer::shared_ptr< InsertionSiteTestMover const > InsertionSiteTestMoverCOP;


}
}

#endif
