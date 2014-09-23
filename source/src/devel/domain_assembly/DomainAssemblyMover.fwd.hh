// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/domain_assembly/DomainAssemblyMover
/// @brief


#ifndef INCLUDED_devel_domain_assembly_DomainAssemblyMover_FWD_HH
#define INCLUDED_devel_domain_assembly_DomainAssemblyMover_FWD_HH

// Utility Headers
#include <utility/pointer/owning_ptr.hh>


namespace devel{
namespace domain_assembly{

class DomainAssemblyMover;

typedef utility::pointer::shared_ptr< DomainAssemblyMover > DomainAssemblyMoverOP;
typedef utility::pointer::shared_ptr< DomainAssemblyMover const > DomainAssemblyMoverCOP;

}
}


#endif
