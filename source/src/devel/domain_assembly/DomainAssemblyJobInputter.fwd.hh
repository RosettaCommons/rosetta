// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/domain_assembly/DomainAssemblyJobInputter.fwd.hh
/// @brief  Forward declaration of DomainAssemblyJobInputter which creates a fusion pose from two PDB files
/// @author Robert Lindner <rlindner@mpimf-heidelberg.mpg.de>

#ifndef INCLUDED_devel_domain_assembly_DomainAssemblyJobInputter_fwd_hh
#define INCLUDED_devel_domain_assembly_DomainAssemblyJobInputter_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace domain_assembly {

class DomainAssemblyJobInputter;
typedef utility::pointer::shared_ptr< DomainAssemblyJobInputter > DomainAssemblyJobInputterOP;
typedef utility::pointer::shared_ptr< DomainAssemblyJobInputter const > DomainAssemblyJobInputterCOP;

}// domain_assembly
}// devel

#endif //INCLUDED_devel_domain_assembly_DomainAssemblyJobInputter_fwd_hh
