// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/symmetry/SymmetricConformation.fwd.hh
/// @author Ingemar Andre

#ifndef INCLUDED_core_conformation_symmetry_SymmetricConformation_fwd_hh
#define INCLUDED_core_conformation_symmetry_SymmetricConformation_fwd_hh


#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace conformation {
namespace symmetry {

class SymmetricConformation;
typedef utility::pointer::shared_ptr< SymmetricConformation > SymmetricConformationOP;
typedef utility::pointer::shared_ptr< SymmetricConformation const > SymmetricConformationCOP;

} // symmetry
} // conformation
} // core

#ifdef USEBOOSTSERIALIZE
#include <boost/serialization/base_object.hpp>
#endif

#endif
