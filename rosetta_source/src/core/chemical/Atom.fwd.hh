// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/Atom.fwd.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_Atom_fwd_hh
#define INCLUDED_core_chemical_Atom_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

namespace core {
namespace chemical {

class Atom;

typedef  utility::pointer::owning_ptr< Atom >  AtomOP;
typedef  utility::pointer::owning_ptr< Atom const >  AtomCOP;
typedef  utility::vector1< AtomOP >  AtomOPs;
typedef  utility::vector1< AtomCOP >  AtomCOPs;

} // chemical
} // core



#endif // INCLUDED_core_chemical_Atom_FWD_HH
