// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/chemical/Atom.fwd.hh
/// @brief  Class declaration for chemical::Atom
/// @note   not to be confused with conformation::Atom
/// @author Phil Bradley


#ifndef INCLUDED_core_chemical_Atom_FWD_HH
#define INCLUDED_core_chemical_Atom_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/vector1.fwd.hh>


namespace core {
namespace chemical {

/// @brief  A class for defining atom parameters, known as atom_types, with properties specific to the atom.
class Atom;

typedef  utility::pointer::shared_ptr< Atom >  AtomOP;
typedef  utility::pointer::shared_ptr< Atom const >  AtomCOP;
typedef  utility::pointer::weak_ptr< Atom >  AtomAP;
typedef  utility::vector1< AtomOP >  AtomOPs;
typedef  utility::vector1< AtomCOP >  AtomCOPs;
typedef  utility::vector1< AtomAP >  AtomAPs;

}  // chemical
}  // core

#endif  // INCLUDED_core_chemical_Atom_FWD_HH
