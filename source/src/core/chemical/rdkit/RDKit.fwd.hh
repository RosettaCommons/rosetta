// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/rdkit/RDKit.fwd.hh
/// @brief Forward declarations which make working with RDKit nicer
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_rdkit_RDKit_fwd_hh
#define INCLUDED_core_chemical_rdkit_RDKit_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// C++ headers

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>

namespace ForceFields {

class ForceField;

}

namespace RDKit {

class RWMol;
class ROMol;
class Atom;
class ChemicalReaction;

// These owning pointers works because we're now using extrinsic owning pointers (Boost/C++11)
//
// Be careful about passing the pointers from this OP back to RDKit
// Some functions take ownership, which is problematic as OPs can never give up ownership

typedef ::RDKit::RWMOL_SPTR  RWMolOP;
typedef ::RDKit::ROMOL_SPTR  ROMolOP;
typedef utility::pointer::shared_ptr< Atom > AtomOP;

typedef utility::pointer::shared_ptr< ChemicalReaction >  ChemicalReactionOP;
typedef utility::pointer::shared_ptr< ChemicalReaction const >  ChemicalReactionCOP;

typedef utility::pointer::shared_ptr< ::ForceFields::ForceField > ForceFieldOP;

}

#endif
