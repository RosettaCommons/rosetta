// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/RotamerSetFactory.hh
/// @brief  Residue Set Factory class for symmetric packing
/// @author  Ingemar Andre

#ifndef INCLUDED_core_pack_rotamer_set_symmetry_SymmetricRotamerSetFactory_hh
#define INCLUDED_core_pack_rotamer_set_symmetry_SymmetricRotamerSetFactory_hh

// Unit headers
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSetFactory.fwd.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>

// Project Headers
#include <core/conformation/Residue.fwd.hh>


// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace rotamer_set {
namespace symmetry {

class SymmetricRotamerSetFactory : public RotamerSetFactory
{
public:
	virtual RotamerSetOP create_rotamer_set( conformation::Residue const & );

};

}
}
}
}

#endif
