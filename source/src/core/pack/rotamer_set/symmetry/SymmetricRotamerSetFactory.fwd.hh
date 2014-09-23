// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file   core/pack/rotamer_set/RotamerSetFactory.fwd.hh
/// @brief  Residue Set Facotry class forward declaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_pack_rotamer_set_symmetry_SymmetricRotamerSetFactory_fwd_hh
#define INCLUDED_core_pack_rotamer_set_symmetry_SymmetricRotamerSetFactory_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace rotamer_set {
namespace symmetry {

class SymmetricRotamerSetFactory;

typedef utility::pointer::shared_ptr< SymmetricRotamerSetFactory > SymmetricRotamerSetFactoryOP;

}
}
}
}

#endif
