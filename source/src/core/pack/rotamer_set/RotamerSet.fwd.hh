// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/RotamerSet.fwd.hh
/// @brief  Residue set class forward declarations
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_pack_rotamer_set_RotamerSet_fwd_hh
#define INCLUDED_core_pack_rotamer_set_RotamerSet_fwd_hh

#include <core/conformation/Residue.fwd.hh>

// utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace core {
namespace pack {
namespace rotamer_set {

class RotamerSet;

typedef utility::vector1< conformation::ResidueOP > Rotamers;
typedef utility::pointer::shared_ptr< RotamerSet > RotamerSetOP;
typedef utility::pointer::shared_ptr< RotamerSet const > RotamerSetCOP;


} // namespace rotamer_set
} // namespace pack
} // namespace core


#endif // INCLUDED_core_pack_RotamerSet_RotamerSet_fwd_HH
