// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/rotamer_set/bb_independent_rotamers.cc
/// @brief  build rotamers for residue type given neutral beckbone
/// @author Florian Ricter (floric@uw.edu)

#ifndef INCLUDED_core_pack_rotamer_set_bb_independent_rotamers_hh
#define INCLUDED_core_pack_rotamer_set_bb_independent_rotamers_hh

#include <core/conformation/Residue.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <utility/vector1.hh>

namespace core {
namespace pack {
namespace rotamer_set {


/// @brief a crude function to spit out a list of rotamers
/// given the residue type only, independent of backbone
/// currently there is no proper way of doing this, since
/// the Dunbrack bbind library is not implemented in rosetta.
/// this function tries to circumvent that by constructing
/// a one residue pose and then using the regular dunbrack
/// library, which will use neutral phi/psi for the only
/// residue in the pose
/// the bool ignore_cmdline can be used if someone only
/// wants base inverse rotamers but use the full set
/// in packing
utility::vector1< core::conformation::ResidueCOP >
bb_independent_rotamers(
	core::chemical::ResidueTypeCOP rot_restype,
	bool ignore_cmdline = false
);


} // namespace
} // namespace
} // namespace


#endif // INCLUDED_core_pack_rotamer_set_bb_independent_rotamers_hh
