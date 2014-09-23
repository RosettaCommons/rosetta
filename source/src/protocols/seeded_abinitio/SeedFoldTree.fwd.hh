// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/seeded_abinitio/SeedFoldTree.fwd.hh
/// @author Eva-Maria Strauch (evas01@uw.edu)

#ifndef INCLUDED_PROTOCOLS_SEEDED_ABINITIO_SEEDFOLDTREE_FWD_HH
#define INCLUDED_PROTOCOLS_SEEDED_ABINITIO_SEEDFOLDTREE_FWD_HH

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace seeded_abinitio {

class SeedFoldTree;
typedef utility::pointer::shared_ptr<SeedFoldTree> SeedFoldTreeOP;
typedef utility::pointer::shared_ptr<SeedFoldTree const> SeedFoldTreeCOP;

}  // namespace seeded_abinitio
}  // namespace protocols

#endif  // INCLUDED_PROTOCOLS_SEEDED_ABINITIO_SEEDFOLDTREE_FWD_HH
