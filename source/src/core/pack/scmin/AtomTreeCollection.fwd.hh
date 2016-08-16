// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/scmin/AtomTreeCollection.hh
/// @brief  Forward declaration for the classes holding sets of AtomTrees used during variuos packing+minimizing schemes
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_core_pack_scmin_AtomTreeCollection_fwd_hh
#define INCLUDED_core_pack_scmin_AtomTreeCollection_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace scmin {

class ResidueAtomTreeCollectionMomento;
typedef utility::pointer::shared_ptr< ResidueAtomTreeCollectionMomento > ResidueAtomTreeCollectionMomentoOP;
typedef utility::pointer::shared_ptr< ResidueAtomTreeCollectionMomento const > ResidueAtomTreeCollectionMomentoCOP;

class ResidueAtomTreeCollection;

typedef utility::pointer::shared_ptr< ResidueAtomTreeCollection > ResidueAtomTreeCollectionOP;
typedef utility::pointer::shared_ptr< ResidueAtomTreeCollection const > ResidueAtomTreeCollectionCOP;

class AtomTreeCollection;

typedef utility::pointer::shared_ptr< AtomTreeCollection > AtomTreeCollectionOP;
typedef utility::pointer::shared_ptr< AtomTreeCollection const > AtomTreeCollectionCOP;

} // namespace scmin
} // namespace pack
} // namespace core

#endif
