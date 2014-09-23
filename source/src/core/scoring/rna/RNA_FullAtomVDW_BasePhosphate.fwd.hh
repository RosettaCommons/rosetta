// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rna/RNA_FullAtomVDW_BasePhosphate.fwd.hh
/// @author Parin Sripakdeevong, Rhiju Das (rhiju@stanford.edu)


#ifndef INCLUDED_core_scoring_rna_RNA_FullAtomVDW_BasePhosphate_FWD_HH
#define INCLUDED_core_scoring_rna_RNA_FullAtomVDW_BasePhosphate_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace rna {

class RNA_FullAtomVDW_BasePhosphate;

typedef utility::pointer::shared_ptr< RNA_FullAtomVDW_BasePhosphate > RNA_FullAtomVDW_BasePhosphateOP;
typedef utility::pointer::shared_ptr< RNA_FullAtomVDW_BasePhosphate const > RNA_FullAtomVDW_BasePhosphateCOP;

} //rna
} //scoring
} //core


#endif
