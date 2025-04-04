// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/rna/MC_RNA_OneJump.fwd.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_recces_sampler_rna_MC_RNA_OneJump_FWD_HH
#define INCLUDED_protocols_recces_sampler_rna_MC_RNA_OneJump_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace recces {
namespace sampler {
namespace rna {

class MC_RNA_OneJump;
typedef utility::pointer::shared_ptr< MC_RNA_OneJump > MC_RNA_OneJumpOP;
typedef utility::pointer::shared_ptr< MC_RNA_OneJump const > MC_RNA_OneJumpCOP;

} //rna
} //sampler
} //recces
} //protocols

#endif
