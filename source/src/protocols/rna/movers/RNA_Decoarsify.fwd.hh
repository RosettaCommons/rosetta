// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/movers/RNA_Decoarsify.fwd.hh
/// @brief Make a coarse RNA pose into a fullatom representation.
/// @author Andy Watkins (andy.watkins2@gmail.com)

#ifndef INCLUDED_protocols_rna_movers_RNA_Decoarsify_fwd_hh
#define INCLUDED_protocols_rna_movers_RNA_Decoarsify_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace rna {
namespace movers {

class RNA_Decoarsify;

using RNA_DecoarsifyOP = utility::pointer::shared_ptr< RNA_Decoarsify >;
using RNA_DecoarsifyCOP = utility::pointer::shared_ptr< RNA_Decoarsify const >;

} //movers
} //rna
} //protocols

#endif //INCLUDED_protocols_rna_movers_RNA_Decoarsify_fwd_hh
