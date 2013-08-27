// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/RNA_TorsionScreener.fwd.hh
/// @brief Screener checking whether the rna torsions are resonable
/// @author Fang-Chieh Chou

#ifndef INCLUDED_protocols_swa_rna_RNA_TorsionScreener_fwd_HH
#define INCLUDED_protocols_swa_rna_RNA_TorsionScreener_fwd_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace swa {
namespace rna {

class RNA_TorsionScreener;
typedef utility::pointer::owning_ptr< RNA_TorsionScreener > RNA_TorsionScreenerOP;

}
}
}
#endif
