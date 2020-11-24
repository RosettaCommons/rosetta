// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide/crosslinker/1_4_BBMB_Helper.fwd.hh
/// @brief A derived class of the CrosslinkerMoverHelper base class, used to set up
/// the 1,3,5-tris(bromomethyl)benzene (BBMB) cross-linker.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_cyclic_peptide_crosslinker_1_4_BBMB_Helper_fwd_hh
#define INCLUDED_protocols_cyclic_peptide_crosslinker_1_4_BBMB_Helper_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace cyclic_peptide {
namespace crosslinker {

class One_Four_BBMB_Helper;

typedef utility::pointer::shared_ptr< One_Four_BBMB_Helper > One_Four_BBMB_HelperOP;
typedef utility::pointer::shared_ptr< One_Four_BBMB_Helper const > One_Four_BBMB_HelperCOP;

} //crosslinker
} //protocols
} //cyclic_peptide

#endif //INCLUDED_protocols_cyclic_peptide_crosslinker_1_4_BBMB_Helper_fwd_hh
