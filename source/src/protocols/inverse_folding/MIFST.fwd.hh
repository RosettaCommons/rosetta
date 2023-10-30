// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/inverse_folding/MIFST.fwd.hh
/// @brief A class for using the MIF-ST model developed by Yang et al.
/// @author Moritz Ertelt (moritz.ertelt@gmail.com)

#ifndef INCLUDED_protocols_inverse_folding_MIFST_fwd_hh
#define INCLUDED_protocols_inverse_folding_MIFST_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace inverse_folding {

class MIFST;

using MIFSTOP = utility::pointer::shared_ptr< MIFST >;
using MIFSTCOP = utility::pointer::shared_ptr< MIFST const >;

} //inverse_folding
} //protocols

#endif //INCLUDED_protocols_inverse_folding_MIFST_fwd_hh
