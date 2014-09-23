// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/fldsgn/CircularPermutation.fwd.hh
/// @brief  forward declaration for CircularPermutation
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#ifndef INCLUDED_protocols_fldsgn_CircularPermutation_fwd_hh
#define INCLUDED_protocols_fldsgn_CircularPermutation_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>


namespace protocols {
namespace fldsgn {

/// @brief forward declaration for CircularPermutation
class CircularPermutation;


/// @brief CircularPermutation owning pointer
typedef utility::pointer::shared_ptr< CircularPermutation > CircularPermutation_OP;


/// @brief CircularPermutation const owning pointer
typedef utility::pointer::shared_ptr< CircularPermutation const > CircularPermutation_COP;


/// @brief CircularPermutation access pointer
typedef utility::pointer::weak_ptr< CircularPermutation > CircularPermutation_AP;


/// @brief CircularPermutation const access pointer
typedef utility::pointer::weak_ptr< CircularPermutation const > CircularPermutation_CAP;


} // namespace fldsgn
} // namespace protocols


#endif /* INCLUDED_protocols_fldsgn_CircularPermutation_FWD_HH */
