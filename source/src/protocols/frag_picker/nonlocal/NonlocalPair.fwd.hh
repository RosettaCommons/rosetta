// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/nonlocal/NonlocalPair.fwd.hh
/// @brief  forward declaration for NonlocalPair
/// @author David E. Kim (dekim@u.washington.edu)

#ifndef INCLUDED_protocols_frag_picker_nonlocal_NonlocalPair_fwd_hh
#define INCLUDED_protocols_frag_picker_nonlocal_NonlocalPair_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace frag_picker {
namespace nonlocal {

/// @brief forward declaration for NonlocalPair
class NonlocalPair;

typedef utility::pointer::shared_ptr<NonlocalPair> NonlocalPairOP;
typedef utility::pointer::shared_ptr<NonlocalPair const> NonlocalPairCOP;

} // nonlocal
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_nonlocal_NonlocalPair_FWD_HH */
