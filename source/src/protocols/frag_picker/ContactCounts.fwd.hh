// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/nonlocal/ContactCounts.fwd.hh
/// @brief  forward declaration for a ContactCounts
/// @author David E. Kim (dekim@u.washington.edu)

#ifndef INCLUDED_protocols_frag_picker_ContactCounts_fwd_hh
#define INCLUDED_protocols_frag_picker_ContactCounts_fwd_hh

// utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace frag_picker {

/// @brief forward declaration for ContactCounts
class ContactCounts;

typedef utility::pointer::shared_ptr<ContactCounts> ContactCountsOP;
typedef utility::pointer::shared_ptr<ContactCounts const> ContactCountsCOP;

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_ContactCounts_FWD_HH */
