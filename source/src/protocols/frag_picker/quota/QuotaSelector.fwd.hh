// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/quota/QuotaSelector.fwd.hh
/// @brief
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_quota_QuotaSelector_fwd_hh
#define INCLUDED_protocols_frag_picker_quota_QuotaSelector_fwd_hh

// utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace frag_picker {
namespace quota {

class QuotaSelector;

typedef utility::pointer::shared_ptr<QuotaSelector> QuotaSelectorOP;
typedef utility::pointer::shared_ptr<QuotaSelector const> QuotaSelectorCOP;

} // quota
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_quota_QuotaSelector_FWD_HH */
