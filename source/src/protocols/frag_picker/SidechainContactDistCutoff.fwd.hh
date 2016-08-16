// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/frag_picker/SidechainContactDistCutoff.fwd.hh
/// @author David Kim (dekim@u.washington.edu)

#ifndef INCLUDED_PROTOCOLS_FRAG_PICKER_SIDECHAINCONTACTDISTCUTOFF_FWD_HH
#define INCLUDED_PROTOCOLS_FRAG_PICKER_SIDECHAINCONTACTDISTCUTOFF_FWD_HH

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace frag_picker {

class SidechainContactDistCutoff;
typedef utility::pointer::shared_ptr<SidechainContactDistCutoff> SidechainContactDistCutoffOP;
typedef utility::pointer::shared_ptr<SidechainContactDistCutoff const> SidechainContactDistCutoffCOP;

}  // namespace frag_picker
}  // namespace protocols

#endif  // PROTOCOLS_FRAG_PICKER_SIDECHAIN_CONTACT_DIST_CUTOFF_FWD_HH_
