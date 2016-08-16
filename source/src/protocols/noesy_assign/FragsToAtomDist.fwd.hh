// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragments/FragstoAtomDist.fwd.hh
/// @brief  simulate fragments and recored atom distances
/// @author Zaiyong Zhang (zaiyong.zhang@tum.de)
/// @date   Thr NOV 22 13:22:31 2012


#ifndef INCLUDED_protocols_noesy_assign_FragstoATomDist_fwd_hh
#define INCLUDED_protocols_noesy_assign_FragstoATomDist_fwd_hh

// Utility headers
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace noesy_assign {

// Forward
class FragsToAtomDist;

// Types
typedef  utility::pointer::shared_ptr< FragsToAtomDist >  FragsToAtomDistOP;
typedef  utility::pointer::shared_ptr< FragsToAtomDist const >  FragsToAtomDistCOP;

typedef  utility::pointer::weak_ptr< FragsToAtomDist >  FragsToAtomDistAP;
typedef  utility::pointer::weak_ptr< FragsToAtomDist const >  FragsToAtomDistCAP;

} // namespace noesy_assign
} // namespace protocols

#endif
