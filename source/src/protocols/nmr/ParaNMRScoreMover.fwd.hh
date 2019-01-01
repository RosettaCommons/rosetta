// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/ParaNMRScoreMover.fwd.hh
/// @brief   forward declaration for ParaNMRScoreMover
/// @details last Modified: 11/15/18
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_ParaNMRScoreMover_FWD_HH
#define INCLUDED_protocols_nmr_ParaNMRScoreMover_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace nmr {

class ParaNMRScoreMover;

typedef utility::pointer::shared_ptr< ParaNMRScoreMover > ParaNMRScoreMoverOP;
typedef utility::pointer::shared_ptr< ParaNMRScoreMover const > ParaNMRScoreMoverCOP;
typedef utility::pointer::weak_ptr< ParaNMRScoreMover > ParaNMRScoreMoverAP;
typedef utility::pointer::weak_ptr< ParaNMRScoreMover const > ParaNMRScoreMoverCAP;

} // namespace nmr
} // namespace protocols

#endif // INCLUDED_protocols_nmr_ParaNMRScoreMover_FWD_HH
