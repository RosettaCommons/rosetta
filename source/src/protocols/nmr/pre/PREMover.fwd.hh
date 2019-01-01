// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pre/PREMover.fwd.hh
/// @brief   forward declaration for PREMover
/// @details last Modified: 05/30/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_pre_PREMover_FWD_HH
#define INCLUDED_protocols_nmr_pre_PREMover_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace protocols {
namespace nmr {
namespace pre {

class PREDistanceRecord;
class PREMover;

typedef utility::pointer::shared_ptr< PREDistanceRecord > PREDistanceRecordOP;
typedef utility::pointer::shared_ptr< PREDistanceRecord const > PREDistanceRecordCOP;
typedef utility::pointer::weak_ptr< PREDistanceRecord > PREDistanceRecordAP;
typedef utility::pointer::weak_ptr< PREDistanceRecord const > PREDistanceRecordCAP;

typedef utility::pointer::shared_ptr< PREMover > PREMoverOP;
typedef utility::pointer::shared_ptr< PREMover const > PREMoverCOP;
typedef utility::pointer::weak_ptr< PREMover > PREMoverAP;
typedef utility::pointer::weak_ptr< PREMover const > PREMoverCAP;


} // namespace pre
} // namespace nmr
} // namespace protocols

#endif // INCLUDED_protocols_nmr_pre_PREMover_FWD_HH
