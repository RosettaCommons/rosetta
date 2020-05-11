// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/epr_deer/DEERIO.fwd.hh
/// @brief  IO class for data obtained with double electron-electron resonance (DEER)
/// @details This class is only called once at the beginning of the application by the
///      energy method. It then parses the input file and the electron coordinate file.
///      The input file is then used to create a "map" (essentially a fake graph) which
///      then becomes the basis of the DEERData object.

/// @author  Diego del Alamo ( del.alamo@vanderbilt.edu )

#ifndef INCLUDED_core_scoring_epr_deer_DEERIO_fwd_hh
#define INCLUDED_core_scoring_epr_deer_DEERIO_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace epr_deer {

class DEERDistanceIO;

typedef utility::pointer::shared_ptr< DEERDistanceIO > DEERIOOP;
typedef utility::pointer::shared_ptr< DEERDistanceIO const > DEERIOCOP;

} // namespace epr_deer
} // namespace scoring
} // namespace core


#endif // INCLUDED_core_scoring_epr_deer_distance_DEERDistanceIO_fwd_hh
