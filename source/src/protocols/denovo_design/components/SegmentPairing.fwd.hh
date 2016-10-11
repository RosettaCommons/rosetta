// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/SegmentPairing.fwd.hh
/// @brief Handles user-specified pairing between/among segments
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_components_SegmentPairing_fwd_hh
#define INCLUDED_protocols_denovo_design_components_SegmentPairing_fwd_hh

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.fwd.hh>

// Forward
namespace protocols {
namespace denovo_design {
namespace components {

/// @brief Individual strands are oriented pointing either "UP" or "DOWN"
///        If two adjacent strands have the same orientation, they are parallel
///        If two adjacent strands have different orientation, they are antiparallel
enum StrandOrientation {
	UP = 1,
	DOWN = 2,
	ORIENTATIONS_END = 3
};
typedef utility::vector1< StrandOrientation > StrandOrientations;

typedef long int RegisterShift;
typedef utility::vector1< RegisterShift > RegisterShifts;

typedef std::pair< core::Size, core::Size > ResiduePair;
typedef utility::vector1< ResiduePair > ResiduePairs;

class SegmentPairing;
typedef utility::pointer::shared_ptr< SegmentPairing > SegmentPairingOP;
typedef utility::pointer::shared_ptr< SegmentPairing const > SegmentPairingCOP;
typedef utility::vector1< SegmentPairingOP > SegmentPairingOPs;
typedef utility::vector1< SegmentPairingCOP > SegmentPairingCOPs;

class HelixPairing;
typedef utility::pointer::shared_ptr< HelixPairing > HelixPairingOP;
typedef utility::pointer::shared_ptr< HelixPairing const > HelixPairingCOP;
typedef utility::vector1< HelixPairingOP > HelixPairingOPs;
typedef utility::vector1< HelixPairingCOP > HelixPairingCOPs;

class StrandPairing;
typedef utility::pointer::shared_ptr< StrandPairing > StrandPairingOP;
typedef utility::pointer::shared_ptr< StrandPairing const > StrandPairingCOP;
typedef utility::vector1< StrandPairingOP > StrandPairingOPs;
typedef utility::vector1< StrandPairingCOP > StrandPairingCOPs;

class HelixSheetPairing;
typedef utility::pointer::shared_ptr< HelixSheetPairing > HelixSheetPairingOP;
typedef utility::pointer::shared_ptr< HelixSheetPairing const > HelixSheetPairingCOP;
typedef utility::vector1< HelixSheetPairingOP > HelixSheetPairingOPs;
typedef utility::vector1< HelixSheetPairingCOP > HelixSheetPairingCOPs;


} //protocols
} //denovo_design
} //components

#endif //INCLUDED_protocols_denovo_design_components_SegmentPairing_fwd_hh

