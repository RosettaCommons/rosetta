// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/NMRDummySpinlabelVoxelGrid.fwd.hh
/// @brief   forward declaration for NMRDummySpinlabelVoxelGrid.hh
/// @details last Modified: 02/11/17
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_NMRDummySpinlabelVoxelGrid_FWD_HH
#define INCLUDED_core_scoring_nmr_NMRDummySpinlabelVoxelGrid_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/access_ptr.hh>

namespace core {
namespace scoring {
namespace nmr {

class VoxelGridPoint;
class NMRDummySpinlabelAtom;
class VoxelGridPoint_AA;
class NMRDummySpinlabelVoxelGrid;

typedef utility::pointer::shared_ptr< VoxelGridPoint > VoxelGridPointOP;
typedef utility::pointer::shared_ptr< VoxelGridPoint const > VoxelGridPointCOP;
typedef utility::pointer::weak_ptr< VoxelGridPoint > VoxelGridPointAP;
typedef utility::pointer::weak_ptr< VoxelGridPoint const > VoxelGridPointCAP;

typedef utility::pointer::shared_ptr< NMRDummySpinlabelAtom > NMRDummySpinlabelAtomOP;
typedef utility::pointer::shared_ptr< NMRDummySpinlabelAtom const > NMRDummySpinlabelAtomCOP;
typedef utility::pointer::weak_ptr< NMRDummySpinlabelAtom > NMRDummySpinlabelAtomAP;
typedef utility::pointer::weak_ptr< NMRDummySpinlabelAtom const > NMRDummySpinlabelAtomCAP;

typedef utility::pointer::shared_ptr< VoxelGridPoint_AA > VoxelGridPoint_AAOP;
typedef utility::pointer::shared_ptr< VoxelGridPoint_AA const > VoxelGridPoint_AACOP;
typedef utility::pointer::weak_ptr< VoxelGridPoint_AA > VoxelGridPoint_AAAP;
typedef utility::pointer::weak_ptr< VoxelGridPoint_AA const > VoxelGridPoint_AACAP;

typedef utility::pointer::shared_ptr< NMRDummySpinlabelVoxelGrid > NMRDummySpinlabelVoxelGridOP;
typedef utility::pointer::shared_ptr< NMRDummySpinlabelVoxelGrid const > NMRDummySpinlabelVoxelGridCOP;
typedef utility::pointer::weak_ptr< NMRDummySpinlabelVoxelGrid > NMRDummySpinlabelVoxelGridAP;
typedef utility::pointer::weak_ptr< NMRDummySpinlabelVoxelGrid const > NMRDummySpinlabelVoxelGridCAP;


} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_NMRDummySpinlabelVoxelGrid_FWD_HH
