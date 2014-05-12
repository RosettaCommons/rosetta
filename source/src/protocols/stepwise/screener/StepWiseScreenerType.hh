// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/StepWiseScreenerType.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_StepWiseScreenerType_HH
#define INCLUDED_protocols_stepwise_StepWiseScreenerType_HH

namespace protocols {
namespace stepwise {
namespace screener {

	enum StepWiseScreenerType{
		NONE,
		ANCHOR_SUGAR,
    ATR_REP,
		BASE_BIN_MAP,
		BASE_CENTROID,
		BULGE_APPLIER,
		CENTROID_DISTANCE,
		FAST_FORWARD_TO_NEXT_RIGID_BODY,
		FAST_FORWARD_TO_NEXT_RESIDUE_ALTERNATIVE,
		INTEGRATION_TEST,
		NATIVE_RMSD,
		O2PRIME_PACK,
		PACK,
		PARTITION_CONTACT,
		PHOSPHATE_PACK,
		POSE_SELECTION,
		PROTEIN_ATR_REP,
		RESIDUE_APPLIER,
		RESIDUE_CONTACT,
		PROTEIN_CCD_CLOSURE,
		RNA_CHAIN_CLOSURE,
		RNA_CHAIN_CLOSABLE_GEOMETRY,
		RNA_CHAIN_CLOSABLE_GEOMETRY_RESIDUE_BASED,
		RNA_CHAIN_CLOSABLE_GEOMETRY_STUB_BASED,
		SAMPLE_APPLIER,
		SAMPLE_APPLIER_WITH_RESIDUE_LIST,
		SCORER,
		SIMPLE_POSE_SELECTION,
		SIMPLE_RMSD,
		STUB_APPLIER,
		STUB_DISTANCE,
		SUGAR_INSTANTIATOR,
		TAG_DEFINITION,
		VDW_ATR_REP,
		VDW_BIN
	};

} //screener
} //stepwise
} //protocols

#endif
