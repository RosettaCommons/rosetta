// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerType.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_RotamerType_HH
#define INCLUDED_protocols_rotamer_sampler_RotamerType_HH

namespace protocols {
namespace rotamer_sampler {

	enum RotamerType {
		NONE,
		ANY,
		COPY_DOF,
		MC_ONE_TORSION,
		ONE_TORSION,
		ONE_VALUE,
		ONE_VALUE_COMB,
		RESIDUE_LIST,
		RESIDUE_ALTERNATIVE,
		RESIDUE_ALTERNATIVE_COMB,
		RIGID_BODY,
		RIGID_BODY_WITH_RESIDUE_LIST,
		RIGID_BODY_WITH_RESIDUE_ALTERNATIVES,
		RNA_CHI,
		RNA_KIC,
		RNA_KINEMATIC_CLOSER,
		RNA_MC_MULTI_SUITE,
		RNA_MC_SUGAR,
		RNA_MC_SUITE,
		RNA_NUCLEOSIDE,
		RNA_SUGAR,
		RNA_SUITE,
		SIZED,
		SIZED_ANY,
		SIZED_COMB
	};

} //rotamer_sampler
} //protocols

#endif
