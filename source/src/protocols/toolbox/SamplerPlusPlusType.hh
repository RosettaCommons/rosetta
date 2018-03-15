// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/toolbox/SamplerPlusPlusType.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_sampler_plus_plus_SamplerPlusPlusType_HH
#define INCLUDED_sampler_plus_plus_SamplerPlusPlusType_HH

namespace protocols {
namespace toolbox {

enum SamplerPlusPlusType {
	NO_ROTAMER_TYPE,
	ANY,
	COMB,
	COPY_DOF,
	INPUT_STREAM,
	JUMP,
	ONE_TORSION,
	ONE_VALUE,
	ONE_VALUE_COMB,
	PROTEIN_FRAGMENT,
	PROTEIN_MAIN_CHAIN,
	PROTEIN_BETA_ANTIPARALLEL,
	RESIDUE_LIST,
	RESIDUE_ALTERNATIVE,
	RESIDUE_ALTERNATIVE_COMB,
	RING_CONFORMERS,
	RIGID_BODY,
	RIGID_BODY_WITH_RESIDUE_LIST,
	RIGID_BODY_WITH_RESIDUE_ALTERNATIVES,
	RNA_CHI,
	RNA_KIC,
	RNA_KINEMATIC_CLOSER,
	MC_ANY,
	MC_COMB,
	MC_ONE_TORSION,
	MC_LOOP,
	MC_RNA_KIC,
	MC_RNA_MULTI_SUITE,
	MC_RNA_ONE_JUMP,
	MC_RNA_SUGAR,
	MC_RNA_SUITE,
	RNA_NUCLEOSIDE,
	RNA_SUGAR,
	RNA_SUITE,
	SIZED,
	SIZED_ANY,
	SIZED_COMB
};

} //toolbox
} //protocols

#endif
