// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/rms/rms_enum.hh
/// @brief  RMS Enums using in SimpleMetrics for RMSD
/// @uathor Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_scoring_rms_enum_HH
#define INCLUDED_core_scoring_rms_enum_HH



namespace core {
namespace scoring {

enum rmsd_atoms{
	rmsd_protein_bb_heavy = 1,
	rmsd_protein_bb_heavy_including_O,
	rmsd_protein_bb_ca,
	rmsd_sc_heavy,
	rmsd_sc,
	rmsd_all_heavy,
	rmsd_all,
	rmsd_atom_total = rmsd_all
};

} // end namespace scoring
} // end namespace core

#endif
