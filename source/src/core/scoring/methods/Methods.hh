// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/Methods.hh
/// @brief  Energy Method enumeration
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_methods_Methods_hh
#define INCLUDED_core_scoring_methods_Methods_hh

namespace core {
namespace scoring {
namespace methods {


enum EnergyMethods
{
	etable_method = 1, // first method begins at one for one-based indexing
	dunbrack_method,
	hbond_method,
	elec_method,
	lkball_method,
	mm_lj_energy_inter_method,
	pair_e_method, // give this a new name!
	reference_e_method,
	vdw_method,
	ramachandran_method,
	n_energy_methods = ramachandran_method // keep this guy last
};


// Eventually, gen born and classic coulombic energies will go here
enum LongRangeEnergyType {
	constraints_lr = 1,
	gen_born_lr,
	multipole_elec_lr,
	PB_elec_lr,
	cart_bonded_lr,
	rama2b_lr,
	rna_suite_lr,
	DFIRE,
	sym_bonus_lr,
	elec_dens_energy,
	elec_dens_fast_energy,
	elec_dens_cen_energy,
	elec_dens_allatom_cen_energy,
	elec_dens_atomwise_energy,
	patterson_corr_energy,
	fa_disulfide_energy,
	disulfide_matching_energy,
	centroid_disulfide_energy,
	n_long_range_types = centroid_disulfide_energy // keep this guy last
};

} // namespace methods
} // namespace scoring
} // namespace core

#endif
