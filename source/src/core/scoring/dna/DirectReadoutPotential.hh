// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/DirectReadoutPotential.hh
/// @brief  1st pass implementation of Kono + Sarai's protein-DNA interaction potential
/// @details  Needs polishing, converting to mini standards in some respects, but still in trial stage.
/// @author Amy Ticoll


#ifndef INCLUDED_core_scoring_dna_DirectReadoutPotential_hh
#define INCLUDED_core_scoring_dna_DirectReadoutPotential_hh

// Unit Headers
#include <core/scoring/dna/DirectReadoutPotential.fwd.hh>

// Package headers

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <string>

namespace core {
namespace scoring {
namespace dna {

/// @brief  1st pass implementation of Kono + Sarai's protein-DNA interaction potential
/// @details  Needs polishing, converting to mini standards in some respects, but still in trial stage.

class DirectReadoutPotential : public utility::pointer::ReferenceCount {

public:
	typedef std::string string;

public:

	/// @brief ctor, reads data file
	DirectReadoutPotential();

	// get score method
	Real rsd_rsd_energy(conformation::Residue const & rsd1,conformation::Residue const & rsd2) const;

private:
	Real score[9][9][4][20][4];
	int num_pairs[4][20];
	//  int pair_base[80];
	//  int pair_aa[80];
	//  int pair_num[80];

	string A_bins[9][9][4];
	string G_bins[9][9][4];
	string T_bins[9][9][4];
	string C_bins[9][9][4];

	int aas_at_grid;
	const Real wt;
	const Real RT;

	void fill_bins(string (&my_array)[9][9][4], char const base);
	void get_pairs();
	int get_xy_bin(Real coord) const;
	int get_z_bin(Real coord) const;

};

}
}
}

#endif
