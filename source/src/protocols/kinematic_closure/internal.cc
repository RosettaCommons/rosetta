// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief  Source file for internal helpers.
/// @author Kale Kundert (kale.kundert@ucsf.edu)

#include <protocols/kinematic_closure/internal.hh>

namespace protocols {
namespace kinematic_closure {

const Real IdealParameters::c_n_ca_angle = 121.7;
const Real IdealParameters::n_ca_c_angle = 111.2;
const Real IdealParameters::ca_c_n_angle = 116.2;
const Real IdealParameters::c_n_length = 1.32869;
const Real IdealParameters::n_ca_length = 1.458;
const Real IdealParameters::ca_c_length = 1.52326;
const Real IdealParameters::omega_dihedral = 179.8;
const Real IdealParameters::mean_n_ca_c_angle = 110.86;
const Real IdealParameters::std_n_ca_c_angle = 2.48;
const Real IdealParameters::min_n_ca_c_angle = 105.90;
const Real IdealParameters::max_n_ca_c_angle = 118.94;

Size num_rama_filter_fails = 0;
Size num_bump_filter_fails = 0;

} // end namespace kinematic_closure
} // end namespace protocols




