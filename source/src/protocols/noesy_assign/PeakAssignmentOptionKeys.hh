// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file FragmentSampler.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @details
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_noesy_assign_PeakAssignmentOptionKeys_HH
#define INCLUDED_protocols_noesy_assign_PeakAssignmentOptionKeys_HH

// Utility headers
#include <basic/options/option_macros.hh>

//Auto Headers
OPT_1GRP_KEY( Real, noesy_weights, chemshift )
//OPT_1GRP_KEY( Real, noesy_weights, network_high )

OPT_1GRP_KEY( Real, noesy_weights, Vmin )
OPT_1GRP_KEY( Real, noesy_weights, symmetry )
OPT_1GRP_KEY( Real, noesy_weights, covalent )
//OPT_1GRP_KEY( Real, noesy_weights, decoys )
OPT_1GRP_KEY( Real, noesy_weights, Smax )

OPT_1GRP_KEY( Real, noesy_weights, network_min )
OPT_1GRP_KEY( Real, noesy_weights, network_atom_min )
OPT_1GRP_KEY( Real, noesy_weights, dcalibrate )
OPT_1GRP_KEY( Real, noesy_weights, calibration_target )

OPT_1GRP_KEY( Boolean, noesy, atom_dependent_calibration  )
OPT_1GRP_KEY( Boolean, noesy, ignore_resonancefile_tolerances  )
OPT_1GRP_KEY( Boolean, noesy, ignore_resonancefile_intensities  )

OPT_2GRP_KEY( RealVector, noesy_weights, defaults, Vmin )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, symmetry )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, covalent )
//OPT_2GRP_KEY( RealVector, noesy_weights, defaults, decoys )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, Smax )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, dcut )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, network_min )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, network_atom_min )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, dcalibrate )
OPT_2GRP_KEY( RealVector, noesy_weights, defaults, calibration_target )


OPT_1GRP_KEY( Real, noesy_weights, centroid_padding )
OPT_1GRP_KEY( Boolean, noesy, map_to_cen_atom )

OPT_2GRP_KEY( Integer, noesy, elim, max_assign )
OPT_2GRP_KEY( Real, noesy, elim, vmin )
OPT_2GRP_KEY( Real, noesy, elim, dist_viol )
OPT_2GRP_KEY( Real, noesy, elim, dcut )


OPT_1GRP_KEY( Real, noesy_weights, cst_strength  )
OPT_1GRP_KEY( Integer, noesy_weights, cycle  )
OPT_1GRP_KEY( Real, noesy_weights, min_symmetry_reinforcement )
OPT_1GRP_KEY( Boolean, noesy, no_network  )
OPT_2GRP_KEY( Real, noesy, network, vmin  )
OPT_2GRP_KEY( Real, noesy, network, vmax  )
OPT_2GRP_KEY( Real, noesy, network, reswise_high )
OPT_2GRP_KEY( Real, noesy, network, reswise_min )
OPT_2GRP_KEY( Real, noesy, network, atomwise_min )
OPT_2GRP_KEY( Boolean, noesy, network, include_reverse_dir  )
OPT_2GRP_KEY( Boolean, noesy, network, allow_same_residue_connect  )
OPT_2GRP_KEY( Boolean, noesy, network, use_all_covalent_atoms  )
OPT_2GRP_KEY( String, noesy, network, mode  )

OPT_1GRP_KEY( Boolean, noesy, use_local_distviol )
OPT_2GRP_KEY( Real, noesy, local_distviol, range )
OPT_2GRP_KEY( Real, noesy, local_distviol, global_buffer )
OPT_2GRP_KEY( Real, noesy, local_distviol, global_factor )
OPT_2GRP_KEY( Real, noesy, local_distviol, cutoff )
OPT_2GRP_KEY( Real, noesy, local_distviol, cutoff_buffer )


OPT_2GRP_KEY( Real, noesy, calibration, convergence )
OPT_2GRP_KEY( Real, noesy, calibration, max_noe_dist )
OPT_2GRP_KEY( Real, noesy, calibration, start_nudging )
OPT_2GRP_KEY( Real, noesy, calibration, stop_nudging )
OPT_2GRP_KEY( Real, noesy, calibration, max_nudging )
OPT_2GRP_KEY( Boolean, noesy, calibration, ignore_eliminated_peaks )

OPT_2GRP_KEY( Boolean, noesy, calibration, eliminate )
OPT_2GRP_KEY( Boolean, noesy, calibration, use_median )
OPT_2GRP_KEY( Integer, noesy, calibration, cycles )


OPT_3GRP_KEY( RealVector, noesy, prob, sigmoid, m )
OPT_3GRP_KEY( RealVector, noesy, prob, sigmoid, tau )
OPT_3GRP_KEY( RealVector, noesy, prob, sigmoid, w )
OPT_2GRP_KEY( RealVector, noesy, prob, level )
#endif
