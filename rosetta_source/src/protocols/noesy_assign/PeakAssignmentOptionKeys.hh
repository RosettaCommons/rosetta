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
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


#ifndef INCLUDED_protocols_noesy_assign_PeakAssignmentOptionKeys_HH
#define INCLUDED_protocols_noesy_assign_PeakAssignmentOptionKeys_HH

// Utility headers
#include <basic/options/option_macros.hh>

#ifdef DEFINE_OPTIONS_NOW
#define EXPORT( optmacro ) optmacro
#else
#define EXPORT( optmacro ) EXTERN_##optmacro
#endif


//Auto Headers
EXPORT( OPT_1GRP_KEY( Real, noesy_weights, chemshift  ) )
//EXPORT( OPT_1GRP_KEY( Real, noesy_weights, network_high  ) )

EXPORT( OPT_1GRP_KEY( Real, noesy_weights, Vmin  ) )
EXPORT( OPT_1GRP_KEY( Real, noesy_weights, symmetry  ) )
EXPORT( OPT_1GRP_KEY( Real, noesy_weights, covalent  ) )
EXPORT( OPT_1GRP_KEY( Real, noesy_weights, decoys  ) )
EXPORT( OPT_1GRP_KEY( Real, noesy_weights, Smax  ) )
EXPORT( OPT_1GRP_KEY( Real, noesy_weights, dcut  ) )
EXPORT( OPT_1GRP_KEY( Real, noesy_weights, network_min  ) )
EXPORT( OPT_1GRP_KEY( Real, noesy_weights, network_atom_min  ) )
EXPORT( OPT_1GRP_KEY( Real, noesy_weights, dcalibrate  ) )
EXPORT( OPT_1GRP_KEY( Real, noesy_weights, calibration_target  ) )

EXPORT( OPT_1GRP_KEY( Boolean, noesy, atom_dependent_calibration  ) )
EXPORT( OPT_1GRP_KEY( Boolean, noesy, ignore_resonancefile_tolerances  ) )

EXPORT( OPT_2GRP_KEY( RealVector, noesy_weights, defaults, Vmin  ) )
EXPORT( OPT_2GRP_KEY( RealVector, noesy_weights, defaults, symmetry  ) )
EXPORT( OPT_2GRP_KEY( RealVector, noesy_weights, defaults, covalent  ) )
EXPORT( OPT_2GRP_KEY( RealVector, noesy_weights, defaults, decoys  ) )
EXPORT( OPT_2GRP_KEY( RealVector, noesy_weights, defaults, Smax  ) )
EXPORT( OPT_2GRP_KEY( RealVector, noesy_weights, defaults, dcut  ) )
EXPORT( OPT_2GRP_KEY( RealVector, noesy_weights, defaults, network_min  ) )
EXPORT( OPT_2GRP_KEY( RealVector, noesy_weights, defaults, network_atom_min  ) )
EXPORT( OPT_2GRP_KEY( RealVector, noesy_weights, defaults, dcalibrate  ) )
EXPORT( OPT_2GRP_KEY( RealVector, noesy_weights, defaults, calibration_target  ) )


EXPORT( OPT_1GRP_KEY( Real, noesy_weights, elim_dist_viol  ) )
EXPORT( OPT_1GRP_KEY( Real, noesy_weights, centroid_padding  ) )
EXPORT( OPT_1GRP_KEY( Boolean, noesy, map_to_cen_atom  ) )


EXPORT( OPT_1GRP_KEY( Real, noesy_weights, cst_strength  ) )
EXPORT( OPT_1GRP_KEY( Integer, noesy_weights, cycle  ) )

EXPORT( OPT_1GRP_KEY( Boolean, noesy, no_network  ) )
EXPORT( OPT_2GRP_KEY( Boolean, noesy, network, include_reverse_dir  ) )
EXPORT( OPT_2GRP_KEY( Boolean, noesy, network, allow_same_residue_connect  ) )
EXPORT( OPT_2GRP_KEY( Boolean, noesy, network, use_all_covalent_atoms  ) )

EXPORT( OPT_1GRP_KEY( Boolean, noesy, use_local_distviol  ) )
EXPORT( OPT_2GRP_KEY( Real, noesy, local_distviol, range  ) )
EXPORT( OPT_2GRP_KEY( Real, noesy, local_distviol, global_buffer  ) )

#endif
