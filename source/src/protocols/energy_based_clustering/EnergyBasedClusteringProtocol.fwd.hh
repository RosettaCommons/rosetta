// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cluster/energy_based_clustering/EnergyBasedClusteringProtocol.fwd.hh
/// @brief Performs the work done by the energy_based_clustering app.  Uses an energy-biased cookie-cutter approach to
/// cluster a large number of structures without generating an all-by-all RMSD matrix.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_protocols_cluster_energy_based_clustering_EnergyBasedClusteringProtocol_fwd_hh
#define INCLUDED_protocols_cluster_energy_based_clustering_EnergyBasedClusteringProtocol_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

// Forward
namespace protocols {
namespace energy_based_clustering {

class EnergyBasedClusteringProtocol;

typedef utility::pointer::shared_ptr< EnergyBasedClusteringProtocol > EnergyBasedClusteringProtocolOP;
typedef utility::pointer::shared_ptr< EnergyBasedClusteringProtocol const > EnergyBasedClusteringProtocolCOP;

} //energy_based_clustering
} //protocols

#endif //INCLUDED_protocols_cluster_energy_based_clustering_EnergyBasedClusteringProtocol_fwd_hh
