// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/public/clustering/energy_based_clustering.cc
/// @brief Cluster a set of structures using an energy-biased cluster centre selection algorithm.
/// @author Vikram K. Mulligan (vmullig@uw.edu)
/// @details This file started off as the bettercluster.cc pilot app in src/apps/pilot/vmullig/.
///
///     This uses a "cookie-cutter" based approach.  The basic algorithm is as follows:
///         1.  Score all input structures, and store minimal information needed for clustering for each structure in an unclustered pool.
///         2.  Select the lowest-energy structure in the pool, and transfer it from the pool into a new cluster.  This is now the centre of the new
///         cluster.  Calculate the RMSD of every structure remaining in the unclustered pool to the new cluster centre.  (If necessary, consider
///         circular permutations as well).  Transfer any structure with an RMSD below the threshold from the unclustered pool to the new cluster.
///         3.  Repeat step 2 with the structures remaining in the unclustered pool, generating new clusters until no structures remain.
///
///     This algorithm has the advantage of exactness and non-stochasticity without requiring a full RMSD matrix of all structures to all structures.
///     It works well for very large datasets.
///
/// Original description: A better clustering algorithm for generating conformational libraries for
/// explicit multistate design.
///
/// History:
/// --File created 6 May 2013 by Vikram K. Mulligan, Baker Laboratory.
/// --Modified 3 Jun 2013 to make Cartesian-based clustering much faster.  (No
/// more rebuilding poses).
/// --Modified 5 Jun 2013 to allow clustering of backbone-cyclized peptides.
/// --Modified 28 Aug 2013 to allow N-offset cyclic permutations to be clustered.  (e.g. Offsets
///  are incremented by 2 residues, or by 3 residues, or whatever.)
/// --Modified 4 Sept 2013 to allow clustering of peptides with beta-amino acid residues:
///   --check for beta residues based on first structure loaded (they must be consistent -- DONE).
///   --align CM atoms, if present (DONE).
///   --update information that is stored for cartesian or dihedral clustering (DONE).
///   --update PCA file output --> list of backbone dihedrals must include theta (DONE, though
///   other apps must also be updated to READ these properly).
/// --Modified 9 Oct 2013 to allow clustering of homooligomers (with swapping of oligomer subunits
///  during RMSD calculation).
///    TODO:
///    --store jumps for PCA analysis
///    --calculate RMSD for all possible permutations of homooligomer subunits and keep lowest
/// --Modified 21 May 2014 to allow silent file output.
/// --Modified 22 May 2014:
///  --Added support for constraints files.
///  --Added support for a user-defined list of additional atoms to use in RMSD
///  calculation.
/// --Modified 18 June 2014:
///  --Took out all references to the Alglib library, because of the asenine Rosettacommons rules
///  about third-party libraries.
///  --Added a PCA function to the numeric library, and linked this to that.
/// --Modified 11 Aug 2014:
///  --Added an option to ignore entire chains in the RMSD calculation.
///  --Added an option to skip PCA analysis.
/// --Modified 21 Sept 2015:
///     --Added an option to dump out only the first N clusters.
/// --Modified 30 Jan 2017:
///     --Added an option to filter out structures that aren't symmetric.
///     --Modified 18 July 2017:
///         --Converted to a true public application, added tests, and cleaned up implementation.  Temporarily disabled PCA analysis.

// Rosetta headers:
#include <devel/init.hh>
#include <protocols/energy_based_clustering/EnergyBasedClusteringProtocol.hh>
#include <basic/Tracer.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers:
#include <stdio.h>

static basic::Tracer TR( "apps.public.clustering.energy_based_clustering" );

/// @brief Minimial main() function for energy_based_clustering application.  Most of the heavy lifting is done by
/// the protocols::energy_based_clustering::EnergyBasedClusteringProtocol class.
int main( int argc, char * argv [] ) {
	using namespace protocols::energy_based_clustering;
	try {

		EnergyBasedClusteringProtocol::register_options();
		devel::init(argc, argv);

		TR << "Starting energy_based_clustering application." << std::endl;
		TR << "This application was first called \"bettercluster\", and was created on 6 May 2013 by Vikram K. Mulligan (vmullig@uw.edu)." << std::endl;
		TR << "It was upgraded to a public application on 18 July 2017." << std::endl;
		TR << "\n************************************************************************************\nIF YOU USE THIS APPLICATION, PLEASE CITE:\nHosseinzadeh P., Bhardwaj G., Mulligan V., et al.  (2017).  Manuscript under review.\n************************************************************************************" << std::endl;

		EnergyBasedClusteringOptions options(true);
		EnergyBasedClusteringProtocol cluster_protocol( options );
		cluster_protocol.go();

		TR << "Terminating energy_based_clustering application with exit code 0 (no errors)." << std::endl;
	} catch (utility::excn::Exception& excn ) {
		std::cerr << "Terminating energy_based_clustering application with errors:" << std::endl;
		excn.show( std::cerr );
		return -1;
	}

	return 0;
}

