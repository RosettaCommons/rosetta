// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  apps/pilot/membrane/membrane_ddG.cc
///
/// @brief  Membrane Framework Application: ddG of Mutation in Membranes
/// @details The membrane framework currently supports computing ddG of mutations:
///    Uses membrane representation, energ functions, and current tools for Rosetta
///    design for computing ddGs.
///    Last Modified: 7/26/14
///
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/membrane/ddG/MembraneDDGMover.hh>
#include <protocols/membrane/ddG/Mutation.hh>

// Project Headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/internal_util.hh>

#include <core/chemical/AA.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

using namespace protocols::membrane::ddG;

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.ralford.membrane_ddG" );

/// @brief ompLA Task
utility::vector1< MutationOP >
ompLA_task() {

	using namespace protocols::membrane::ddG;

	// Create a new list of mutations
	utility::vector1< MutationOP > mutations;
	mutations.resize(20);

	// Append a list
	mutations[1] = new Mutation( 181, core::chemical::aa_ala );
	mutations[2] = new Mutation( 181, core::chemical::aa_cys );
	mutations[3] = new Mutation( 181, core::chemical::aa_asp );
	mutations[4] = new Mutation( 181, core::chemical::aa_glu );
	mutations[5] = new Mutation( 181, core::chemical::aa_phe );
	mutations[6] = new Mutation( 181, core::chemical::aa_gly );
	mutations[7] = new Mutation( 181, core::chemical::aa_his );
	mutations[8] = new Mutation( 181, core::chemical::aa_ile );
	mutations[9] = new Mutation( 181, core::chemical::aa_lys );
	mutations[10] = new Mutation( 181, core::chemical::aa_leu );
	mutations[11] = new Mutation( 181, core::chemical::aa_met );
	mutations[12] = new Mutation( 181, core::chemical::aa_asn );
	mutations[13] = new Mutation( 181, core::chemical::aa_pro );
	mutations[14] = new Mutation( 181, core::chemical::aa_gln );
	mutations[15] = new Mutation( 181, core::chemical::aa_arg );
	mutations[16] = new Mutation( 181, core::chemical::aa_ser );
	mutations[17] = new Mutation( 181, core::chemical::aa_thr );
	mutations[18] = new Mutation( 181, core::chemical::aa_val );
	mutations[19] = new Mutation( 181, core::chemical::aa_trp );
	mutations[20] = new Mutation( 181, core::chemical::aa_tyr );

	return mutations;
}

/// @brief Lukas Tamm - Aromatic residues at the interface regions task
/// @remarks protein name = ompA
utility::vector1< MutationOP >
ompA_task() {

	using namespace protocols::membrane::ddG;

	// Create a new list of mutations
	utility::vector1< MutationOP > mutations;
	mutations.resize(12);

	// Append a list
	mutations[1] = new Mutation( 7, core::chemical::aa_ala );
	mutations[2] = new Mutation( 15, core::chemical::aa_ala );
	mutations[3] = new Mutation( 57, core::chemical::aa_ala );
	mutations[4] = new Mutation( 143, core::chemical::aa_ala );
	mutations[5] = new Mutation( 43, core::chemical::aa_ala );
	mutations[6] = new Mutation( 55, core::chemical::aa_ala );
	mutations[7] = new Mutation( 129, core::chemical::aa_ala );
	mutations[8] = new Mutation( 141, core::chemical::aa_ala );
	mutations[9] = new Mutation( 168, core::chemical::aa_ala );
	mutations[10] = new Mutation( 51, core::chemical::aa_ala );
	mutations[11] = new Mutation( 123, core::chemical::aa_ala );
	mutations[12] = new Mutation( 170, core::chemical::aa_ala );

	return mutations;
}

/// @brief Bowie - Bacteriorhodopsin B Helix Mutations
/// @remarks PDBID = 1PY6
utility::vector1< MutationOP >
bacteriorhodpsin_task() {

	using namespace protocols::membrane::ddG;

	// Create a new list of mutations
	utility::vector1< MutationOP > mutations;
	mutations.resize(12);

	// Append a list
	mutations[1] = new Mutation( 31, core::chemical::aa_ala );
	mutations[2] = new Mutation( 32, core::chemical::aa_ala );
	mutations[3] = new Mutation( 33, core::chemical::aa_ala );
	mutations[4] = new Mutation( 34, core::chemical::aa_ala );
	mutations[5] = new Mutation( 36, core::chemical::aa_ala );
	mutations[6] = new Mutation( 37, core::chemical::aa_ala );
	mutations[7] = new Mutation( 38, core::chemical::aa_ala );
	mutations[8] = new Mutation( 39, core::chemical::aa_ala );
	mutations[9] = new Mutation( 41, core::chemical::aa_ala );
	mutations[10] = new Mutation( 42, core::chemical::aa_ala );
	mutations[11] = new Mutation( 45, core::chemical::aa_ala );
	mutations[12] = new Mutation( 44, core::chemical::aa_ala );
	mutations[10] = new Mutation( 45, core::chemical::aa_ala );
	mutations[11] = new Mutation( 46, core::chemical::aa_ala );
	mutations[12] = new Mutation( 48, core::chemical::aa_ala );
	mutations[11] = new Mutation( 50, core::chemical::aa_ala );
	mutations[12] = new Mutation( 51, core::chemical::aa_ala );
	mutations[10] = new Mutation( 52, core::chemical::aa_ala );
	mutations[11] = new Mutation( 53, core::chemical::aa_ala );
	mutations[12] = new Mutation( 54, core::chemical::aa_ala );
	mutations[12] = new Mutation( 55, core::chemical::aa_ala );
	mutations[10] = new Mutation( 56, core::chemical::aa_ala );
	mutations[11] = new Mutation( 57, core::chemical::aa_ala );
	mutations[12] = new Mutation( 58, core::chemical::aa_ala );

	return mutations;
}

/// @brief Main method
int
main( int argc, char * argv [] )
{

	using namespace core::scoring;
	using namespace protocols::membrane::ddG;
	using namespace protocols::jd2;

	try {

		// Devel init factories
		devel::init(argc, argv);

		// Register JD2 options
		protocols::jd2::register_options();

		// Load in mutants for the appropriate task
		utility::vector1< MutationOP > ompA_mutations = ompA_task();
		utility::vector1< MutationOP > ompLA_mutations = ompLA_task();
		utility::vector1< MutationOP > bacteriorhodopsin_mutations = bacteriorhodpsin_task();

		// Make ddG calculation
		MembraneDDGMoverOP ddG = new MembraneDDGMover( ompLA_mutations );
		JobDistributor::get_instance()->go( ddG );

		return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}


