// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   test_sugar_torsion_getters_and_setters.cc
/// @brief  This is simple pilot app for testing carbohydrates.
/// @note   I intend to convert this into unit tests once everything is sound.
/// @author Labonte

// includes
#include <iostream>

#include <devel/init.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>


int
main( int argc, char *argv[] )
{
	using namespace std;

	try {
		using namespace core;
		using namespace import_pose;
		using namespace pose;
		using namespace conformation;

		// Initialize core.
		devel::init( argc, argv );

		// Import test pose.
		Pose pose;
		pose_from_pdb( pose, "../test/core/chemical/carbohydrates/Lex.pdb" );

		cout << "Lewisx:" << endl << pose << endl << endl;


		// Extract residues for ease of use.
		ResidueCOP res1 = pose.residue(1).get_self_ptr();
		ResidueCOP res2 = pose.residue(2).get_self_ptr();
		ResidueCOP res3 = pose.residue(3).get_self_ptr();

		cout << "Getters:" << endl;
		cout << " Residue 1: " << res1->name() << endl;
		cout << "  Phi: " << pose.phi(1) << " -- expected value: undefined" << endl;
		cout << "  Psi: " << pose.psi(1) << " -- expected value: undefined" << endl;
		cout << "  Omega: " << pose.omega(1) << " -- expected value: undefined" << endl << endl;

		cout << " Residue 2: " << res2->name() << endl;
		cout << "  Phi: " << pose.phi(2) << " -- expected value: -85.8°" << endl;
		cout << "  Psi: " << pose.psi(2) << " -- expected value: 135.6°" << endl;
		cout << "  Omega: " << pose.omega(2) << " -- expected value: undefined" << endl << endl;

		cout << " Residue 3: " << res3->name() << endl;
		cout << "  Phi: " << pose.phi(3) << " -- expected value: -76.9°" << endl;
		cout << "  Psi: " << pose.psi(3) << " -- expected value: -97.0°" << endl;
		cout << "  Omega: " << pose.omega(3) << " -- expected value: undefined" << endl << endl;

		cout << "Setters:" << endl;
		cout << " Residue 1: " << res1->name() << endl;
		cout << "  Setting phi to 180. ";
		pose.set_phi(1, 180.0);
		cout << "Phi is now " << pose.phi(1) << endl;
		cout << "  Setting psi to 180. ";
		pose.set_psi(1, 180.0);
		cout << "Psi is now " << pose.psi(1) << endl;
		cout << "  Setting omega to 180. ";
		pose.set_omega(1, 180.0);
		cout << "Omega is now " << pose.omega(1) << endl << endl;

		cout << " Residue 2: " << res2->name() << endl;
		cout << "  Setting phi to 180. ";
		pose.set_phi(2, 180.0);
		cout << "Phi is now " << pose.phi(2) << endl;
		cout << "  Setting psi to 180. ";
		pose.set_psi(2, 180.0);
		cout << "Psi is now " << pose.psi(2) << endl;
		cout << "  Setting omega to 180. ";
		pose.set_omega(2, 180.0);
		cout << "Omega is now " << pose.omega(2) << endl << endl;

		cout << " Residue 3: " << res3->name() << endl;
		cout << "  Setting phi to 180. ";
		pose.set_phi(3, 180.0);
		cout << "Phi is now " << pose.phi(3) << endl;
		cout << "  Setting psi to 180. ";
		pose.set_psi(3, 180.0);
		cout << "Psi is now " << pose.psi(3) << endl;
		cout << "  Setting omega to 180. ";
		pose.set_omega(3, 180.0);
		cout << "Omega is now " << pose.omega(3) << endl << endl;
	} catch ( utility::excn::EXCN_Base const & e ) {
		cerr << "caught exception " << e.msg() << endl;
		return -1;
	}
	return 0;
}
