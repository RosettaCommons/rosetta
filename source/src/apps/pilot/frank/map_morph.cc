// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#include <protocols/viewer/viewers.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>

#include <basic/options/util.hh>
#include <devel/init.hh>

#include <numeric/xyz.functions.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>

#include <ObjexxFCL/string.functions.hh>
#include <core/kinematics/Jump.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <ObjexxFCL/FArray3D.hh>

// C++ headers
#include <iostream>
#include <fstream>
#include <string>


using namespace core;
using namespace basic;
using namespace utility;
using namespace protocols;


void
map_morph()
{
	using namespace protocols::moves;
	using namespace scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	core::Size axis=3; // z

	core::scoring::electron_density::ElectronDensity &edm = core::scoring::electron_density::getDensityMap();

	// get coords of helical axis
	numeric::xyzVector<int> grid = edm.getGrid();
	numeric::xyzVector<core::Real> midpt(1.0+grid[0]/2.0, 1.0+grid[1]/2.0, 1.0+grid[2]/2.0);

	//core::Real theta = 1.3745;
	//core::Real rise = 25.22;

	ObjexxFCL::FArray3D< double > spline;
	core::scoring::electron_density::spline_coeffs( edm.data(), spline );

	// now perturb
	core::scoring::electron_density::ElectronDensity edmNEW = edm;
	for (int i_th=-2; i_th<=2; ++i_th) {
		ObjexxFCL::FArray3D< float > newdata = edm.data();

		core::Real del_theta_per_A_rise = 0.005*i_th;


		for (int z=1; z<=grid[2]; ++z) {
			//std::cerr << z << "/" << grid[2] << std::endl;
			for (int y=1; y<=grid[1]; ++y) {
				for (int x=1; x<=grid[0]; ++x) {
					numeric::xyzVector<core::Real> idx( x, y, z );
					idx = idx-midpt;
					numeric::xyzMatrix<core::Real> R = numeric::rotation_matrix_radians( numeric::xyzVector<core::Real>(0,0,1), del_theta_per_A_rise*idx[2] );
					numeric::xyzVector<core::Real> idxNEW = R*(idx) + midpt;
					core::Real newVal = core::scoring::electron_density::interp_spline( spline, idxNEW );
					newdata( x,y,z ) = newVal;
				}
			}
		}

		edmNEW.set( newdata);

		// dump
		edmNEW.writeMRC( (del_theta_per_A_rise>=0 ? "theta_p" : "theta_m")+utility::to_string( std::abs(del_theta_per_A_rise) )+".mrc" );
	}
}


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		devel::init( argc, argv );
		map_morph();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
