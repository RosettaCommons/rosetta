// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

#include <protocols/viewer/viewers.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>

#include <basic/options/util.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>
#include <core/kinematics/Jump.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ headers
#include <iostream>
#include <string>


using namespace core;
using namespace basic;
using namespace utility;
using namespace protocols;


OPT_1GRP_KEY(File, edensity, alt_mapfile)
OPT_1GRP_KEY(Integer, edensity, nresbins)


///////////////////////////////////////////////////////////////////////////////
void
densityTools()
{
	using namespace protocols::moves;
	using namespace scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// outputs
	Size nresobins = option[ edensity::nresbins ]();
	utility::vector1< core::Real > resobins, mapI, modelI, modelmapFSC;

	// [1] map intensity statistics
	resobins = core::scoring::electron_density::getDensityMap().getResolutionBins(nresobins);
	mapI = core::scoring::electron_density::getDensityMap().getIntensities(nresobins);

	// [2] 2ary map stats (intensity + map v map FSC)
	bool usermap = false;
	if (option[ edensity::alt_mapfile ].user()) {
		usermap = true;
	}

	// [3] model-map stats (intensity + model v map FSC + RSCC + per-res corrleations)
	bool userpose = false;
	pose::Pose pose;
	if (option[ in::file::l ].user() || option[ in::file::s ].user()) {
		userpose = true;
		std::string pdbfile = basic::options::start_file();
		core::import_pose::pose_from_pdb( pose, pdbfile );

		// align to map
		protocols::electron_density::SetupForDensityScoringMoverOP dockindens( new protocols::electron_density::SetupForDensityScoringMover );
		core::scoring::electron_density::getDensityMap().set_nres( pose.total_residue() );
		dockindens->apply( pose );

		modelI = core::scoring::electron_density::getDensityMap().getIntensities( pose, option[ edensity::nresbins ]() );
		modelmapFSC = core::scoring::electron_density::getDensityMap().getFSC( pose, option[ edensity::nresbins ]() );
	}

	for (Size i=1; i<=resobins.size(); ++i) {
		std::cerr << resobins[i] << " " << mapI[i];
		if (userpose) std::cerr << " " << modelI[i] << " " << modelmapFSC[i];
		std::cerr << std::endl;
	}

	// [4] rescale maps to target intensity
	if (userpose) {
		utility::vector1< core::Real > rescale_factor(nresobins,0.0);
		for (Size i=1; i<=nresobins; ++i)
			rescale_factor[i] = modelI[i] / mapI[i];
		core::scoring::electron_density::getDensityMap().scaleIntensities( rescale_factor );
		core::scoring::electron_density::getDensityMap().writeMRC( "scale_modelI.mrc" );
		for (Size i=1; i<=nresobins; ++i)
			rescale_factor[i] *= modelmapFSC[i];
		core::scoring::electron_density::getDensityMap().scaleIntensities( rescale_factor );
		core::scoring::electron_density::getDensityMap().writeMRC( "scale_modelI_FSCwt.mrc" );
	}
}



///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
	// options, random initialization
	NEW_OPT(edensity::alt_mapfile, "alt mapfile", "");
	NEW_OPT(edensity::nresbins, "#reolution bins for statistics", 50);
	devel::init( argc, argv );
	densityTools();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	return 0;
}
