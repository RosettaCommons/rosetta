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
#include <fstream>
#include <string>


using namespace core;
using namespace basic;
using namespace utility;
using namespace protocols;


OPT_1GRP_KEY(File, edensity, alt_mapfile)
OPT_1GRP_KEY(Integer, edensity, nresbins)

// fpd
//   very quickly read heavy atom positions from a PDB for a possibly massive massive PDB file
typedef utility::vector1< std::pair< numeric::xyzVector< core::Real >, std::string > > lightPose;

// map atom names to elements
//   loosely based on openbabel logic
std::string
name2elt( std::string line ) {
	std::string atmid = line.substr(12,4);
	while (!atmid.empty() && atmid[0] == ' ') atmid = atmid.substr(1,atmid.size()-1);
	while (!atmid.empty() && atmid[atmid.size()-1] == ' ') atmid = atmid.substr(0,atmid.size()-1);
	std::string resname = line.substr(17,3);
	while (!resname.empty() && resname[0] == ' ') resname = resname.substr(1,resname.size()-1);
	while (!resname.empty() && resname[resname.size()-1] == ' ') resname = resname.substr(0,resname.size()-1);

	std::string type;

	if (line.substr(0,4) == "ATOM") {
		type = atmid.substr(0,2);
		if (isdigit(type[0])) {
			// sometimes non-standard files have, e.g 11HH
			if (!isdigit(type[1])) type = atmid.substr(1,1);
			else type = atmid.substr(2,1);
		} else if (line[12] == ' ' &&
				type!="Zn" && type!="Fe" && type!="ZN" && type!="FE" ||
				std::isdigit(type[1]))	//type[1] is digit in Platon
			type = atmid.substr(0,1);     // one-character element

			if (resname.substr(0,2) == "AS" || resname[0] == 'N') {
				if (atmid == "AD1") type = "O";
				if (atmid == "AD2") type = "N";
			}
			if (resname.substr(0,3) == "HIS" || resname[0] == 'H') {
				if (atmid == "AD1" || atmid == "AE2") type = "N";
				if (atmid == "AE1" || atmid == "AD2") type = "C";
			}
			if (resname.substr(0,2) == "GL" || resname[0] == 'Q') {
				if (atmid == "AE1") type = "O";
				if (atmid == "AE2") type = "N";
			}
			if (atmid.substr(0,2) == "HH") // ARG
          type = "H";
	} else {
		if (isalpha(atmid[0])) {
			if (atmid.size() > 2 && (atmid[2] == '\0' || atmid[2] == ' '))
				type = atmid.substr(0,2);
			else if (atmid[0] == 'A') // alpha prefix
				type = atmid.substr(1, atmid.size() - 1);
			else
				type = atmid.substr(0,1);
		} else if (atmid[0] == ' ')
			type = atmid.substr(1,1); // one char element
		else
			type = atmid.substr(1,2);

		if (atmid == resname) {
			type = atmid;
			if (type.size() == 2) type[1] = toupper(type[1]);
		} else if (resname == "ADR" || resname == "COA" || resname == "FAD" ||
				resname == "GPG" || resname == "NAD" || resname == "NAL" ||
				resname == "NDP" || resname == "ABA")  {
			if (type.size() > 1) type = type.substr(0,1);
		} else if (isdigit(type[0])) {
			type = type.substr(1,1);
		} else if (type.size() > 1 && isdigit(type[1])) {
			type = type.substr(0,1);
		} else if (type.size() > 1 && isalpha(type[1])) {
			if (type[0] == 'O' && type[1] == 'H')
				type = type.substr(0,1); // no "Oh" element (e.g. 1MBN)
			else if(islower(type[1])) {
				type[1] = toupper(type[1]);
			}
		}
	}
	return type;
}


// quick and dirty PDB read where we only care about heavyatom locations and atom ids
void
readPDBcoords(std::string filename, lightPose &atmlist) {
	std::ifstream inpdb(filename.c_str());
	std::string buf;

	while (std::getline(inpdb, buf ) ) {
		if( buf.substr(0,4)!="ATOM" && buf.substr(0,6)!="HETATM") continue;

		numeric::xyzVector< core::Real > X;
		X[0] = atof(buf.substr(30,8).c_str());
		X[1] = atof(buf.substr(38,8).c_str());
		X[2] = atof(buf.substr(46,8).c_str());

		// horrible hacky logic mapping name->elt (could use PDB fields 76-77 if Rosetta used them by default)
		std::string elt = name2elt( buf );
		if (elt == "H") continue;

		atmlist.push_back( std::make_pair(X,elt) );
	}
}



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
	lightPose pose;
	if (option[ in::file::l ].user() || option[ in::file::s ].user()) {
		userpose = true;
		std::string pdbfile = basic::options::start_file();
		readPDBcoords( pdbfile, pose );

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
