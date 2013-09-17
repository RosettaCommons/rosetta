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
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/import_pose/import_pose.hh>
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
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

// C++ headers
#include <iostream>
#include <fstream>
#include <string>


using namespace core;
using namespace basic;
using namespace utility;
using namespace protocols;


OPT_1GRP_KEY(File, edensity, alt_mapfile)
OPT_1GRP_KEY(Integer, denstools, nresbins)
OPT_1GRP_KEY(Real, denstools, lowres)
OPT_1GRP_KEY(Real, denstools, hires)
OPT_1GRP_KEY(Boolean, denstools, rescale_map)
OPT_1GRP_KEY(Boolean, denstools, verbose)
OPT_1GRP_KEY(Boolean, denstools, dump_map_and_mask)
OPT_1GRP_KEY(Boolean, denstools, nomask)
OPT_1GRP_KEY(Boolean, denstools, perres)


using core::scoring::electron_density::poseCoords;
using core::scoring::electron_density::poseCoord;

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
		} else if ( (line[12] == ' ' && type!="Zn" && type!="Fe" && type!="ZN" && type!="FE")
								|| isdigit(type[1]) )
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
		if (atmid.substr(0,2) == "HD" || atmid.substr(0,2) == "HE" || atmid.substr(0,2) == "HG")
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
readPDBcoords(std::string filename, poseCoords &atmlist) {
	std::ifstream inpdb(filename.c_str());
	std::string buf;

	while (std::getline(inpdb, buf ) ) {
		if( buf.substr(0,4)!="ATOM" && buf.substr(0,6)!="HETATM") continue;
		poseCoord atom_i;

		atom_i.x_[0] = atof(buf.substr(30,8).c_str());
		atom_i.x_[1] = atof(buf.substr(38,8).c_str());
		atom_i.x_[2] = atof(buf.substr(46,8).c_str());
		atom_i.B_ = atof(buf.substr(60,6).c_str());

		atom_i.elt_ = name2elt( buf ); // horrible hacky logic mapping name->elt (could use PDB fields 76-77 if on by default)
		if (atom_i.elt_ == "H") continue;

		atmlist.push_back( atom_i );
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
	Size nresobins = option[ denstools::nresbins ]();
	utility::vector1< core::Real > resobins, mapI, mapIprime, modelI, modelSum, modelmapFSC;
	utility::vector1< core::Real > mapAltI, mapmapFSC;
	utility::vector1< core::Real > perResCC;
	Real rscc, fsc=0, mm_rscc, mm_fsc=0;

	// resolution limits for analysis
	core::Real hires = option[ denstools::hires ]();
	core::Real lowres = option[ denstools::lowres ]();
	if (hires == 0.0) { hires = core::scoring::electron_density::getDensityMap().maxNominalRes(); }

	runtime_assert( lowres > hires );

	hires = 1.0/hires;
	lowres = 1.0/lowres;

	// [0] load model, mask density
	bool userpose = (option[ in::file::s ].user());
	poseCoords pose;
	core::pose::Pose fullpose;  // needed for per-residue stats (if requested)
	std::string pdbfile;
	if (userpose) {
		pdbfile = basic::options::start_file();
		readPDBcoords( pdbfile, pose );
		if (!option[ denstools::nomask ]() ) {
			core::scoring::electron_density::getDensityMap().maskDensityMap( pose, 0 );
		}
	}

	// [1] map intensity statistics
	resobins = core::scoring::electron_density::getDensityMap().getResolutionBins(nresobins, lowres, hires);
	mapI = core::scoring::electron_density::getDensityMap().getIntensities(nresobins, lowres, hires);

	// [2] ALT map stats (intensity + map v map FSC)
	core::scoring::electron_density::ElectronDensity mapAlt;
	bool usermap = option[ edensity::alt_mapfile ].user();
	if (usermap) {
		std::string mapfile = option[ edensity::alt_mapfile ];
		core::Real mapreso = option[ edensity::mapreso ]();
		core::Real mapsampling = option[ edensity::grid_spacing ]();
		std::cerr << "Loading alternate density map " << mapfile << std::endl;
		mapAlt.readMRCandResize( mapfile , mapreso , mapsampling );

		if (userpose && !option[ denstools::nomask ]() ) {
			mapAlt.maskDensityMap( pose, 0 );
		}

		mapAltI = mapAlt.getIntensities(nresobins, lowres, hires);
		mapmapFSC = core::scoring::electron_density::getDensityMap().getFSC( mapAlt.data(), nresobins, lowres, hires );
		mm_rscc = core::scoring::electron_density::getDensityMap().getRSCC( mapAlt.data() );

		for (Size i=1; i<=resobins.size(); ++i)
			mm_fsc+=mapmapFSC[i];
		mm_fsc /= resobins.size();
	}

	// [3] model-map stats (intensity + model v map FSC + RSCC + per-res corrleations)
	if (userpose) {
		if (option[ denstools::perres ]()) {
			core::chemical::ResidueTypeSetCAP rsd_set;
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
			core::import_pose::pose_from_pdb( fullpose, pdbfile );

			core::Size nres = fullpose.total_residue();
			perResCC.resize( nres, 0.0 );

			protocols::electron_density::SetupForDensityScoringMoverOP dockindens
				( new protocols::electron_density::SetupForDensityScoringMover );
			dockindens->apply( fullpose );

			core::scoring::electron_density::getDensityMap().set_nres( nres );
			core::scoring::electron_density::getDensityMap().setScoreWindowContext( true );
			core::scoring::electron_density::getDensityMap().setWindow( 3 );

			core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
			scorefxn->set_weight( core::scoring::elec_dens_window, 1.0 );
			(*scorefxn)(fullpose);

			for (int r=1; r<=nres; ++r) {
				perResCC[r] = core::scoring::electron_density::getDensityMap().matchRes( r , fullpose.residue(r), fullpose, NULL , false);
			}
		}

		core::scoring::electron_density::getDensityMap().getIntensities( pose, nresobins, lowres, hires, modelI );
		modelmapFSC = core::scoring::electron_density::getDensityMap().getFSC( pose, nresobins, lowres, hires );
		rscc = core::scoring::electron_density::getDensityMap().getRSCC( pose );
		for (Size i=1; i<=resobins.size(); ++i)
			fsc+=modelmapFSC[i];
		fsc /= resobins.size();
	}

	// [5] optionally: rescale maps to target intensity
	if (option[ denstools::rescale_map ]()) {
		if (userpose) {
			utility::vector1< core::Real > rescale_factor(nresobins,0.0);
			for (Size i=1; i<=nresobins; ++i)
				rescale_factor[i] = sqrt(modelI[i] / mapI[i]);
			core::scoring::electron_density::getDensityMap().scaleIntensities( rescale_factor, lowres, hires );
			core::scoring::electron_density::getDensityMap().writeMRC( "scale_modelI.mrc" );
			mapIprime = core::scoring::electron_density::getDensityMap().getIntensities(nresobins, lowres, hires);

			//for (Size i=1; i<=nresobins; ++i)
			//	rescale_factor[i] *= modelmapFSC[i];
			//core::scoring::electron_density::getDensityMap().scaleIntensities( rescale_factor, lowres, hires );
			//core::scoring::electron_density::getDensityMap().writeMRC( "scale_modelI_FSCwt.mrc" );
		} else if (usermap) {
			utility::vector1< core::Real > rescale_factor(nresobins,0.0);
			for (Size i=1; i<=nresobins; ++i)
				rescale_factor[i] = sqrt(mapAltI[i] / mapI[i]);
			core::scoring::electron_density::getDensityMap().scaleIntensities( rescale_factor, lowres, hires );
			core::scoring::electron_density::getDensityMap().writeMRC( "scale_altmapI.mrc" );
			mapIprime = core::scoring::electron_density::getDensityMap().getIntensities(nresobins, lowres, hires);
		}
	}

	// verbose
	if ( option[ denstools::verbose ]()) {
			std::cerr << "1/res  I_map1";
			if (usermap)  std::cerr << "  I_map2  FSC_map1_map2" ;
			if (userpose) std::cerr << "  I_model  FSC_map1_model";
			std::cerr << std::endl;
		for (Size i=1; i<=resobins.size(); ++i) {
			std::cerr << resobins[i] << " " << mapI[i];
			if (usermap)  std::cerr << " " << mapAltI[i] << " " << mapmapFSC[i];
			if (userpose) std::cerr << " " << modelI[i] << " " << modelmapFSC[i];
			std::cerr << std::endl;
		}
	}

	// dump
	if (userpose && option[ denstools::dump_map_and_mask ]()) {
		core::scoring::electron_density::getDensityMap().writeMRC( "model_dens.mrc", true, false );
		core::scoring::electron_density::getDensityMap().writeMRC( "model_mask.mrc", false, true  );
	}

	if (userpose && option[ denstools::perres ]()) {
		for (int r=1; r<=perResCC.size(); ++r) {
			std::cerr << "residue " << r << "  cc=" << perResCC[r] << std::endl;
		}
	}

	// compact
	if (userpose) {
		std::cerr << pdbfile << " " << fsc << " " << rscc << std::endl;
	}
	if (usermap) {
		std::cerr << option[ edensity::alt_mapfile ]() << " " << mm_fsc << " " << mm_rscc << std::endl;
	}
}



///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
	// options, random initialization
	NEW_OPT(denstools::lowres, "lowres", 500.0);
	NEW_OPT(denstools::hires, "hires", 0.0);
	NEW_OPT(edensity::alt_mapfile, "alt mapfile", "");
	NEW_OPT(denstools::nresbins, "#resolution bins for statistics", 20);
	NEW_OPT(denstools::rescale_map, "scale map I == model I?", false);
	NEW_OPT(denstools::dump_map_and_mask, "dump_map_and_mask", false);
	NEW_OPT(denstools::nomask, "nomask", false);
	NEW_OPT(denstools::perres, "output per-residue stats", false);
	NEW_OPT(denstools::verbose, "extra output", false);
	devel::init( argc, argv );
	densityTools();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}
	return 0;
}
