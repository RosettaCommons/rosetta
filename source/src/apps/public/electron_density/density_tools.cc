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
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/hybridization/FragmentBiasAssigner.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>

#include <basic/options/util.hh>
#include <devel/init.hh>


#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>

#include <ObjexxFCL/string.functions.hh>
#include <core/kinematics/Jump.hh>
#include <basic/Tracer.hh>

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
OPT_1GRP_KEY(Integer, denstools, rscc_wdw)
OPT_1GRP_KEY(Real, denstools, lowres)
OPT_1GRP_KEY(Real, denstools, hires)
OPT_1GRP_KEY(Real, denstools, truncate_lowres)
OPT_1GRP_KEY(Real, denstools, truncate_hires)
OPT_1GRP_KEY(Boolean, denstools, truncate_map)
OPT_1GRP_KEY(Boolean, denstools, verbose)
OPT_1GRP_KEY(Boolean, denstools, super_verbose)
OPT_1GRP_KEY(Boolean, denstools, dump_map_and_mask)
OPT_1GRP_KEY(Boolean, denstools, mask)
OPT_1GRP_KEY(Boolean, denstools, mask_modelmap)
OPT_1GRP_KEY(Real, denstools, mask_resolution)
OPT_1GRP_KEY(Boolean, denstools, perres)
OPT_1GRP_KEY(Boolean, denstools, bin_squared)
OPT_1GRP_KEY(Boolean, denstools, maskonly)
OPT_1GRP_KEY(Boolean, denstools, cutonly)


using core::scoring::electron_density::poseCoords;
using core::scoring::electron_density::poseCoord;


static basic::Tracer TR( "density_tools" );

// map atom names to elements
//   loosely based on openbabel logic
std::string
name2elt( std::string line ) {
	std::string atmid = line.substr(12,4);
	while ( !atmid.empty() && atmid[0] == ' ' ) atmid = atmid.substr(1,atmid.size()-1);
	while ( !atmid.empty() && atmid[atmid.size()-1] == ' ' ) atmid = atmid.substr(0,atmid.size()-1);
	std::string resname = line.substr(17,3);
	while ( !resname.empty() && resname[0] == ' ' ) resname = resname.substr(1,resname.size()-1);
	while ( !resname.empty() && resname[resname.size()-1] == ' ' ) resname = resname.substr(0,resname.size()-1);

	std::string type;

	if ( line.substr(0,4) == "ATOM" ) {
		type = atmid.substr(0,2);
		if ( isdigit(type[0]) ) {
			// sometimes non-standard files have, e.g 11HH
			if ( !isdigit(type[1]) ) type = atmid.substr(1,1);
			else type = atmid.substr(2,1);
		} else if ( (line[12] == ' ' && type!="Zn" && type!="Fe" && type!="ZN" && type!="FE")
				|| isdigit(type[1]) ) {
			type = atmid.substr(0,1);     // one-character element
		}

		if ( resname.substr(0,2) == "AS" || resname[0] == 'N' ) {
			if ( atmid == "AD1" ) type = "O";
			if ( atmid == "AD2" ) type = "N";
		}
		if ( resname.substr(0,3) == "HIS" || resname[0] == 'H' ) {
			if ( atmid == "AD1" || atmid == "AE2" ) type = "N";
			if ( atmid == "AE1" || atmid == "AD2" ) type = "C";
		}
		if ( resname.substr(0,2) == "GL" || resname[0] == 'Q' ) {
			if ( atmid == "AE1" ) type = "O";
			if ( atmid == "AE2" ) type = "N";
		}
		if ( atmid.substr(0,2) == "HH" ) { // ARG
			type = "H";
		}
		if ( atmid.substr(0,2) == "HD" || atmid.substr(0,2) == "HE" || atmid.substr(0,2) == "HG" ) {
			type = "H";
		}
	} else {
		if ( isalpha(atmid[0]) ) {
			if ( atmid.size() > 2 && (atmid[2] == '\0' || atmid[2] == ' ') ) {
				type = atmid.substr(0,2);
			} else if ( atmid[0] == 'A' ) { // alpha prefix
				type = atmid.substr(1, atmid.size() - 1);
			} else {
				type = atmid.substr(0,1);
			}
		} else if ( atmid[0] == ' ' ) {
			type = atmid.substr(1,1); // one char element
		} else {
			type = atmid.substr(1,2);
		}

		if ( atmid == resname ) {
			type = atmid;
			if ( type.size() == 2 ) type[1] = toupper(type[1]);
		} else if ( resname == "ADR" || resname == "COA" || resname == "FAD" ||
				resname == "GPG" || resname == "NAD" || resname == "NAL" ||
				resname == "NDP" || resname == "ABA" )  {
			if ( type.size() > 1 ) type = type.substr(0,1);
		} else if ( isdigit(type[0]) ) {
			type = type.substr(1,1);
		} else if ( type.size() > 1 && isdigit(type[1]) ) {
			type = type.substr(0,1);
		} else if ( type.size() > 1 && isalpha(type[1]) ) {
			if ( type[0] == 'O' && type[1] == 'H' ) {
				type = type.substr(0,1); // no "Oh" element (e.g. 1MBN)
			} else if ( islower(type[1]) ) {
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

	while ( std::getline(inpdb, buf ) ) {
		if ( buf.substr(0,4)!="ATOM" && buf.substr(0,6)!="HETATM" ) continue;
		poseCoord atom_i;

		atom_i.x_[0] = atof(buf.substr(30,8).c_str());
		atom_i.x_[1] = atof(buf.substr(38,8).c_str());
		atom_i.x_[2] = atof(buf.substr(46,8).c_str());
		atom_i.B_ = atof(buf.substr(60,6).c_str());

		atom_i.elt_ = name2elt( buf ); // horrible hacky logic mapping name->elt (could use PDB fields 76-77 if on by default)
		if ( atom_i.elt_ == "H" ) continue;

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
	bool bin_squared =  option[ denstools::bin_squared ]();
	utility::vector1< core::Real > resobins, mapI, maskedmapI, modelI, maskI;
	utility::vector1< core::Size > resobin_counts;
	utility::vector1< core::Real > perResCC;
	std::map< core::Size, core::Real > perResStrain;

	utility::vector1< core::Real > mapmapFSC, maskedMapMapFSC;
	utility::vector1< core::Real > modelmapFSC, maskedModelMapFSC;

	// resolution limits for analysis
	core::Real hires = option[ denstools::hires ]();
	core::Real lowres = option[ denstools::lowres ]();
	core::Real truncate_hires = option[ denstools::truncate_hires ]();
	core::Real truncate_lowres = option[ denstools::truncate_lowres ]();
	core::Real mask_resolution = option[ denstools::mask_resolution ]();
	bool perres =option[ denstools::perres ]();

	ObjexxFCL::FArray3D< double > rhoC, rhoMask, rhoOmask, rhoO2mask;
	ObjexxFCL::FArray3D< std::complex<double> > FrhoC, FrhoMask, FrhoCmask, FrhoOmask, FrhoO2mask;
	ObjexxFCL::FArray3D< std::complex<double> > FrhoO, FrhoO2;

	if ( hires == 0.0 ) hires = core::scoring::electron_density::getDensityMap().maxNominalRes();
	if ( truncate_hires == 0.0 ) truncate_hires = core::scoring::electron_density::getDensityMap().maxNominalRes();

	runtime_assert( lowres > hires );
	runtime_assert( truncate_lowres > truncate_hires );

	hires = 1.0/hires;
	lowres = 1.0/lowres;
	truncate_lowres = 1.0/truncate_lowres;
	truncate_hires = 1.0/truncate_hires;

	//  fft
	TR.Trace << "Stage 1: FFT rho_obs" << std::endl;
	numeric::fourier::fft3(core::scoring::electron_density::getDensityMap().get_data(), FrhoO);
	core::scoring::electron_density::getDensityMap().getIntensities(FrhoO, nresobins, lowres, hires, mapI, bin_squared);

	// load model, mask density
	TR.Trace << "Stage 2: model processing" << std::endl;
	bool userpose = (option[ in::file::s ].user());

	core::scoring::electron_density::getDensityMap().getResolutionBins(nresobins, lowres, hires, resobins, resobin_counts, bin_squared);

	// Truncate map only
	if ( option[ denstools::truncate_map ] ) {
		utility::vector1< core::Real > trunc = mapI;
		for ( int i=1; i<=(int)mapI.size(); ++i ) {
			if ( resobins[i] <= truncate_lowres || resobins[i] >= truncate_hires ) {
				trunc[i] = 0;
			} else {
				trunc[i] = 1;
			}
		}

		core::scoring::electron_density::getDensityMap().scaleIntensities( trunc, lowres, hires, bin_squared );
		core::scoring::electron_density::getDensityMap().writeMRC( "map_trunc.mrc" );
		TR.Trace << "Wrote truncated map" << std::endl;
		return;
	}

	poseCoords pose;
	core::pose::Pose fullpose;  // needed for per-residue stats (if requested)
	std::string pdbfile;
	if ( userpose ) {
		pdbfile = basic::options::start_file();
		readPDBcoords( pdbfile, pose );

		TR.Trace << "       : calc rho_c" << std::endl;
		core::scoring::electron_density::getDensityMap().calcRhoC( pose, mask_resolution, rhoC, rhoMask );
	}

	// Mask map only
	if ( option[ denstools::maskonly ] ) {
		if ( !userpose ) {
			utility_exit_with_message(" -maskonly given but no PDB file provided!");
		}

		rhoOmask = rhoMask;
		for ( int i=0; i<rhoC.u1()*rhoC.u2()*rhoC.u3(); ++i ) {
			rhoOmask[i] *= core::scoring::electron_density::getDensityMap().get_data()[i];
		}

		core::scoring::electron_density::getDensityMap().set_data(rhoOmask);
		core::scoring::electron_density::getDensityMap().writeMRC( option[ edensity::mapfile]+"_masked.mrc" );
		return;
	}

	// Cut map only
	if ( option[ denstools::cutonly ] ) {
		if ( !userpose ) {
			utility_exit_with_message(" -maskonly given but no PDB file provided!");
		}

		rhoOmask = rhoMask;
		for ( int i=0; i<rhoC.u1()*rhoC.u2()*rhoC.u3(); ++i ) {
			rhoOmask[i] = (1-rhoOmask[i])*core::scoring::electron_density::getDensityMap().get_data()[i];
		}

		core::scoring::electron_density::getDensityMap().set_data(rhoOmask);
		core::scoring::electron_density::getDensityMap().writeMRC( option[ edensity::mapfile]+"_cut.mrc" );
		return;
	}

	// collect masked map statistics
	if ( userpose ) {
		TR.Trace << "       : FFT" << std::endl;
		numeric::fourier::fft3(rhoC, FrhoC);
		numeric::fourier::fft3(rhoMask, FrhoMask);

		// apply mask
		TR.Trace << "       : mask map" << std::endl;
		rhoOmask = rhoMask;
		for ( int i=0; i<rhoC.u1()*rhoC.u2()*rhoC.u3(); ++i ) {
			rhoOmask[i] *= core::scoring::electron_density::getDensityMap().get_data()[i];
		}
		numeric::fourier::fft3(rhoOmask, FrhoOmask);

		// intensities
		TR.Trace << "       : get intensities" << std::endl;
		core::scoring::electron_density::getDensityMap().getIntensities(FrhoC, nresobins, lowres, hires, modelI, bin_squared);
		core::scoring::electron_density::getDensityMap().getIntensities(FrhoMask, nresobins, lowres, hires, maskI, bin_squared);
		core::scoring::electron_density::getDensityMap().getIntensities(FrhoOmask, nresobins, lowres, hires, maskedmapI, bin_squared);
	}


	// map vs map stats
	TR.Trace << "Stage 3: map vs map" << std::endl;
	core::scoring::electron_density::ElectronDensity mapAlt;
	bool usermap = option[ edensity::alt_mapfile ].user();
	if ( usermap ) {
		std::string mapfile = option[ edensity::alt_mapfile ];
		core::Real mapreso = option[ edensity::mapreso ]();
		core::Real mapsampling = option[ edensity::grid_spacing ]();
		TR.Trace << "       : load alt_map + FFT" << std::endl;
		mapAlt.readMRCandResize( mapfile , mapreso , mapsampling );
		numeric::fourier::fft3(mapAlt.get_data(), FrhoO2);

		if ( userpose ) {
			TR.Trace << "       : mask alt_map" << std::endl;
			rhoO2mask = rhoMask;
			for ( int i=0; i<rhoC.u1()*rhoC.u2()*rhoC.u3(); ++i ) {
				rhoO2mask[i] *= mapAlt.get_data()[i];
			}
			numeric::fourier::fft3(rhoO2mask, FrhoO2mask);
		}

		// unmasked
		TR.Trace << "       : get FSCs" << std::endl;
		core::scoring::electron_density::getDensityMap().getFSC(FrhoO, FrhoO2, nresobins, lowres, hires, mapmapFSC, bin_squared);

		// masked
		if ( userpose ) {
			core::scoring::electron_density::getDensityMap().getFSC(FrhoOmask, FrhoO2mask, nresobins, lowres, hires, mapmapFSC, bin_squared);
		}
	}

	// [3] model-map stats (intensity + model v map FSC + RSCC + per-res corrleations)
	Real modelMapFSCsum=0, maskedModelMapFSCsum=0, RSCC=0;
	core::Size bincountsum=0;

	TR.Trace << "Stage 4: model vs map" << std::endl;
	if ( userpose ) {
		if ( perres ) {
			TR.Trace << "       : per-res" << std::endl;
			core::chemical::ResidueTypeSetCAP rsd_set;
			rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
			core::import_pose::pose_from_file( fullpose, pdbfile , core::import_pose::PDB_file);

			core::Size nres = fullpose.size();
			perResCC.resize( nres, 0.0 );
			for (Size i= 1; i <=nres; ++i) perResStrain[i] = 0.0;

			protocols::electron_density::SetupForDensityScoringMoverOP dockindens( new protocols::electron_density::SetupForDensityScoringMover );
			dockindens->apply( fullpose );

			core::scoring::electron_density::getDensityMap().set_nres( nres );
			core::scoring::electron_density::getDensityMap().setScoreWindowContext( true );
			core::scoring::electron_density::getDensityMap().setWindow( option[ denstools::rscc_wdw ]() );

			core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
			scorefxn->set_weight( core::scoring::elec_dens_window, 1.0 );
			(*scorefxn)(fullpose);

			for ( core::uint r = 1; r <= nres; ++r ) {
				perResCC[r] = core::scoring::electron_density::getDensityMap().matchRes( r , fullpose.residue(r), fullpose, nullptr , false);
			}
			protocols::hybridization::FragmentBiasAssigner fa(fullpose);
			fa.automode_scores( fullpose, perResStrain );
		}

		TR.Trace << "       : FSCs" << std::endl;
		core::scoring::electron_density::getDensityMap().getFSC(FrhoO, FrhoC, nresobins, lowres, hires, modelmapFSC, bin_squared);
		core::scoring::electron_density::getDensityMap().getFSC(FrhoOmask, FrhoC, nresobins, lowres, hires, maskedModelMapFSC, bin_squared);

		TR.Trace << "       : RSCC" << std::endl;
		RSCC = core::scoring::electron_density::getDensityMap().getRSCC( rhoC, rhoMask );

		// sum FSC over reso bins
		for ( Size i=1; i<=resobins.size(); ++i ) {
			modelMapFSCsum += resobin_counts[i] * modelmapFSC[i];
			maskedModelMapFSCsum += resobin_counts[i] * maskedModelMapFSC[i];
			bincountsum += resobin_counts[i];
		}
		modelMapFSCsum /= bincountsum;
		maskedModelMapFSCsum /= bincountsum;
	}

	// verbose
	if ( option[ denstools::verbose ]() ) {
		TR << "res count";
		if ( mapI.size() ) TR << " I_map";
		if ( maskedmapI.size() ) TR << " I_map_masked";
		if ( modelI.size() ) TR << " I_model" ;
		if ( maskI.size() ) TR << " I_mask" ;
		if ( mapmapFSC.size() )  TR << " FSC_mapmap";
		if ( modelmapFSC.size() )  TR << " FSC_modelmap";
		if ( maskedModelMapFSC.size() )  TR << " FSC_maskedmodelmap";
		TR << std::endl;
		for ( Size i=1; i<=resobins.size(); ++i ) {
			TR << resobins[i] << " " << resobin_counts[i];
			if ( mapI.size() ) TR << " " << mapI[i];
			if ( maskedmapI.size() ) TR << " " << maskedmapI[i];
			if ( modelI.size() ) TR << " " << modelI[i];
			if ( maskI.size() ) TR << " " << maskI[i];
			if ( mapmapFSC.size() )  TR << " " << mapmapFSC[i];
			if ( modelmapFSC.size() )  TR << " " << modelmapFSC[i];
			if ( maskedModelMapFSC.size() )  TR << " " << maskedModelMapFSC[i];
			TR << std::endl;
		}
	}

	// dump
	if ( userpose && option[ denstools::dump_map_and_mask ]() ) {
		core::scoring::electron_density::getDensityMap().set_data(rhoOmask);
		core::scoring::electron_density::getDensityMap().writeMRC( "map_masked.mrc" );

		core::scoring::electron_density::getDensityMap().set_data(rhoC);
		core::scoring::electron_density::getDensityMap().writeMRC( "model_dens.mrc" );

		core::scoring::electron_density::getDensityMap().set_data(rhoMask);
		core::scoring::electron_density::getDensityMap().writeMRC( "model_mask.mrc" );
	}

	if ( userpose && option[ denstools::perres ]() ) {
		for ( core::uint r = 1; r <= perResCC.size(); ++r ) {
			if ( fullpose.pdb_info() ) {
				core::pose::PDBInfoOP pdbinfo = fullpose.pdb_info();
				TR << "PERRESCC residue " << fullpose.residue(r).name3() << " " << pdbinfo->chain(r) << " " << pdbinfo->number(r) << pdbinfo->icode(r)
					<< " " << perResCC[r] << " " << perResStrain[r] << std::endl;
			} else {
				TR << "PERRESCC residue " << r << " " << perResCC[r] << " " << perResStrain[r] << std::endl;
			}
		}
	}

	// compact
	if ( userpose ) {
		TR << pdbfile << " RSCC/FSC/FSCmask: " << RSCC << " " << modelMapFSCsum << " " << maskedModelMapFSCsum;
		TR << std::endl;
	}
}


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
		// options, random initialization
		NEW_OPT(denstools::lowres, "low res limit", 1000.0);
		NEW_OPT(denstools::hires, "high res limit", 0.0);
		NEW_OPT(edensity::alt_mapfile, "alt mapfile", "");
		NEW_OPT(denstools::nresbins, "# resolution bins for statistics", 200);
		NEW_OPT(denstools::truncate_map, "dump map with map I scaled to to model I?", false);
		NEW_OPT(denstools::truncate_lowres, "low res truncation", 1000.0);
		NEW_OPT(denstools::truncate_hires, "high res truncation", 0.0);
		NEW_OPT(denstools::dump_map_and_mask, "dump mrc of rho_calc and eps_calc", false);
		NEW_OPT(denstools::mask, "mask all calcs", false);
		NEW_OPT(denstools::mask_modelmap, "mask model-map calcs only", false);
		NEW_OPT(denstools::mask_resolution, "radius for masking", 0.0);
		NEW_OPT(denstools::perres, "output per-residue stats", false);
		NEW_OPT(denstools::verbose, "dump extra output", false);
		NEW_OPT(denstools::super_verbose, "dump a lot of extra output", false);
		NEW_OPT(denstools::bin_squared, "bin uniformly in 1/res^2 (default bins 1/res)", false);
		NEW_OPT(denstools::rscc_wdw, "sliding window to calculate rscc", 3);
		NEW_OPT(denstools::maskonly, "mask", false);
		NEW_OPT(denstools::cutonly, "cut", false);

		devel::init( argc, argv );
		densityTools();
	} catch (utility::excn::Exception const & e ) {
		TR.Trace << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
