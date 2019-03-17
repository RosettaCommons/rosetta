// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.

/// @file    core/io/rcsb/ExperimentalTechnique.cc
/// @brief   Helper function definitions for ExperimentalTechnique.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @author  Matthew O'Meara


// Unit Header
#include <core/io/rcsb/ExperimentalTechnique.hh>

// Basic headers
#include <basic/Tracer.hh>


namespace core {
namespace io {
namespace rcsb {

static basic::Tracer TR( "core.io.rcsb.ExperimentalTechnique" );


std::string
experimental_technique_to_string( ExperimentalTechnique technique )
{
	using namespace std;

	string t;
	switch ( technique ) {
	case UNKNOWN_EXPDTA :           t = "UNKNOWN";                  break;
	case X_RAY_DIFFRACTION :        t = "X-RAY DIFFRACTION";        break;
	case FIBER_DIFFRACTION :        t = "FIBER DIFFRACTION";        break;
	case NEUTRON_DIFFRACTION :      t = "NEUTRON DIFFRACTION";      break;
	case ELECTRON_CRYSTALLOGRAPHY : t = "ELECTRON CRYSTALLOGRAPHY"; break;
	case ELECTRON_MICROSCOPY :      t = "ELECTRON MICROSCOPY";      break;
	case SOLID_STATE_NMR :          t = "SOLID-STATE NMR";          break;
	case SOLUTION_NMR :             t = "SOLUTION NMR";             break;
	case SOLUTION_SCATTERING :      t = "SOLUTION SCATTERING";      break;
	case THEORETICAL_MODEL :        t = "THEORETICAL MODEL";        break;

	case EPR :                      t = "EPR";                      break;

	case ELECTRON_DEFRACTION :
		t = "ELECTRON DEFRACTION";
		TR.Warning
			<< "Encountered obsolete experimental technique coding '"
			<< t << "'" << endl;
		break;

	case CRYO_ELECTRON_MICROSCOPY :
		t = "CRYO-ELECTRON MICROSCOPY";
		TR.Warning
			<< "Encountered obsolete experimental technique coding '"
			<< t << "'" << endl;
		break;

	case SOLUTION_SCATTERING_THEORETICAL_MODEL :
		t = "SOLUTION SCATTERING, THEORETICAL MODEL";
		TR.Warning
			<< "Encountered obsolete experimental technique coding '"
			<< t << "'" << endl;
		break;

	case FLORECENCE_TRANSFER :
		t = "FLORECENCE TRANSFER";
		TR.Warning
			<< "Encountered obsolete experimental technique coding '"
			<< t << "'" << endl;
		break;

	case NMR :
		t = "NMR";
		TR.Warning
			<< "Encountered obsolete experimental technique coding '"
			<< t << "'" << endl;
		break;

	default :
		TR.Error
			<< "Unrecognized experimental technique value '"
			<< technique << "'" << endl;
		t = "UNKNOWN";
	}
	return t;
}

ExperimentalTechnique
string_to_experimental_technique( std::string const & technique )
{
	using namespace std;

	if ( technique == "X-RAY DIFFRACTION" )        return X_RAY_DIFFRACTION;
	else if ( technique == "FIBER DIFFRACTION" )   return FIBER_DIFFRACTION;
	else if ( technique == "NEUTRON DIFFRACTION" ) return NEUTRON_DIFFRACTION;
	else if ( technique == "ELECTRON CRYSTALLOGRAPHY" ) {
		return ELECTRON_CRYSTALLOGRAPHY;
	} else if ( technique == "ELECTRON MICROSCOPY" ) return ELECTRON_MICROSCOPY;
	else if ( technique == "SOLID-STATE NMR" )     return SOLID_STATE_NMR;
	else if ( technique == "SOLUTION NMR" )        return SOLUTION_NMR;
	else if ( technique == "SOLUTION SCATTERING" ) return SOLUTION_SCATTERING;
	else if ( technique == "THEORETICAL MODEL" )   return THEORETICAL_MODEL;
	else if ( technique == "EPR" )                 return EPR;

	// Handle obsolete technique strings
	else if ( technique == "ELECTRON DEFRACTION" ) {
		TR.Warning
			<< "Encountered obsolete experimental technique string '"
			<< technique << "'" << endl;
		return ELECTRON_DEFRACTION;
	} else if ( technique == "CRYO-ELECTRON MICROSCOPY" ) {
		TR.Warning
			<< "Encountered obsolete experimental technique string '"
			<< technique << "'" << endl;
		return CRYO_ELECTRON_MICROSCOPY;
	} else if ( technique == "FLORECENCE TRANSFER" ) {
		TR.Warning
			<< "Encountered obsolete experimental technique string '"
			<< technique << "'" << endl;
		return FLORECENCE_TRANSFER;
	} else if ( technique == "NMR" ) {
		TR.Warning
			<< "Encountered obsolete experimental technique string '"
			<< technique << "'" << endl;
		return NMR;
	} else {
		TR.Error
			<< "Unrecognized experimental technique string '"
			<< technique << "'" << endl;
		return UNKNOWN_EXPDTA;
	}
	return UNKNOWN_EXPDTA;
}

}  // namespace rcsb
}  // namespace io
}  // namespace core
