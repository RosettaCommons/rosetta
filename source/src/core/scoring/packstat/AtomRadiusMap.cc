// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packstat/AtomRadiusMap.cc
///
/// @brief
/// @author will sheffler


// Unit header or inline function header
#include <core/scoring/packstat/AtomRadiusMap.hh>
#include <basic/options/option.hh>

// option key includes
#include <basic/options/keys/packstat.OptionKeys.gen.hh>

#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace packstat {

using namespace std;

core::Real HYDROGEN_RADIUS = 0.0;

typedef pair<string const,string const> P;

std::string ToUpper( std::string s ) {
#ifdef WIN32
	for ( size_t i = 0; i < s.size(); ++i ) s[i] = toupper( s[i] );
#else
	for ( char & i : s ) i = std::toupper( i );
#endif
	// std::transform( s.begin(), s.end(), s.begin(), std::toupper );
	return s;
}

void type_map_add( map<P,PackstatReal> & type_map, string type, string res, PackstatReal rad ) {
	type = core::scoring::packstat::ToUpper(type);
	res  = core::scoring::packstat::ToUpper(res );
	type_map[ P(type,res) ] = rad;
}

AtomRadiusMap::AtomRadiusMap(
	int include_water // -1 (default) means read from global options
) {

	if ( include_water < 0 ) {
		using namespace basic::options;
		include_water = option[ basic::options::OptionKeys::packstat::include_water ];
	}

	// complete guesses
	type_map_add( type_map_, "UNX", "*", 2 ); // unknown...
	type_map_add( type_map_, "UNK", "*", 2 ); // say 2?

	type_map_add( type_map_, "XD", "*", 0 ); // deuterium?
	type_map_add( type_map_, "D", "*",  0 ); //
	type_map_add( type_map_, "1D", "*", 0 ); //
	type_map_add( type_map_, "2D", "*", 0 ); //
	type_map_add( type_map_, "3D", "*", 0 ); //
	type_map_add( type_map_, "*DO5","*",0 ); //

	// from Chimera
	type_map_add( type_map_, "al", "*", 0.54 );
	type_map_add( type_map_, "as", "*", 0.58 );
	type_map_add( type_map_, "au", "*", 1.37 );
	type_map_add( type_map_, "ba", "*", 1.35 );
	type_map_add( type_map_, "be", "*", 0.45 );
	type_map_add( type_map_, "bi", "*", 1.03 );
	type_map_add( type_map_, "ca", "*", 1.0  );
	type_map_add( type_map_, "cd", "*", 0.95 );
	type_map_add( type_map_, "co", "*", 0.65 );
	type_map_add( type_map_, "cr", "*", 0.73 );
	type_map_add( type_map_, "cs", "*", 1.67 );
	type_map_add( type_map_, "cu", "*", 0.73 );
	type_map_add( type_map_, "ga", "*", 0.62 );
	type_map_add( type_map_, "ge", "*", 0.73 );
	type_map_add( type_map_, "hg", "*", 1.02 );
	type_map_add( type_map_, "li", "*", 0.76 );
	type_map_add( type_map_, "mg", "*", 0.72 );
	type_map_add( type_map_, "mn", "*", 0.83 );
	type_map_add( type_map_, "mo", "*", 0.69 );
	type_map_add( type_map_, "na", "*", 1.02 );
	type_map_add( type_map_, "ni", "*", 0.69 );
	type_map_add( type_map_, "pb", "*", 1.19 );
	type_map_add( type_map_, "pd", "*", 0.86 );
	type_map_add( type_map_, "pt", "*", 0.8  );
	type_map_add( type_map_, "rb", "*", 1.52 );
	type_map_add( type_map_, "sb", "*", 0.76 );
	type_map_add( type_map_, "sc", "*", 0.75 );
	type_map_add( type_map_, "sn", "*", 0.69 );
	type_map_add( type_map_, "sr", "*", 1.18 );
	type_map_add( type_map_, "tc", "*", 0.65 );
	type_map_add( type_map_, "ti", "*", 0.86 );
	type_map_add( type_map_, "v" , "*", 0.79 );
	type_map_add( type_map_, "zn", "*", 0.74 );
	type_map_add( type_map_, "zr", "*", 0.72 );

	// random web searches
	type_map_add( type_map_, "b" , "*", 1.7  ); // boron?
	type_map_add( type_map_, "re", "*", 1.37 ); // Rhenium?
	type_map_add( type_map_, "rh", "*", 1.34 ); // Rhodium ?
	type_map_add( type_map_, "ru", "*", 1.32 ); // Ruthenium?
	type_map_add( type_map_, "tb", "*", 1.76 ); // terbium
	type_map_add( type_map_, "tl", "*", 1.96 ); // thallium
	type_map_add( type_map_, "xe", "*", 2.16 ); // xenon
	type_map_add( type_map_, "y",  "*", 1.77 ); // Yttrium
	type_map_add( type_map_, "lu", "*", 1.71 ); // Lutetium
	type_map_add( type_map_, "la", "*", 1.87 ); // Lanthanum
	type_map_add( type_map_, "U" , "*", 1.86 ); // Uranium!
	type_map_add( type_map_, "au", "*", 1.66 ); // gold
	type_map_add( type_map_, "ta", "*", 1.43 ); // tantalum
	type_map_add( type_map_, "eu", "*", 1.99 ); // europium
	type_map_add( type_map_, "gd", "*", 1.78 ); // Gadolinium

	// from ccdc.cam.ac.uk quest
	type_map_add( type_map_, "Ag", "*", 1.72 );
	type_map_add( type_map_, "Ar", "*", 1.88 );
	type_map_add( type_map_, "As", "*", 1.85 );
	type_map_add( type_map_, "Au", "*", 1.66 );
	type_map_add( type_map_, "Br", "*", 1.85 );
	type_map_add( type_map_, "Cd", "*", 1.58 );
	type_map_add( type_map_, "Cl", "*", 1.75 );
	type_map_add( type_map_, "Cu", "*", 1.40 );
	type_map_add( type_map_, "Ga", "*", 1.87 );
	type_map_add( type_map_, "He", "*", 1.40 );
	type_map_add( type_map_, "Hg", "*", 1.55 );
	type_map_add( type_map_, "In", "*", 1.93 );
	type_map_add( type_map_, "Kr", "*", 2.02 );
	type_map_add( type_map_, "Li", "*", 1.82 );
	type_map_add( type_map_, "Mg", "*", 1.73 );
	type_map_add( type_map_, "Na", "*", 2.27 );
	type_map_add( type_map_, "Ne", "*", 1.54 );
	type_map_add( type_map_, "Ni", "*", 1.63 );
	type_map_add( type_map_, "Pb", "*", 2.02 );
	type_map_add( type_map_, "Pd", "*", 1.63 );
	type_map_add( type_map_, "Pt", "*", 1.72 );
	type_map_add( type_map_, "Se", "*", 1.90 );
	type_map_add( type_map_, "Si", "*", 2.10 );
	type_map_add( type_map_, "Sn", "*", 2.17 );
	type_map_add( type_map_, "Te", "*", 2.06 );
	type_map_add( type_map_, "Tl", "*", 1.96 );
	type_map_add( type_map_, "U ", "*", 1.86 );
	type_map_add( type_map_, "Xe", "*", 2.16 );
	type_map_add( type_map_, "Zn", "*", 1.39 );

	// std rosetta stuff
	type_map_add( type_map_, "C"   ,"*", 2.0000 );
	type_map_add( type_map_, "N"   ,"*", 1.7500 );
	type_map_add( type_map_, "O"   ,"*", 1.5500 );
	type_map_add( type_map_, "S"   ,"*", 1.9000 );
	type_map_add( type_map_, "N"   ,"*", 1.7500 );
	type_map_add( type_map_, "C"   ,"*", 2.0000 );
	type_map_add( type_map_, "O"   ,"*", 1.5500 );
	type_map_add( type_map_, "P"   ,"*", 2.1500 );
	type_map_add( type_map_, "F"   ,"*", 1.7100 );
	type_map_add( type_map_, "CL"  ,"*", 2.0700 );
	type_map_add( type_map_, "BR"  ,"*", 2.2200 );
	type_map_add( type_map_, "I"   ,"*", 2.3600 );
	type_map_add( type_map_, "ZN"  ,"*", 1.0900 );
	type_map_add( type_map_, "FE"  ,"*", 0.7800 );
	type_map_add( type_map_, "FE"  ,"*", 0.6500 );
	type_map_add( type_map_, "MG"  ,"*", 1.1850 );
	type_map_add( type_map_, "CA"  ,"*", 1.3670 );
	type_map_add( type_map_, "NA"  ,"*", 1.3638 );
	type_map_add( type_map_, "K"   ,"*", 1.7638 );
	type_map_add( type_map_, "CNH2","*", 2.0000 );
	type_map_add( type_map_, "COO" ,"*", 2.0000 );
	type_map_add( type_map_, "CH1" ,"*", 2.0000 );
	type_map_add( type_map_, "CH2" ,"*", 2.0000 );
	type_map_add( type_map_, "CH3" ,"*", 2.0000 );
	type_map_add( type_map_, "aroC","*", 2.0000 );
	type_map_add( type_map_, "Ntrp","*", 1.7500 );
	type_map_add( type_map_, "Nhis","*", 1.7500 );
	type_map_add( type_map_, "NH2O","*", 1.7500 );
	type_map_add( type_map_, "Nlys","*", 1.7500 );
	type_map_add( type_map_, "Narg","*", 1.7500 );
	type_map_add( type_map_, "Npro","*", 1.7500 );
	type_map_add( type_map_, "OH"  ,"*", 1.5500 );
	type_map_add( type_map_, "ONH2","*", 1.5500 );
	type_map_add( type_map_, "OOC" ,"*", 1.5500 );
	type_map_add( type_map_, "S"   ,"*", 1.9000 );
	type_map_add( type_map_, "Nbb" ,"*", 1.7500 );
	type_map_add( type_map_, "CAbb","*", 2.0000 );
	type_map_add( type_map_, "CObb","*", 2.0000 );
	type_map_add( type_map_, "OCbb","*", 1.5500 );
	type_map_add( type_map_, "Phos","*", 2.1500 );
	type_map_add( type_map_, "F"   ,"*", 1.7100 );
	type_map_add( type_map_, "Cl"  ,"*", 2.0700 );
	type_map_add( type_map_, "Br"  ,"*", 2.2200 );
	type_map_add( type_map_, "I"   ,"*", 2.3600 );
	type_map_add( type_map_, "Zn","*", 1.0900 );
	type_map_add( type_map_, "Fe","*", 0.7800 );
	type_map_add( type_map_, "Mg","*", 1.1850 );
	type_map_add( type_map_, "Ca","*", 1.3670 );
	type_map_add( type_map_, "Na","*", 1.3638 );
	type_map_add( type_map_, "K" ,"*", 1.7638 );
	type_map_add( type_map_, "Cu" ,"*", 1.4 );

	if ( include_water ) {
		type_map_add( type_map_, "O" ,"HOH",  1.4 );
		type_map_add( type_map_, "*" ,"HOH",  1.4 );
		type_map_add( type_map_, "WA", "ADW", 1.4 ); // water?
		type_map_add( type_map_, "WB", "ADW", 1.4 ); // water?
		type_map_add( type_map_, "W", "*",    1.4 ); // water?
	} else {
		// no water!
		type_map_add( type_map_, "O" ,"HOH",  0 );
		type_map_add( type_map_, "*" ,"HOH",  0 );
		type_map_add( type_map_, "WA", "ADW", 0 ); // water?
		type_map_add( type_map_, "WB", "ADW", 0 ); // water?
		type_map_add( type_map_, "W", "*",    0 ); // water?
	}

	// virtual residues
	type_map_add( type_map_, "*" ,"SCK", 0 );
	type_map_add( type_map_, "VIRT" ,"*", 0 );
	type_map_add( type_map_, "*"    ,"X", 0 );

	// hydrogen
	type_map_add( type_map_, "H"   ,"*", HYDROGEN_RADIUS );
	type_map_add( type_map_, "Hpol","*", HYDROGEN_RADIUS );
	type_map_add( type_map_, "Hapo","*", HYDROGEN_RADIUS );
	type_map_add( type_map_, "Haro","*", HYDROGEN_RADIUS );
	type_map_add( type_map_, "HNbb","*", HYDROGEN_RADIUS );
}

PackstatReal
AtomRadiusMap::get_radius(
	std::string type,
	std::string res
) const {
	type = core::scoring::packstat::ToUpper(type);
	res  = core::scoring::packstat::ToUpper(res );

	if (      type_map_.count( P(type,res) ) ) {
		return type_map_.find(  P(type,res) )->second;
	}
	if (      type_map_.count( P(type,"*") ) ) {
		return type_map_.find(  P(type,"*") )->second;
	}
	if (      type_map_.count( P("*" ,res) ) ) {
		return type_map_.find(  P("*" ,res) )->second;
	}

	if ( 'H' == type[1] && '0' <= type[0] && type[0] <= '9' ) return HYDROGEN_RADIUS;

	if ( type.size() > 1 ) {
		type = type.substr(0,type.size()-1);
		return get_radius(type,res);
	}

	return -1;
}

} // namespace packstat
} // namespace scoring
} // namespace core
