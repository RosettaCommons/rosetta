// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/pilot/rmoretti/bootstrap_gasteiger_types.cc
///
/// @brief  BCL uses internalized function calls to initialize various types.
/// This application is here to convert those to the database file format of Rosetta.
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>
#include <core/chemical/Element.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/gasteiger/util.hh>

#include <core/types.hh>

#include <devel/init.hh>

#include <utility/vector1.hh>
#include <numeric/util.hh>
#include <utility/excn/EXCN_Base.hh>

#include <fstream>

utility::vector1< core::chemical::ElementCOP >
initialize_elements () {
	using namespace core::chemical;

	utility::vector1< core::chemical::ElementCOP > out;

    // 1st period                     atomic number  period  grOUP  Symbol    Element Name        mass            adiu  vdw radius  melting po                              sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d  f       pymol 3
 	out.push_back( new Element(              1,      1,     1, "H"     , "Hydrogen"     ,    1.01,           0.32,       1.20,         ElectronConfiguration(                             1,                   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Hydrogen
	out.push_back( new Element(              2,      1,     8, "He"    , "Helium"       ,       4,           0.93,       1.22,         ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Helium
	// 2st period                     atomic number  period  grOUp  Symbol    Element Name        mass            adiu  vdw radius  melting po                              sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d  f
	out.push_back( new Element(              3,      2,     1, "Li"    , "Lithium"      ,    6.94,           1.23,       1.52,         ElectronConfiguration(                             1,                   1,   2,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Lithium
	out.push_back( new Element(              4,      2,     2, "Be"    , "Beryllium"    ,    9.01,            0.9,       1.70,         ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Beryllium
	out.push_back( new Element(              5,      2,     3, "B"     , "Boron"        ,   10.81,           0.82,       2.08,         ElectronConfiguration(                             3,                   3,   2,   0,   0,   0,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Boron
	out.push_back( new Element(              6,      2,     4, "C"     , "Carbon"       ,   12.01,           0.77,       1.85,         ElectronConfiguration(                             4,                   4,   2,   0,   0,   0,   2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Carbon
	out.push_back( new Element(              7,      2,     5, "N"     , "Nitrogen"     ,   14.01,           0.75,       1.54,         ElectronConfiguration(                             5,                   5,   2,   0,   0,   0,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Nitrogen
	out.push_back( new Element(              8,      2,     6, "O"     , "Oxygen"       ,      16,           0.73,       1.40,         ElectronConfiguration(                             6,                   6,   2,   0,   0,   0,   2,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Oxygen
	out.push_back( new Element(              9,      2,     7, "F"     , "Fluorine"     ,      19,           0.72,       1.35,         ElectronConfiguration(                             7,                   7,   2,   0,   0,   0,   2,   5,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Fluorine
	out.push_back( new Element(             10,      2,     8, "Ne"    , "Neon"         ,   20.18,           0.71,       1.60,         ElectronConfiguration(                             8,                   8,   2,   0,   0,   0,   2,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Neon
	// 3st period                     atomic number  period  grOUp  Symbol    Element Name        mass            adiu  vdw radius  melting po                              sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d  f
	out.push_back( new Element(             11,      3,     1, "Na"    , "Sodium"       ,   22.99,           1.54,       2.31,         ElectronConfiguration(                             1,                   1,   2,   0,   0,   0,   2,   6,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Sodium
	out.push_back( new Element(             12,      3,     2, "Mg"    , "Magnesium"    ,   24.31,           1.36,       1.73,         ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Magnesium
	out.push_back( new Element(             13,      3,     3, "Al"    , "Aluminum"     ,   26.98,           1.18,       2.05,         ElectronConfiguration(                             3,                   3,   2,   0,   0,   0,   2,   6,   0,   0,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Aluminum
	out.push_back( new Element(             14,      3,     4, "Si"    , "Silicon"      ,   28.09,           1.11,       2.00,         ElectronConfiguration(                             4,                   4,   2,   0,   0,   0,   2,   6,   0,   0,   2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Silicon
	out.push_back( new Element(             15,      3,     5, "P"     , "Phosphorus"   ,   30.97,           1.06,       1.90,         ElectronConfiguration(                             5,                   5,   2,   0,   0,   0,   2,   6,   0,   0,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Phosphorus
	out.push_back( new Element(             16,      3,     6, "S"     , "Sulfur"       ,   32.07,           1.02,       1.85,         ElectronConfiguration(                             6,                   6,   2,   0,   0,   0,   2,   6,   0,   0,   2,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Sulfur
	out.push_back( new Element(             17,      3,     7, "Cl"    , "Chlorine"     ,   35.45,           0.99,       1.81,         ElectronConfiguration(                             7,                   7,   2,   0,   0,   0,   2,   6,   0,   0,   2,   5,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Chlorine
	out.push_back( new Element(             18,      3,     8, "Ar"    , "Argon"        ,   39.95,           0.98,       1.91,         ElectronConfiguration(                             8,                   8,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Argon
	// 4st period                     atomic number  period  grOUp  Symbol    Element Name        mass            adiu  vdw radius  melting po                              sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d  f
	out.push_back( new Element(             19,      4,     1, "K"     , "Potassium"    ,    39.1,           2.03,       2.31,         ElectronConfiguration(                             1,                   1,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Potassium
	out.push_back( new Element(             20,      4,     2, "Ca"    , "Calcium"      ,   40.08,           1.74,       1.97,         ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Calcium
	out.push_back( new Element(             21,      4,     0, "Sc"    , "Scandium"     ,   44.96,           1.44,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   3,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   1,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Scandium
	out.push_back( new Element(             22,      4,     0, "Ti"    , "Titanium"     ,   47.88,           1.32,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   4,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   2,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Titanium
	out.push_back( new Element(             23,      4,     0, "V"     , "Vanadium"     ,   50.94,           1.22,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   5,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   3,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Vanadium
	out.push_back( new Element(             24,      4,     0, "Cr"    , "Chromium"     ,      52,           1.18,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   6,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   4,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Chromium
	out.push_back( new Element(             25,      4,     0, "Mn"    , "Manganese"    ,   54.94,           1.17,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   7,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   5,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Manganese
	out.push_back( new Element(             26,      4,     0, "Fe"    , "Iron"         ,   55.85,           1.17,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   8,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   6,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Iron
	out.push_back( new Element(             27,      4,     0, "Co"    , "Cobalt"       ,   58.93,           1.16,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   9,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   7,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Cobalt
	out.push_back( new Element(             28,      4,     0, "Ni"    , "Nickel"       ,   58.69,           1.15,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  10,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   8,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Nickel
	out.push_back( new Element(             29,      4,     0, "Cu"    , "Copper"       ,   63.55,           1.17,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  11,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Copper
	out.push_back( new Element(             30,      4,     0, "Zn"    , "Zinc"         ,   65.39,           1.25,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  12,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Zinc
	out.push_back( new Element(             31,      4,     3, "Ga"    , "Gallium"      ,   69.72,           1.26,       1.70,         ElectronConfiguration(                             3,                  13,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Gallium
	out.push_back( new Element(             32,      4,     4, "Ge"    , "Germanium"    ,   72.61,           1.22,       1.70,         ElectronConfiguration(                             4,                  14,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Germanium
	out.push_back( new Element(             33,      4,     5, "As"    , "Arsenic"      ,   74.92,            1.2,       2.00,         ElectronConfiguration(                             5,                  15,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Arsenic
	out.push_back( new Element(             34,      4,     6, "Se"    , "Selenium"     ,   78.96,           1.16,       2.00,         ElectronConfiguration(                             6,                  16,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Selenium
	out.push_back( new Element(             35,      4,     7, "Br"    , "Bromine"      ,    79.9,           1.14,       2.10,         ElectronConfiguration(                             7,                  17,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   5,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Bromine
	out.push_back( new Element(             36,      4,     8, "Kr"    , "Krypton"      ,    83.8,           1.12,       2.10,         ElectronConfiguration(                             8,                  18,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Krypton
	// 5st period                     atomic number  period  grOUp  Symbol    Element Name        mass            adiu  vdw radius  melting po                              sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d  f
	out.push_back( new Element(             37,      5,     1, "Rb"    , "Rubidium"     ,   85.47,           2.16,       1.70,         ElectronConfiguration(                             1,                   1,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Rubidium
	out.push_back( new Element(             38,      5,     2, "Sr"    , "Strontium"    ,   87.62,           1.91,       1.70,         ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Strontium
	out.push_back( new Element(             39,      5,     0, "Y"     , "Yttrium"      ,   88.91,           1.62,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   3,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   1,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Yttrium
	out.push_back( new Element(             40,      5,     0, "Zr"    , "Zirconium"    ,   91.22,           1.45,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   4,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   2,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Zirconium
	out.push_back( new Element(             41,      5,     0, "Nb"    , "Niobium"      ,   92.91,           1.34,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   5,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   4,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Niobium
	out.push_back( new Element(             42,      5,     0, "Mo"    , "Molybdenum"   ,   95.94,            1.3,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   6,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   5,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Molybdenum
	out.push_back( new Element(             43,      5,     0, "Tc"    , "Technetium"   ,      98,           1.27,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   7,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   6,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Technetium
	out.push_back( new Element(             44,      5,     0, "Ru"    , "Ruthenium"    ,  101.07,           1.25,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   8,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   7,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Ruthenium
	out.push_back( new Element(             45,      5,     0, "Rh"    , "Rhodium"      ,  102.91,           1.25,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   9,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   8,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Rhodium
	out.push_back( new Element(             46,      5,     0, "Pd"    , "Palladium"    ,  106.42,           1.28,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  10,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Palladium
	out.push_back( new Element(             47,      5,     0, "Ag"    , "Silver"       ,  107.87,           1.34,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  11,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Silver
	out.push_back( new Element(             48,      5,     0, "Cd"    , "Cadmium"      ,  112.41,           1.48,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  12,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Cadmium
	out.push_back( new Element(             49,      5,     3, "In"    , "Indium"       ,  114.82,           1.44,       1.70,         ElectronConfiguration(                             3,                  13,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Indium
	out.push_back( new Element(             50,      5,     4, "Sn"    , "Tin"          ,  118.71,           1.41,       1.70,         ElectronConfiguration(                             4,                  14,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Tin
	out.push_back( new Element(             51,      5,     5, "Sb"    , "Antimony"     ,  121.76,            1.4,       2.20,         ElectronConfiguration(                             5,                  15,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Antimony
	out.push_back( new Element(             52,      5,     6, "Te"    , "Tellurium"    ,   127.6,           1.36,       2.20,         ElectronConfiguration(                             6,                  16,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Tellurium
	out.push_back( new Element(             53,      5,     7, "I"     , "Iodine"       ,   126.9,           1.33,       2.15,         ElectronConfiguration(                             7,                  17,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   5,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Iodine
	out.push_back( new Element(             54,      5,     8, "Xe"    , "Xenon"        ,  131.29,           1.31,       2.16,         ElectronConfiguration(                             8,                  18,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 )));// Xenon
	// 6st period                     atomic number  period  grOUp  Symbol    Element Name        mass            adiu  vdw radius  melting po                              sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d  f
	out.push_back( new Element(             55,      6,     1, "Cs"    , "Cesium"       ,  132.91,           2.35,       1.70,         ElectronConfiguration(                             1,                   1,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   6,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0 )));// Cesium
	out.push_back( new Element(             56,      6,     2, "Ba"    , "Barium"       ,  137.33,           1.98,       1.70,         ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Barium
	out.push_back( new Element(             57,      6,     0, "La"    , "Lanthanum"    ,  138.91,           1.69,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   3,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   6,   1,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Lanthanum
	out.push_back( new Element(             58,      6,     0, "Ce"    , "Cerium"       ,  140.12,           1.65,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   4,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   1,   2,   6,   1,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Cerium
	out.push_back( new Element(             59,      6,     0, "Pr"    , "Praseodymium" ,  140.91,           1.65,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   5,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   3,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Praseodymium
	out.push_back( new Element(             60,      6,     0, "Nd"    , "Neodymium"    ,  144.24,           1.64,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   6,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   4,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Neodymium
	out.push_back( new Element(             61,      6,     0, "Pm"    , "Promethium"   ,     145,           1.63,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   7,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   5,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Promethium
	out.push_back( new Element(             62,      6,     0, "Sm"    , "Samarium"     ,  150.36,           1.62,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   8,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   6,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Samarium
	out.push_back( new Element(             63,      6,     0, "Eu"    , "Europium"     ,  151.97,           1.85,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   9,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   7,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Europium
	out.push_back( new Element(             64,      6,     0, "Gd"    , "Gadolinium"   ,  157.25,           1.61,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  10,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   7,   2,   6,   1,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Gadolinium
	out.push_back( new Element(             65,      6,     0, "Tb"    , "Terbium"      ,  158.93,           1.59,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  11,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   9,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Terbium
	out.push_back( new Element(             66,      6,     0, "Dy"    , "Dysprosium"   ,   162.5,           1.59,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  12,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  10,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Dysprosium
	out.push_back( new Element(             67,      6,     0, "Ho"    , "Holmium"      ,  164.93,           1.58,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  13,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  11,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Holmium
	out.push_back( new Element(             68,      6,     0, "Er"    , "Erbium"       ,  167.26,           1.57,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  14,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  12,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Erbium
	out.push_back( new Element(             69,      6,     0, "Tm"    , "Thulium"      ,  168.93,           1.56,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  15,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  13,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Thulium
	out.push_back( new Element(             70,      6,     0, "Yb"    , "Ytterbium"    ,  173.04,           1.74,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  16,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Ytterbium
	out.push_back( new Element(             71,      6,     0, "Lu"    , "Lutetium"     ,  174.97,           1.56,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  17,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   1,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Lutetium
	out.push_back( new Element(             72,      6,     0, "Hf"    , "Hafnium"      ,  178.49,           1.44,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  18,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   2,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Hafnium
	out.push_back( new Element(             73,      6,     0, "Ta"    , "Tantalum"     ,  180.95,           1.34,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  19,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   3,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Tantalum
	out.push_back( new Element(             74,      6,     0, "W"     , "Tungsten"     ,  183.85,            1.3,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  20,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   4,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Tungsten
	out.push_back( new Element(             75,      6,     0, "Re"    , "Rhenium"      ,  186.21,           1.28,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  21,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   5,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Rhenium
	out.push_back( new Element(             76,      6,     0, "Os"    , "Osmium"       ,   190.2,           1.26,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  22,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   6,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Osmium
	out.push_back( new Element(             77,      6,     0, "Ir"    , "Iridium"      ,  192.22,           1.27,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  23,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   7,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Iridium
	out.push_back( new Element(             78,      6,     0, "Pt"    , "Platinum"     ,  195.08,            1.3,       1.72,         ElectronConfiguration( numeric::get_undefined_size(),                  24,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   9,   0,   1,   0,   0,   0,   0,   0,   0,   0 )));// Platinum
	out.push_back( new Element(             79,      6,     0, "Au"    , "Gold"         ,  196.97,           1.34,       1.66,         ElectronConfiguration( numeric::get_undefined_size(),                  25,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   1,   0,   0,   0,   0,   0,   0,   0 )));// Gold
	out.push_back( new Element(             80,      6,     0, "Hg"    , "Mercury"      ,  200.59,           1.49,       1.55,         ElectronConfiguration( numeric::get_undefined_size(),                  26,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   0,   0,   0,   0,   0,   0,   0 )));// Mercury
	out.push_back( new Element(             81,      6,     3, "Tl"    , "Thallium"     ,  204.38,           1.48,       1.70,         ElectronConfiguration(                             3,                  27,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   1,   0,   0,   0,   0,   0,   0 )));// Thallium
	out.push_back( new Element(             82,      6,     4, "Pb"    , "Lead"         ,   207.2,           1.47,       1.70,         ElectronConfiguration(                             4,                  28,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   2,   0,   0,   0,   0,   0,   0 )));// Lead
	out.push_back( new Element(             83,      6,     5, "Bi"    , "Bismuth"      ,  208.98,           1.46,       1.70,         ElectronConfiguration(                             5,                  29,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   3,   0,   0,   0,   0,   0,   0 )));// Bismuth
	out.push_back( new Element(             84,      6,     6, "Po"    , "Polonium"     ,     209,           1.46,       1.70,         ElectronConfiguration(                             6,                  30,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   4,   0,   0,   0,   0,   0,   0 )));// Polonium
	out.push_back( new Element(             85,      6,     7, "At"    , "Astatine"     ,     210,           1.45,       1.70,         ElectronConfiguration(                             7,                  31,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   5,   0,   0,   0,   0,   0,   0 )));// Astatine
	out.push_back( new Element(             86,      6,     8, "Rn"    , "Radon"        ,     222,              0,       1.70,         ElectronConfiguration(                             8,                  32,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   6,   0,   0,   0,   0,   0,   0 )));// Radon
	// 7st period                     atomic number  period  grOUp  Symbol    Element Name        mass            adiu  vdw radius  melting po                              sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d  f
	out.push_back( new Element(             87,      7,     1, "Fr"    , "Francium"     ,     223,              0,       1.70,         ElectronConfiguration(                             1,                   1,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   6,   0,   0,   1,   0,   0,   0 )));// Francium
	out.push_back( new Element(             88,      7,     2, "Ra"    , "Radium"       ,  226.03,              0,       1.70,         ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   6,   0,   0,   2,   0,   0,   0 )));// Radium
	out.push_back( new Element(             89,      7,     0, "Ac"    , "Actinium"     ,     227,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   3,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   1,   6,   1,   0,   2,   0,   0,   0 )));// Actinium
	out.push_back( new Element(             90,      7,     0, "Th"    , "Thorium"      ,  232.04,           1.65,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   4,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   6,   2,   0,   2,   0,   0,   0 )));// Thorium
	out.push_back( new Element(             91,      7,     0, "Pa"    , "Protactinium" ,  213.04,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   5,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   2,   1,   6,   1,   0,   2,   0,   0,   0 )));// Protactinium
	out.push_back( new Element(             92,      7,     0, "U"     , "Uranium"      ,  238.03,           1.42,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   6,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   3,   1,   6,   1,   0,   2,   0,   0,   0 )));// Uranium
	out.push_back( new Element(             93,      7,     0, "Np"    , "Neptunium"    ,  237.05,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   7,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   4,   1,   6,   1,   0,   2,   0,   0,   0 )));// Neptunium
	out.push_back( new Element(             94,      7,     0, "Pu"    , "Plutonium"    ,     244,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   8,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   6,   2,   6,   0,   0,   2,   0,   0,   0 )));// Plutonium
	out.push_back( new Element(             95,      7,     0, "Am"    , "Americium"    ,     243,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                   9,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   7,   2,   6,   0,   0,   2,   0,   0,   0 )));// Americium
	out.push_back( new Element(             96,      7,     0, "Cm"    , "Curium"       ,     247,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  10,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   7,   1,   6,   1,   0,   2,   0,   0,   0 )));// Curium
	out.push_back( new Element(             97,      7,     0, "Bk"    , "Berkelium"    ,     247,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  11,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   9,   2,   6,   0,   0,   2,   0,   0,   0 )));// Berkelium
	out.push_back( new Element(             98,      7,     0, "Cf"    , "Californium"  ,     251,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  12,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  10,   2,   6,   0,   0,   2,   0,   0,   0 )));// Californium
	out.push_back( new Element(             99,      7,     0, "Es"    , "Einsteinium"  ,     252,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  13,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  11,   2,   6,   0,   0,   2,   0,   0,   0 )));// Einsteinium
	out.push_back( new Element(            100,      7,     0, "Fm"    , "Fermium"      ,     257,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  14,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  12,   2,   6,   0,   0,   2,   0,   0,   0 )));// Fermium
	out.push_back( new Element(            101,      7,     0, "Md"    , "Mendelevium"  ,     258,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  15,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  13,   2,   6,   0,   0,   2,   0,   0,   0 )));// Mendelevium
	out.push_back( new Element(            102,      7,     0, "No"    , "Nobelium"     ,     259,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  16,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   0,   0,   2,   0,   0,   0 )));// Nobelium
	out.push_back( new Element(            103,      7,     0, "Lr"    , "Lawrencium"   ,     260,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  17,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   1,   0,   2,   0,   0,   0 )));// Lawrencium
	out.push_back( new Element(            104,      7,     0, "Rf"    , "Rutherfordium",     261,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  18,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   2,   0,   2,   0,   0,   0 )));// Rutherfordiu
	out.push_back( new Element(            105,      7,     0, "Db"    , "Dubnium"      ,     262,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  19,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   3,   0,   2,   0,   0,   0 )));// Dubnium
	out.push_back( new Element(            106,      7,     0, "Sg"    , "Seaborgium"   ,     263,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  20,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   4,   0,   2,   0,   0,   0 )));// Seaborgium
	out.push_back( new Element(            107,      7,     0, "Bh"    , "Bohrium"      ,     262,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  21,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   5,   0,   2,   0,   0,   0 )));// Bohrium
	out.push_back( new Element(            108,      7,     0, "Hs"    , "Hassium"      ,     265,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  22,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   6,   0,   2,   0,   0,   0 )));// Hassium
	out.push_back( new Element(            109,      7,     0, "Mt"    , "Meitnerium"   ,     266,              0,       1.70,         ElectronConfiguration( numeric::get_undefined_size(),                  23,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   7,   0,   2,   0,   0,   0 )));// Meitnerium

	return out;

/* Copy - paste from BCL
    // 1st period              element name                       atomic number  period  group  Symbol    Element Name        mass                  gyromagnetic ratio  covalent radiu  vdw radius  melting poin  boiling poi                     electro negativity  ionization poten                                  sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d   7f       pymol rgb r       pymol rgb g        pymol rgb b                               CromerMann A1    CromerMann A2    CromerMann A3    CromerMann A4    CromerMann B1    CromerMann B2    CromerMann B3    CromerMann B4     CromerMann C   Displaced Volume A^3
    e_Hydrogen(      AddEnum(  "Hydrogen"     , Element(              1,      1,     1, "H"     , "Hydrogen"     ,    1.01,                          267510000,           0.32,       1.20,        14.03,       20.27,                                   2.2,             13.6, ElectronConfiguration(                             1,                   1,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.900000000f,      0.900000000f,     0.900000000f, ElementStructureFactor(       0.493002,        0.322912,        0.140191,        0.040810,       10.510900,       26.125700,        3.142360,       57.799700,        0.003038,  5.1495025)))),// Hydrogen
    e_Helium(        AddEnum(  "Helium"       , Element(              2,      1,     8, "He"    , "Helium"       ,       4,                          203780000,           0.93,       1.22,         0.95,        4.22,         util::GetUndefined< double>(),            24.59, ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.850980392f,      1.000000000f,     1.000000000f, ElementStructureFactor(       0.873400,        0.630900,        0.311200,        0.178000,        9.103700,        3.356800,       22.927600,        0.982100,        0.064000            )))),// Helium
    // 2nd period              element name                       atomic number  period  group  Symbol    Element Name        mass                  gyromagnetic ratio  covalent radiu  vdw radius  melting poin  boiling poi                     electro negativity  ionization poten                                  sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d   7f
    e_Lithium(       AddEnum(  "Lithium"      , Elem#include <devel/init.hh>entTypeData(              3,      2,     1, "Li"    , "Lithium"      ,    6.94,                          103960000,           1.23,       1.52,        453.7,        1615,                                     1,             5.39, ElectronConfiguration(                             1,                   1,   2,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.800000000f,      0.501960784f,     1.000000000f, ElementStructureFactor(       1.128200,        0.750800,        0.617500,        0.465300,        3.954600,        1.052400,       85.390500,      168.261000,        0.037700            )))),// Lithium
    e_Beryllium(     AddEnum(  "Beryllium"    , Element(              4,      2,     2, "Be"    , "Beryllium"    ,    9.01,                           37590000,            0.9,       1.70,         1551,        3243,                                   1.6,             9.32, ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.760784314f,      1.000000000f,     0.000000000f, ElementStructureFactor(       1.591900,        1.127800,        0.539100,        0.702900,       43.642700,        1.862300,      103.483000,        0.542000,        0.038500            )))),// Beryllium
    e_Boron(         AddEnum(  "Boron"        , Element(              5,      2,     3, "B"     , "Boron"        ,   10.81,                           85796000,           0.82,       2.08,         2573,        4275,                                     2,              8.3, ElectronConfiguration(                             3,                   3,   2,   0,   0,   0,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    1.000000000f,      0.709803922f,     0.709803922f, ElementStructureFactor(       2.054500,        1.332600,        1.097900,        0.706800,       23.218500,        1.021000,       60.349800,        0.140300,       -0.193200            )))),// Boron
    e_Carbon(        AddEnum(  "Carbon"       , Element(              6,      2,     4, "C"     , "Carbon"       ,   12.01,                           67263000,           0.77,       1.85,         3773,        5100,                                   2.6,            11.26, ElectronConfiguration(                             4,                   4,   2,   0,   0,   0,   2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.200000000f,      1.000000000f,     0.200000000f, ElementStructureFactor(       2.310000,        1.020000,        1.588600,        0.865000,       20.843900,       10.207500,        0.568700,       51.651200,        0.215600, 16.4451800)))),// Carbon
    e_Nitrogen(      AddEnum(  "Nitrogen"     , Element(              7,      2,     5, "N"     , "Nitrogen"     ,   14.01,                          -27117000,           0.75,       1.54,        63.14,       77.35,                                     3,            14.53, ElectronConfiguration(                             5,                   5,   2,   0,   0,   0,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.200000000f,      0.200000000f,     1.000000000f, ElementStructureFactor(      12.212600,        3.132200,        2.012500,        1.166300,        0.005700,        9.893300,       28.997500,        0.582600,      -11.529000,  2.4916940)))),// Nitrogen
    e_Oxygen(        AddEnum(  "Oxygen"       , Element(              8,      2,     6, "O"     , "Oxygen"       ,      16,                           36264000,           0.73,       1.40,        50.35,       90.18,                                   3.4,            13.62, ElectronConfiguration(                             6,                   6,   2,   0,   0,   0,   2,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    1.000000000f,      0.300000000f,     0.300000000f, ElementStructureFactor(       3.048500,        2.286800,        1.546300,        0.867000,       13.277100,        5.701100,        0.323900,       32.908900,        0.250800,  9.1362130)))),// Oxygen
    e_Fluorine(      AddEnum(  "Fluorine"     , Element(              9,      2,     7, "F"     , "Fluorine"     ,      19,                          251710000,           0.72,       1.35,        53.48,       84.95,                                     4,            17.42, ElectronConfiguration(                             7,                   7,   2,   0,   0,   0,   2,   5,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.701960784f,      1.000000000f,     1.000000000f, ElementStructureFactor(       3.539200,        2.641200,        1.517000,        1.024300,       10.282500,        4.294400,        0.261500,       26.147600,        0.277600            )))),// Fluorine
    e_Neon(          AddEnum(  "Neon"         , Element(             10,      2,     8, "Ne"    , "Neon"         ,   20.18,                           21117000,           0.71,       1.60,        24.55,        27.1,         util::GetUndefined< double>(),            21.56, ElectronConfiguration(                             8,                   8,   2,   0,   0,   0,   2,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.701960784f,      0.890196078f,     0.960784314f, ElementStructureFactor(       3.955300,        3.112500,        1.454600,        1.125100,        8.404200,        3.426200,        0.230600,       21.718400,        0.351500            )))),// Neon
    // 3rd period              element name                       atomic number  period  group  Symbol    Element Name        mass                  gyromagnetic ratio  covalent radiu  vdw radius  melting poin  boiling poi                     electro negativity  ionization poten                                  sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d   7f
    Sodium(        AddEnum(  "Sodium"       , Element(             11,      3,     1, "Na"    , "Sodium"       ,   22.99,                           70761000,           1.54,       2.31,          371,        1156,                                   0.9,             5.14, ElectronConfiguration(                             1,                   1,   2,   0,   0,   0,   2,   6,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.670588235f,      0.360784314f,     0.949019608f, ElementStructureFactor(       4.762600,        3.173600,        1.267400,        1.112800,        3.285000,        8.842200,        0.313600,      129.424000,        0.676000            )))),// Sodium
    e_Magnesium(     AddEnum(  "Magnesium"    , Element(             12,      3,     2, "Mg"    , "Magnesium"    ,   24.31,                           16377000,           1.36,       1.73,          922,        1363,                                   1.3,             7.65, ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.541176471f,      1.000000000f,     0.000000000f, ElementStructureFactor(       5.420400,        2.173500,        1.226900,        2.307300,        2.827500,       79.261100,        0.380800,        7.193700,        0.858400            )))),// Magnesium
    e_Aluminum(      AddEnum(  "Aluminum"     , Element(             13,      3,     3, "Al"    , "Aluminum"     ,   26.98,                           69705000,           1.18,       2.05,       933.25,        2740,                                   1.6,             5.99, ElectronConfiguration(                             3,                   3,   2,   0,   0,   0,   2,   6,   0,   0,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.749019608f,      0.650980392f,     0.650980392f, ElementStructureFactor(       6.420200,        1.900200,        1.593600,        1.964600,        3.038700,        0.742600,       31.547200,       85.088600,        1.115100            )))),// Aluminum
    Silicon(       AddEnum(  "Silicon"      , Element(             14,      3,     4, "Si"    , "Silicon"      ,   28.09,                           53146000,           1.11,       2.00,         1683,        2628,                                   1.9,             8.15, ElectronConfiguration(                             4,                   4,   2,   0,   0,   0,   2,   6,   0,   0,   2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.941176471f,      0.784313725f,     0.627450980f, ElementStructureFactor(       6.291500,        3.035300,        1.989100,        1.541000,        2.438600,       32.333700,        0.678500,       81.693700,        1.140700            )))),// Silicon
    e_Phosphorus(    AddEnum(  "Phosphorus"   , Element(             15,      3,     5, "P"     , "Phosphorus"   ,   30.97,                          108290000,           1.06,       1.90,        317.3,         553,                                   2.2,            10.49, ElectronConfiguration(                             5,                   5,   2,   0,   0,   0,   2,   6,   0,   0,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    1.000000000f,      0.501960784f,     0.000000000f, ElementStructureFactor(       6.434500,        4.179100,        1.780000,        1.490800,        1.906700,       27.157000,        0.526000,       68.164500,        1.114900            )))),// Phosphorus
    Sulfur(        AddEnum(  "Sulfur"       , Element(             16,      3,     6, "S"     , "Sulfur"       ,   32.07,                           20534000,           1.02,       1.85,       388.36,      717.75,                                   2.6,            10.36, ElectronConfiguration(                             6,                   6,   2,   0,   0,   0,   2,   6,   0,   0,   2,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.900000000f,      0.775000000f,     0.250000000f, ElementStructureFactor(       6.905300,        5.203400,        1.437900,        1.586300,        1.467900,       22.215100,        0.253600,       56.172000,        0.866900,  25.747510)))),// Sulfur
    e_Chlorine(      AddEnum(  "Chlorine"     , Element(             17,      3,     7, "Cl"    , "Chlorine"     ,   35.45,                           26211000,           0.99,       1.81,       172.16,       239.1,                                   3.2,            12.97, ElectronConfiguration(                             7,                   7,   2,   0,   0,   0,   2,   6,   0,   0,   2,   5,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.121568627f,      0.941176471f,     0.121568627f, ElementStructureFactor(      11.460400,        7.196400,        6.255600,        1.645500,        0.010400,        1.166200,       18.519400,       47.778400,       -9.557400            )))),// Chlorine
    e_Argon(         AddEnum(  "Argon"        , Element(             18,      3,     8, "Ar"    , "Argon"        ,   39.95,      util::GetUndefined< double>(),           0.98,       1.91,        83.81,        87.3,         util::GetUndefined< double>(),            15.76, ElectronConfiguration(                             8,                   8,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.501960784f,      0.819607843f,     0.890196078f, ElementStructureFactor(       7.484500,        6.772300,        0.653900,        1.644200,        0.907200,       14.840700,       43.898300,       33.392900,        1.444500            )))),// Argon
    // 4th period              element name                       atomic number  period  group  Symbol    Element Name        mass                  gyromagnetic ratio  covalent radiu  vdw radius  melting poin  boiling poi                     electro negativity  ionization poten                                  sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d   7f
    e_Potassium(     AddEnum(  "Potassium"    , Element(             19,      4,     1, "K"     , "Potassium"    ,    39.1,                           12482000,           2.03,       2.31,       336.35,        1032,                                   0.8,             4.34, ElectronConfiguration(                             1,                   1,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.560784314f,      0.250980392f,     0.831372549f, ElementStructureFactor(       8.218600,        7.439800,        1.051900,        0.865900,       12.794900,        0.774800,      213.187000,       41.684100,        1.422800            )))),// Potassium
    e_Calcium(       AddEnum(  "Calcium"      , Element(             20,      4,     2, "Ca"    , "Calcium"      ,   40.08,                           18001000,           1.74,       1.97,         1112,        1757,                                     1,             6.11, ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.239215686f,      1.000000000f,     0.000000000f, ElementStructureFactor(       8.626600,        7.387300,        1.589900,        1.021100,       10.442100,        0.659900,       85.748400,      178.437000,        1.375100            )))),// Calcium
    Scandium(      AddEnum(  "Scandium"     , Element(             21,      4,     0, "Sc"    , "Scandium"     ,   44.96,                           64984000,           1.44,       1.70,         1812,        3104,                                   1.4,             6.54, ElectronConfiguration( util::GetUndefined< size_t>(),                   3,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   1,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.901960784f,      0.901960784f,     0.901960784f, ElementStructureFactor(       9.189000,        7.367900,        1.640900,        1.468000,        9.021300,        0.572900,      136.108000,       51.353100,        1.332900            )))),// Scandium
    e_Titanium(      AddEnum(  "Titanium"     , Element(             22,      4,     0, "Ti"    , "Titanium"     ,   47.88,                           15081000,           1.32,       1.70,         1933,        3560,                                   1.5,             6.82, ElectronConfiguration( util::GetUndefined< size_t>(),                   4,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   2,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.749019608f,      0.760784314f,     0.780392157f, ElementStructureFactor(       9.759500,        7.355800,        1.699100,        1.902100,        7.850800,        0.500000,       35.633800,      116.105000,        1.280700            )))),// Titanium
    e_Vanadium(      AddEnum(  "Vanadium"     , Element(             23,      4,     0, "V"     , "Vanadium"     ,   50.94,                           70363000,           1.22,       1.70,         2175,        3682,                                   1.6,             6.74, ElectronConfiguration( util::GetUndefined< size_t>(),                   5,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   3,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.650980392f,      0.650980392f,     0.670588235f, ElementStructureFactor(      10.297100,        7.351100,        2.070300,        2.057100,        6.865700,        0.438500,       26.893800,      102.478000,        1.219900            )))),// Vanadium
    e_Chromium(      AddEnum(  "Chromium"     , Element(             24,      4,     0, "Cr"    , "Chromium"     ,      52,                           15119000,           1.18,       1.70,         2130,        2945,                                   1.7,             6.77, ElectronConfiguration( util::GetUndefined< size_t>(),                   6,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   4,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.541176471f,      0.600000000f,     0.780392157f, ElementStructureFactor(      10.640600,        7.353700,        3.324000,        1.492200,        6.103800,        0.392000,       20.262600,       98.739900,        1.183200            )))),// Chromium
    e_Manganese(     AddEnum(  "Manganese"    , Element(             25,      4,     0, "Mn"    , "Manganese"    ,   54.94,                           66195000,           1.17,       1.70,         1517,        2235,                                   1.6,             7.44, ElectronConfiguration( util::GetUndefined< size_t>(),                   7,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   5,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.611764706f,      0.478431373f,     0.780392157f, ElementStructureFactor(      11.281900,        7.357300,        3.019300,        2.244100,        5.340900,        0.343200,       17.867400,       83.754300,        1.089600            )))),// Manganese
    e_Iron(          AddEnum(  "Iron"         , Element(             26,      4,     0, "Fe"    , "Iron"         ,   55.85,                            8661800,           1.17,       1.70,         1808,        3023,                                   1.8,             7.87, ElectronConfiguration( util::GetUndefined< size_t>(),                   8,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   6,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.878431373f,      0.400000000f,     0.200000000f, ElementStructureFactor(      11.769500,        7.357300,        3.522200,        2.304500,        4.761100,        0.307200,       15.353500,       76.880500,        1.036900            )))),// Iron
    e_Cobalt(        AddEnum(  "Cobalt"       , Element(             27,      4,     0, "Co"    , "Cobalt"       ,   58.93,                           63472000,           1.16,       1.70,         1768,        3143,                                   1.9,             7.86, ElectronConfiguration( util::GetUndefined< size_t>(),                   9,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   7,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.941176471f,      0.564705882f,     0.627450980f, ElementStructureFactor(      12.284100,        7.340900,        4.003400,        2.348800,        4.279100,        0.278400,       13.535900,       71.169200,        1.011800            )))),// Cobalt
    e_Nickel(        AddEnum(  "Nickel"       , Element(             28,      4,     0, "Ni"    , "Nickel"       ,   58.69,                           23905000,           1.15,       1.70,         1726,        3005,                                   1.9,             7.64, ElectronConfiguration( util::GetUndefined< size_t>(),                  10,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,   8,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.313725490f,      0.815686275f,     0.313725490f, ElementStructureFactor(      12.837600,        7.292000,        4.443800,        2.380000,        3.878500,        0.256500,       12.176300,       66.342100,        1.034100            )))),// Nickel
    e_Copper(        AddEnum(  "Copper"       , Element(             29,      4,     0, "Cu"    , "Copper"       ,   63.55,                           70965000,           1.17,       1.70,       1357.6,        2840,                                   1.9,             7.73, ElectronConfiguration( util::GetUndefined< size_t>(),                  11,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.784313725f,      0.501960784f,     0.200000000f, ElementStructureFactor(      13.338000,        7.167600,        5.615800,        1.673500,        3.582800,        0.247000,       11.396600,       64.812600,        1.191000            )))),// Copper
    e_Zinc(          AddEnum(  "Zinc"         , Element(             30,      4,     0, "Zn"    , "Zinc"         ,   65.39,                           16738000,           1.25,       1.70,       692.73,        1180,                                   1.7,             9.39, ElectronConfiguration( util::GetUndefined< size_t>(),                  12,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.490196078f,      0.501960784f,     0.690196078f, ElementStructureFactor(      14.074300,        7.031800,        5.165200,        2.410000,        3.265500,        0.233300,       10.316300,       58.709700,        1.304100            )))),// Zinc
    e_Gallium(       AddEnum(  "Gallium"      , Element(             31,      4,     3, "Ga"    , "Gallium"      ,   69.72,                           81583000,           1.26,       1.70,        302.9,        2676,                                   1.8,                6, ElectronConfiguration(                             3,                  13,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.760784314f,      0.560784314f,     0.560784314f, ElementStructureFactor(      15.235400,        6.700600,        4.359100,        2.962300,        3.066900,        0.241200,       10.780500,       61.413500,        1.718900            )))),// Gallium
    e_Germanium(     AddEnum(  "Germanium"    , Element(             32,      4,     4, "Ge"    , "Germanium"    ,   72.61,                            9331200,           1.22,       1.70,       1210.4,        3103,                                     2,              7.9, ElectronConfiguration(                             4,                  14,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.400000000f,      0.560784314f,     0.560784314f, ElementStructureFactor(      16.081600,        6.374700,        3.706800,        3.683000,        2.850900,        0.251600,       11.446800,       54.762500,        2.131300            )))),// Germanium
    e_Arsenic(       Ad
    dEnum(  "Arsenic"      , Element(             33,      4,     5, "As"    , "Arsenic"      ,   74.92,                           45806000,            1.2,       2.00,         1081,         876,                                   2.2,             9.81, ElectronConfiguration(                             5,                  15,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.741176471f,      0.501960784f,     0.890196078f, ElementStructureFactor(      16.672300,        6.070100,        3.431300,        4.277900,        2.634500,        0.264700,       12.947900,       47.797200,        2.531000            )))),// Arsenic
    Selenium(      AddEnum(  "Selenium"     , Element(             34,      4,     6, "Se"    , "Selenium"     ,   78.96,                           51020000,           1.16,       2.00,          494,         958,                                   2.6,             9.75, ElectronConfiguration(                             6,                  16,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    1.000000000f,      0.631372549f,     0.000000000f, ElementStructureFactor(      17.000600,        5.819600,        3.973100,        4.354300,        2.409800,        0.272600,       15.237200,       43.816300,        2.840900,  28.730910)))),// Selenium
    e_Bromine(       AddEnum(  "Bromine"      , Element(             35,      4,     7, "Br"    , "Bromine"      ,    79.9,                           72246000,           1.14,       2.10,        265.9,      332.25,                                     3,            11.81, ElectronConfiguration(                             7,                  17,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   5,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.650980392f,      0.160784314f,     0.160784314f, ElementStructureFactor(      17.178900,        5.235800,        5.637700,        3.985100,        2.172300,       16.579600,        0.260900,       41.432800,        2.955700            )))),// Bromine
    e_Krypton(       AddEnum(  "Krypton"      , Element(             36,      4,     8, "Kr"    , "Krypton"      ,    83.8,      util::GetUndefined< double>(),           1.12,       2.10,       115.78,       119.8,         util::GetUndefined< double>(),               14, ElectronConfiguration(                             8,                  18,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.360784314f,      0.721568627f,     0.819607843f, ElementStructureFactor(      17.355500,        6.728600,        5.549300,        3.537500,        1.938400,       16.562300,        0.226100,       39.397200,        2.825000            )))),// Krypton
    // 5th period              element name                       atomic number  period  group  Symbol    Element Name        mass                  gyromagnetic ratio  covalent radiu  vdw radius  melting poin  boiling poi                     electro negativity  ionization poten                                  sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d   7f
    e_Rubidium(      AddEnum(  "Rubidium"     , Element(             37,      5,     1, "Rb"    , "Rubidium"     ,   85.47,                           87532000,           2.16,       1.70,       312.64,         961,                                   0.8,             4.18, ElectronConfiguration(                             1,                   1,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.439215686f,      0.180392157f,     0.690196078f, ElementStructureFactor(      17.178400,        9.643500,        5.139900,        1.529200,        1.788800,       17.315100,        0.274800,      164.934000,        3.487300            )))),// Rubidium
    Strontium(     AddEnum(  "Strontium"    , Element(             38,      5,     2, "Sr"    , "Strontium"    ,   87.62,                           11594000,           1.91,       1.70,         1042,        1657,                                     1,              5.7, ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.000000000f,      1.000000000f,     0.000000000f, ElementStructureFactor(      17.566300,        9.818400,        5.422000,        2.669400,        1.556400,       14.098800,        0.166400,      132.376000,        2.506400            )))),// Strontium
    e_Yttrium(       AddEnum(  "Yttrium"      , Element(             39,      5,     0, "Y"     , "Yttrium"      ,   88.91,                           13108000,           1.62,       1.70,         1799,        3611,                                   1.2,             6.38, ElectronConfiguration( util::GetUndefined< size_t>(),                   3,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   1,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.580392157f,      1.000000000f,     1.000000000f, ElementStructureFactor(      17.776000,       10.294600,        5.726290,        3.265880,        1.402900,       12.800600,        0.125599,      104.354000,        1.912130            )))),// Yttrium
    e_Zirconium(     AddEnum(  "Zirconium"    , Element(             40,      5,     0, "Zr"    , "Zirconium"    ,   91.22,                           24868000,           1.45,       1.70,         2125,        4650,                                   1.3,             6.84, ElectronConfiguration( util::GetUndefined< size_t>(),                   4,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   2,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.580392157f,      0.878431373f,     0.878431373f, ElementStructureFactor(      17.876500,       10.948000,        5.417320,        3.657210,        1.276180,       11.916000,        0.117622,       87.662700,        2.069290            )))),// Zirconium
    e_Niobium(       AddEnum(  "Niobium"      , Element(             41,      5,     0, "Nb"    , "Niobium"      ,   92.91,                           65476000,           1.34,       1.70,         2741,        5017,                                   1.6,             6.88, ElectronConfiguration( util::GetUndefined< size_t>(),                   5,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   4,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.450980392f,      0.760784314f,     0.788235294f, ElementStructureFactor(      17.614200,       12.014400,        4.041830,        3.533460,        1.188650,       11.766000,        0.204785,       69.795700,        3.755910            )))),// Niobium
    e_Molybdenum(    AddEnum(  "Molybdenum"   , Element(             42,      5,     0, "Mo"    , "Molybdenum"   ,   95.94,                           17433000,            1.3,       1.70,         2890,        4885,                                   2.2,              7.1, ElectronConfiguration( util::GetUndefined< size_t>(),                   6,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   5,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.329411765f,      0.709803922f,     0.709803922f, ElementStructureFactor(       3.702500,       17.235600,       12.887600,        3.742900,        0.277200,        1.095800,       11.004000,       61.658400,        4.387500            )))),// Molybdenum
    e_Technetium(    AddEnum(  "Technetium"   , Element(             43,      5,     0, "Tc"    , "Technetium"   ,      98,                           60211000,           1.27,       1.70,         2473,        5150,                                   1.9,             7.28, ElectronConfiguration( util::GetUndefined< size_t>(),                   7,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   6,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.231372549f,      0.619607843f,     0.619607843f, ElementStructureFactor(      19.130100,       11.094800,        4.649010,        2.712630,        0.864132,        8.144870,       21.570700,       86.847200,        5.404280            )))),// Technetium
    e_Ruthenium(     AddEnum(  "Ruthenium"    , Element(             44,      5,     0, "Ru"    , "Ruthenium"    ,  101.07,                           13833000,           1.25,       1.70,         2523,        4173,                                   2.2,             7.37, ElectronConfiguration( util::GetUndefined< size_t>(),                   8,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   7,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.141176471f,      0.560784314f,     0.560784314f, ElementStructureFactor(      19.267400,       12.918200,        4.863370,        1.567560,        0.808520,        8.434670,       24.799700,       94.292800,        5.378740            )))),// Ruthenium
    e_Rhodium(       AddEnum(  "Rhodium"      , Element(             45,      5,     0, "Rh"    , "Rhodium"      ,  102.91,                            8520100,           1.25,       1.70,         2239,        4000,                                   2.3,             7.46, ElectronConfiguration( util::GetUndefined< size_t>(),                   9,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,   8,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.039215686f,      0.490196078f,     0.549019608f, ElementStructureFactor(      19.295700,       14.350100,        4.734250,        1.289180,        0.751536,        8.217580,       25.874900,       98.606200,        5.328000            )))),// Rhodium
    e_Palladium(     AddEnum(  "Palladium"    , Element(             46,      5,     0, "Pd"    , "Palladium"    ,  106.42,                           12241000,           1.28,       1.70,         1825,        3237,                                   2.2,             8.34, ElectronConfiguration( util::GetUndefined< size_t>(),                  10,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.000000000f,      0.411764706f,     0.521568627f, ElementStructureFactor(      19.331900,       15.501700,        5.295370,        0.605844,        0.698655,        7.989290,       25.205200,       76.898600,        5.265930            )))),// Palladium
    Silver(        AddEnum(  "Silver"       , Element(             47,      5,     0, "Ag"    , "Silver"       ,  107.87,                           12450000,           1.34,       1.70,         1234,        2436,                                   1.9,             7.58, ElectronConfiguration( util::GetUndefined< size_t>(),                  11,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.752941176f,      0.752941176f,     0.752941176f, ElementStructureFactor(      19.280800,       16.688500,        4.804500,        1.046300,        0.644600,        7.472600,       24.660500,       99.815600,        5.179000            )))),// Silver
    e_Cadmium(       AddEnum(  "Cadmium"      , Element(             48,      5,     0, "Cd"    , "Cadmium"      ,  112.41,                           59329000,           1.48,       1.70,       594.18,        1038,                                   1.7,             8.99, ElectronConfiguration( util::GetUndefined< size_t>(),                  12,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    1.000000000f,      0.850980392f,     0.560784314f, ElementStructureFactor(      19.221400,       17.644400,        4.461000,        1.602900,        0.594600,        6.908900,       24.700800,       87.482500,        5.069400            )))),// Cadmium
    e_Indium(        AddEnum(  "Indium"       , Element(             49,      5,     3, "In"    , "Indium"       ,  114.82,                           58619000,           1.44,       1.70,       429.76,        2346,                                   1.8,             5.79, ElectronConfiguration(                             3,                  13,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   1,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.650980392f,      0.458823529f,     0.450980392f, ElementStructureFactor(      19.162400,       18.559600,        4.294800,        2.039600,        0.547600,        6.377600,       25.849900,       92.802900,        4.939100            )))),// Indium
    e_Tin(           AddEnum(  "Tin"          , Element(             50,      5,     4, "Sn"    , "Tin"          ,  118.71,                           99757000,           1.41,       1.70,       505.06,        2543,                                     2,             7.34, ElectronConfiguration(                             4,                  14,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   2,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.400000000f,      0.501960784f,     0.501960784f, ElementStructureFactor(      19.188900,       19.100500,        4.458500,        2.466300,        5.830300,        0.503100,       26.890900,       83.957100,        4.782100            )))),// Tin
    e_Antimony(      AddEnum(  "Antimony"     , Element(             51,      5,     5, "Sb"    , "Antimony"     ,  121.76,                           64018000,            1.4,       2.20,        903.9,        1860,                                   2.1,             8.64, ElectronConfiguration(                             5,                  15,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   3,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.619607843f,      0.388235294f,     0.709803922f, ElementStructureFactor(      19.641800,       19.045500,        5.037100,        2.682700,        5.303400,        0.460700,       27.907400,       75.282500,        4.590900            )))),// Antimony
    e_Tellurium(     AddEnum(  "Tellurium"    , Element(             52,      5,     6, "Te"    , "Tellurium"    ,   127.6,                           84399000,           1.36,       2.20,       722.65,        1261,                                   2.1,             9.01, ElectronConfiguration(                             6,                  16,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   4,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.831372549f,      0.478431373f,     0.000000000f, ElementStructureFactor(      19.964400,       19.013800,        6.144870,        2.523900,        4.817420,        0.420885,       28.528400,       70.840300,        4.352000            )))),// Tellurium
    e_Iodine(        AddEnum(  "Iodine"       , Element(             53,      5,     7, "I"     , "Iodine"       ,   126.9,                           53526000,           1.33,       2.15,        386.7,       458.4,                                   2.7,            10.45, ElectronConfiguration(                             7,                  17,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   5,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.580392157f,      0.000000000f,     0.580392157f, ElementStructureFactor(      20.147200,       18.994900,        7.513800,        2.273500,        4.347000,        0.381400,       27.766000,       66.877600,        4.071200            )))),// Iodine
    e_Xenon(         AddEnum(  "Xenon"        , Element(             54,      5,     8, "Xe"    , "Xenon"        ,  131.29,                           73988000,           1.31,       2.16,        161.3,      165.03,         util::GetUndefined< double>(),            12.13, ElectronConfiguration(                             8,                  18,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   6,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0 ),    0.258823529f,      0.619607843f,     0.690196078f, ElementStructureFactor(      20.293300,       19.029800,        8.976700,        1.990000,        3.928200,        0.344000,       26.465900,       64.265800,        3.711800            )))),// Xenon
    // 6th period              element name                       atomic number  period  group  Symbol    Element Name        mass                  gyromagnetic ratio  covalent radiu  vdw radius  melting poin  boiling poi                     electro negativity  ionization poten                                  sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d   7f
    e_Cesium(        AddEnum(  "Cesium"       , Element(             55,      6,     1, "Cs"    , "Cesium"       ,  132.91,                                  0,           2.35,       1.70,       301.55,         944,                                   0.8,             3.89, ElectronConfiguration(                             1,                   1,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   6,   0,   0,   1,   0,   0,   0,   0,   0,   0,   0 ),    0.341176471f,      0.090196078f,     0.560784314f, ElementStructureFactor(      20.389200,       19.106200,       10.662000,        1.495300,        3.569000,        0.310700,       24.387900,      213.904000,        3.335200            )))),// Cesium
    e_Barium(        AddEnum(  "Barium"       , Element(             56,      6,     2, "Ba"    , "Barium"       ,  137.33,                                  0,           1.98,       1.70,         1002,        2171,                                   0.9,             5.21, ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.000000000f,      0.788235294f,     0.000000000f, ElementStructureFactor(      20.336100,       19.297000,       10.888000,        2.695900,        3.216000,        0.275600,       20.207300,      167.202000,        2.773100            )))),// Barium
    e_Lanthanum(     AddEnum(  "Lanthanum"    , Element(             57,      6,     0, "La"    , "Lanthanum"    ,  138.91,                                  0,           1.69,       1.70,         1193,        3730,                                   1.1,             5.58, ElectronConfiguration( util::GetUndefined< size_t>(),                   3,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   0,   2,   6,   1,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.439215686f,      0.831372549f,     1.000000000f, ElementStructureFactor(      20.578000,       19.599000,       11.372700,        3.287190,        2.948170,        0.244475,       18.772600,      133.124000,        2.146780            )))),// Lanthanum
    e_Cerium(        AddEnum(  "Cerium"       , Element(             58,      6,     0, "Ce"    , "Cerium"       ,  140.12,                                  0,           1.65,       1.70,         1071,        3699,                                   1.1,             5.54, ElectronConfiguration( util::GetUndefined< size_t>(),                   4,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   1,   2,   6,   1,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    1.000000000f,      1.000000000f,     0.780392157f, ElementStructureFactor(      21.167100,       19.769500,       11.851300,        3.330490,        2.812190,        0.226836,       17.608300,      127.113000,        1.862640            )))),// Cerium
    e_Praseodymium(  AddEnum(  "Praseodymium" , Element(             59,      6,     0, "Pr"    , "Praseodymium" ,  140.91,                                  0,           1.65,       1.70,         1204,        3785,                                   1.1,             5.46, ElectronConfiguration( util::GetUndefined< size_t>(),                   5,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   3,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.850980392f,      1.000000000f,     0.780392157f, ElementStructureFactor(      22.044000,       19.669700,       12.385600,        2.824280,        2.773930,        0.222087,       16.766900,      143.644000,        2.058300            )))),// Praseodymium
    e_Neodymium(     AddEnum(  "Neodymium"    , Element(             60,      6,     0, "Nd"    , "Neodymium"    ,  144.24,                                  0,           1.64,       1.70,         1289,        3341,                                   1.1,             5.53, ElectronConfiguration( util::GetUndefined< size_t>(),                   6,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   4,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.780392157f,      1.000000000f,     0.780392157f, ElementStructureFactor(      22.684500,       19.684700,       12.774000,        2.851370,        2.662480,        0.210628,       15.885000,      137.903000,        1.984860            )))),// Neodymium
    e_Promethium(    AddEnum(  "Promethium"   , Element(             61,      6,     0, "Pm"    , "Promethium"   ,     145,                                  0,           1.63,       1.70,         1204,        3785,                                   1.1,             5.55, ElectronConfiguration( util::GetUndefined< size_t>(),                   7,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   5,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.639215686f,      1.000000000f,     0.780392157f, ElementStructureFactor(      23.340500,       19.609500,       13.123500,        2.875160,        2.562700,        0.202088,       15.100900,      132.721000,        2.028760            )))),// Promethium
    Samarium(      AddEnum(  "Samarium"     , Element(             62,      6,     0, "Sm"    , "Samarium"     ,  150.36,                                  0,           1.62,       1.70,         1345,        2064,                                   1.2,             5.64, ElectronConfiguration( util::GetUndefined< size_t>(),                   8,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   6,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.560784314f,      1.000000000f,     0.780392157f, ElementStructureFactor(      24.004200,       19.425800,       13.439600,        2.896040,        2.472740,        0.196451,       14.399600,      128.007000,        2.209630            )))),// Samarium
    e_Europium(      AddEnum(  "Europium"     , Element(             63,      6,     0, "Eu"    , "Europium"     ,  151.97,                                  0,           1.85,       1.70,         1095,        1870,                                   1.2,             5.67, ElectronConfiguration( util::GetUndefined< size_t>(),                   9,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   7,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.380392157f,      1.000000000f,     0.780392157f, ElementStructureFactor(      24.627400,       19.088600,       13.760300,        2.922700,        2.387900,        0.194200,       13.754600,      123.174000,        2.574500            )))),// Europium
    e_Gadolinium(    AddEnum(  "Gadolinium"   , Element(             64,      6,     0, "Gd"    , "Gadolinium"   ,  157.25,                                  0,           1.61,       1.70,         1585,        3539,                                   1.2,             6.15, ElectronConfiguration( util::GetUndefined< size_t>(),                  10,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   7,   2,   6,   1,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.270588235f,      1.000000000f,     0.780392157f, ElementStructureFactor(      25.070900,       19.079800,       13.851800,        3.545450,        2.253410,        0.181951,       12.933100,      101.398000,        2.419600            )))),// Gadolinium
    e_Terbium(       AddEnum(  "Terbium"      , Element(             65,      6,     0, "Tb"    , "Terbium"      ,  158.93,                                  0,           1.59,       1.70,         1630,        3296,                                   1.2,             5.86, ElectronConfiguration( util::GetUndefined< size_t>(),                  11,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,   9,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.188235294f,      1.000000000f,     0.780392157f, ElementStructureFactor(      25.897600,       18.218500,       14.316700,        2.953540,        2.242560,        0.196143,       12.664800,      115.362000,        3.583240            )))),// Terbium
    e_Dysprosium(    AddEnum(  "Dysprosium"   , Element(             66,      6,     0, "Dy"    , "Dysprosium"   ,   162.5,                                  0,           1.59,       1.70,         1685,        2835,                                   1.2,             5.94, ElectronConfiguration( util::GetUndefined< size_t>(),                  12,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  10,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.121568627f,      1.000000000f,     0.780392157f, ElementStructureFactor(      26.507000,       17.638300,       14.559600,        2.965770,        2.180200,        0.202172,       12.189900,      111.874000,        4.297280            )))),// Dysprosium
    e_Holmium(       AddEnum(  "Holmium"      , Element(             67,      6,     0, "Ho"    , "Holmium"      ,  164.93,                                  0,           1.58,       1.70,         1743,        2968,                                   1.2,             6.02, ElectronConfiguration( util::GetUndefined< size_t>(),                  13,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  11,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.000000000f,      1.000000000f,     0.611764706f, ElementStructureFactor(      26.904900,       17.294000,       14.558300,        3.638370,        2.070510,        0.197940,       11.440700,       92.656600,        4.567960            )))),// Holmium
    e_Erbium(        AddEnum(  "Erbium"       , Element(             68,      6,     0, "Er"    , "Erbium"       ,  167.26,                                  0,           1.57,       1.70,         1795,        3136,                                   1.2,              6.1, ElectronConfiguration( util::GetUndefined< size_t>(),                  14,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  12,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.000000000f,      0.901960784f,     0.458823529f, ElementStructureFactor(      27.656300,       16.428500,       14.977900,        2.982330,        2.073560,        0.223545,       11.360400,      105.703000,        5.920460            )))),// Erbium
    e_Thulium(       AddEnum(  "Thulium"      , Element(             69,      6,     0, "Tm"    , "Thulium"      ,  168.93,                                  0,           1.56,       1.70,         1818,        2220,                                   1.3,             6.18, ElectronConfiguration( util::GetUndefined< size_t>(),                  15,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  13,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.000000000f,      0.831372549f,     0.321568627f, ElementStructureFactor(      28.181900,       15.885100,       15.154200,        2.987060,        2.028590,        0.238849,       10.997500,      102.961000,        6.756210            )))),// Thulium
    e_Ytterbium(     AddEnum(  "Ytterbium"    , Element(             70,      6,     0, "Yb"    , "Ytterbium"    ,  173.04,                                  0,           1.74,       1.70,         1097,        1467,                                   1.1,             6.25, ElectronConfiguration( util::GetUndefined< size_t>(),                  16,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   0,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.000000000f,      0.749019608f,     0.219607843f, ElementStructureFactor(      28.664100,       15.434500,       15.308700,        2.989630,        1.988900,        0.257119,       10.664700,      100.417000,        7.566720            )))),// Ytterbium
    e_Lutetium(      AddEnum(  "Lutetium"     , Element(             71,      6,     0, "Lu"    , "Lutetium"     ,  174.97,                                  0,           1.56,       1.70,         1936,        3668,                                   1.3,             5.43, ElectronConfiguration( util::GetUndefined< size_t>(),                  17,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   1,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.000000000f,      0.670588235f,     0.141176471f, ElementStructureFactor(      28.947600,       15.220800,       15.100000,        3.716010,        1.901820,        9.985190,        0.261033,       84.329800,        7.976280            )))),// Lutetium
    e_Hafnium(       AddEnum(  "Hafnium"      , Element(             72,      6,     0, "Hf"    , "Hafnium"      ,  178.49,                                  0,           1.44,       1.70,         2500,        4876,                                   1.3,             6.65, ElectronConfiguration( util::GetUndefined< size_t>(),                  18,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   2,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.301960784f,      0.760784314f,     1.000000000f, ElementStructureFactor(      29.144000,       15.172600,       14.758600,        4.300130,        1.832620,        9.599900,        0.275116,       72.029000,        8.581540            )))),// Hafnium
    e_Tantalum(      AddEnum(  "Tantalum"     , Element(             73,      6,     0, "Ta"    , "Tantalum"     ,  180.95,                                  0,           1.34,       1.70,         3269,        5698,                                   1.5,             7.89, ElectronConfiguration( util::GetUndefined< size_t>(),                  19,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   3,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.301960784f,      0.650980392f,     1.000000000f, ElementStructureFactor(      29.202400,       15.229300,       14.513500,        4.764920,        1.773330,        9.370460,        0.295977,       63.364400,        9.243540            )))),// Tantalum
    e_Tungsten(      AddEnum(  "Tungsten"     , Element(             74,      6,     0, "W"     , "Tungsten"     ,  183.85,                                  0,            1.3,       1.70,         3680,        5928,                                   2.4,             7.98, ElectronConfiguration( util::GetUndefined< size_t>(),                  20,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   4,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.129411765f,      0.580392157f,     0.839215686f, ElementStructureFactor(      29.081800,       15.430000,       14.432700,        5.119820,        1.720290,        9.225900,        0.321703,       57.056000,        9.887500            )))),// Tungsten
    e_Rhenium(       AddEnum(  "Rhenium"      , Element(             75,      6,     0, "Re"    , "Rhenium"      ,  186.21,                                  0,           1.28,       1.70,         3453,        5900,                                   1.9,             7.88, ElectronConfiguration( util::GetUndefined< size_t>(),                  21,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   5,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.149019608f,      0.490196078f,     0.670588235f, ElementStructureFactor(      28.762100,       15.718900,       14.556400,        5.441740,        1.671910,        9.092270,        0.350500,       52.086100,       10.472000            )))),// Rhenium
    e_Osmium(        AddEnum(  "Osmium"       , Element(             76,      6,     0, "Os"    , "Osmium"       ,   190.2,                                  0,           1.26,       1.70,         3300,        5285,                                   2.2,              8.7, ElectronConfiguration( util::GetUndefined< size_t>(),                  22,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   6,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.149019608f,      0.400000000f,     0.588235294f, ElementStructureFactor(      28.189400,       16.155000,       14.930500,        5.675890,        1.629030,        8.979480,        0.382661,       48.164700,       11.000500            )))),// Osmium
    e_Iridium(       AddEnum(  "Iridium"      , Element(             77,      6,     0, "Ir"    , "Iridium"      ,  192.22,                                  0,           1.27,       1.70,         2716,        4701,                                   2.2,              9.1, ElectronConfiguration( util::GetUndefined< size_t>(),                  23,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   7,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.090196078f,      0.329411765f,     0.529411765f, ElementStructureFactor(      27.304900,       16.729600,       15.611500,        5.833770,        1.592790,        8.865530,        0.417916,       45.001100,       11.472200            )))),// Iridium
    e_Platinum(      AddEnum(  "Platinum"     , Element(             78,      6,     0, "Pt"    , "Platinum"     ,  195.08,                                  0,            1.3,       1.72,         2045,        4100,                                   2.3,                9, ElectronConfiguration( util::GetUndefined< size_t>(),                  24,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,   9,   0,   1,   0,   0,   0,   0,   0,   0,   0 ),    0.815686275f,      0.815686275f,     0.878431373f, ElementStructureFactor(      27.005900,       17.763900,       15.713100,        5.783700,        1.512930,        8.811740,        0.424593,       38.610300,       11.688300            )))),// Platinum
    e_Gold(          AddEnum(  "Gold"         , Element(             79,      6,     0, "Au"    , "Gold"         ,  196.97,                                  0,           1.34,       1.66,      1337.58,        3080,                                   2.5,             9.23, ElectronConfiguration( util::GetUndefined< size_t>(),                  25,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   1,   0,   0,   0,   0,   0,   0,   0 ),    1.000000000f,      0.819607843f,     0.137254902f, ElementStructureFactor(      16.881900,       18.591300,       25.558200,        5.860000,        0.461100,        8.621600,        1.482600,       36.395600,       12.065800            )))),// Gold
    e_Mercury(       AddEnum(  "Mercury"      , Element(             80,      6,     0, "Hg"    , "Mercury"      ,  200.59,                                  0,           1.49,       1.55,       234.28,         630,                                     2,            10.44, ElectronConfiguration( util::GetUndefined< size_t>(),                  26,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   0,   0,   0,   0,   0,   0,   0 ),    0.721568627f,      0.721568627f,     0.815686275f, ElementStructureFactor(      20.680900,       19.041700,       21.657500,        5.967600,        0.545000,        8.448400,        1.572900,       38.324600,       12.608900            )))),// Mercury
    e_Thallium(      AddEnum(  "Thallium"     , Element(             81,      6,     3, "Tl"    , "Thallium"     ,  204.38,                                  0,           1.48,       1.70,          577,        1746,                                     2,             6.11, ElectronConfiguration(                             3,                  27,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   1,   0,   0,   0,   0,   0,   0 ),    0.650980392f,      0.329411765f,     0.301960784f, ElementStructureFactor(      27.544600,       19.158400,       15.538000,        5.525930,        0.655150,        8.707510,        1.963470,       45.814900,       13.174600            )))),// Thallium
    e_Lead(          AddEnum(  "Lead"         , Element(             82,      6,     4, "Pb"    , "Lead"         ,   207.2,                                  0,           1.47,       1.70,        600.6,        2013,                                   2.3,             7.42, ElectronConfiguration(                             4,                  28,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   2,   0,   0,   0,   0,   0,   0 ),    0.341176471f,      0.349019608f,     0.380392157f, ElementStructureFactor(      31.061700,       13.063700,       18.442000,        5.969600,        0.690200,        2.357600,        8.618000,       47.257900,       13.411800            )))),// Lead
    e_Bismuth(       AddEnum(  "Bismuth"      , Element(             83,      6,     5, "Bi"    , "Bismuth"      ,  208.98,                                  0,           1.46,       1.70,       544.52,        1837,                                     2,             7.29, ElectronConfiguration(                             5,                  29,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   3,   0,   0,   0,   0,   0,   0 ),    0.619607843f,      0.309803922f,     0.709803922f, ElementStructureFactor(      33.368900,       12.951000,       16.587700,        6.469200,        0.704000,        2.923800,        8.793700,       48.009300,       13.578200            )))),// Bismuth
    e_Polonium(      AddEnum(  "Polonium"     , Element(             84,      6,     6, "Po"    , "Polonium"     ,     209,                                  0,           1.46,       1.70,          527,        1235,                                     2,             8.42, ElectronConfiguration(                             6,                  30,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   4,   0,   0,   0,   0,   0,   0 ),    0.670588235f,      0.360784314f,     0.000000000f, ElementStructureFactor(      34.672600,       15.473300,       13.113800,        7.025880,        0.700999,        3.550780,        9.556420,       47.004500,       13.677000            )))),// Polonium
    e_Astatine(      AddEnum(  "Astatine"     , Element(             85,      6,     7, "At"    , "Astatine"     ,     210,                                  0,           1.45,       1.70,          575,         610,                                   2.2,             9.65, ElectronConfiguration(                             7,                  31,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   5,   0,   0,   0,   0,   0,   0 ),    0.458823529f,      0.309803922f,     0.270588235f, ElementStructureFactor(      35.316300,       19.021100,        9.498870,        7.425180,        0.685870,        3.974580,       11.382400,       45.471500,       13.710800            )))),// Astatine
    e_Radon(         AddEnum(  "Radon"        , Element(             86,      6,     8, "Rn"    , "Radon"        ,     222,                                  0,              0,       1.70,          202,         211,         util::GetUndefined< double>(),            10.75, ElectronConfiguration(                             8,                  32,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   6,   0,   0,   0,   0,   0,   0 ),    0.258823529f,      0.509803922f,     0.588235294f, ElementStructureFactor(      35.563100,       21.281600,        8.003700,        7.443300,        0.663100,        4.069100,       14.042200,       44.247300,       13.690500            )))),// Radon
    // 7th period              element name                       atomic number  period  group  Symbol    Element Name        mass                  gyromagnetic ratio  covalent radiu  vdw radius  melting poin  boiling poi                     electro negativity  ionization poten                                  sp valence electrons  spd valence electro   1s   1p   1d   1f   2s   2p   2d   2f   3s   3p   3d   3f   4s   4p   4d   4f   5s   5p   5d   5f   6s   6p   6d   6f   7s   7p   7d   7f
    e_Francium(      AddEnum(  "Francium"     , Element(             87,      7,     1, "Fr"    , "Francium"     ,     223,                                  0,              0,       1.70,          300,         950,                                   0.7,             3.83, ElectronConfiguration(                             1,                   1,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   6,   0,   0,   1,   0,   0,   0 ),    0.258823529f,      0.000000000f,     0.400000000f, ElementStructureFactor(      35.929900,       23.054700,       12.143900,        2.112530,        0.646453,        4.176190,       23.105200,      150.645000,       13.724700            )))),// Francium
    e_Radium(        AddEnum(  "Radium"       , Element(             88,      7,     2, "Ra"    , "Radium"       ,  226.03,                                  0,              0,       1.70,          973,        1809,                                   0.9,             5.28, ElectronConfiguration(                             2,                   2,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   6,   0,   0,   2,   0,   0,   0 ),    0.000000000f,      0.490196078f,     0.000000000f, ElementStructureFactor(      35.763000,       22.906400,       12.473900,        3.210970,        0.616341,        3.871350,       19.988700,      142.325000,       13.621100            )))),// Radium
    e_Actinium(      AddEnum(  "Actinium"     , Element(             89,      7,     0, "Ac"    , "Actinium"     ,     227,                                  0,              0,       1.70,         1323,        3473,                                   1.1,             5.17, ElectronConfiguration( util::GetUndefined< size_t>(),                   3,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   1,   6,   1,   0,   2,   0,   0,   0 ),    0.439215686f,      0.670588235f,     0.980392157f, ElementStructureFactor(      35.659700,       23.103200,       12.597700,        4.086550,        0.589092,        3.651550,       18.599000,      117.020000,       13.526600            )))),// Actinium
    e_Thorium(       AddEnum(  "Thorium"      , Element(             90,      7,     0, "Th"    , "Thorium"      ,  232.04,                                  0,           1.65,       1.70,         2028,        5061,                                   1.3,             6.08, ElectronConfiguration( util::GetUndefined< size_t>(),                   4,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   0,   2,   6,   2,   0,   2,   0,   0,   0 ),    0.000000000f,      0.729411765f,     1.000000000f, ElementStructureFactor(      35.564500,       23.421900,       12.747300,        4.807030,        0.563359,        3.462040,       17.830900,       99.172200,       13.431400            )))),// Thorium
    e_Protactinium(  AddEnum(  "Protactinium" , Element(             91,      7,     0, "Pa"    , "Protactinium" ,  213.04,                                  0,              0,       1.70,         2113,        4300,                                   1.5,             5.89, ElectronConfiguration( util::GetUndefined< size_t>(),                   5,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   2,   1,   6,   1,   0,   2,   0,   0,   0 ),    0.000000000f,      0.631372549f,     1.000000000f, ElementStructureFactor(      35.884700,       23.294800,       14.189100,        4.172870,        0.547751,        3.415190,       16.923500,      105.251000,       13.428700            )))),// Protactinium
    Uranium(       AddEnum(  "Uranium"      , Element(             92,      7,     0, "U"     , "Uranium"      ,  238.03,                                  0,           1.42,       1.70,         1405,        4407,                                   1.4,             6.05, ElectronConfiguration( util::GetUndefined< size_t>(),                   6,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   3,   1,   6,   1,   0,   2,   0,   0,   0 ),    0.000000000f,      0.560784314f,     1.000000000f, ElementStructureFactor(      36.022800,       23.412800,       14.949100,        4.188000,        0.529300,        3.325300,       16.092700,      100.613000,       13.396600            )))),// Uranium
    e_Neptunium(     AddEnum(  "Neptunium"    , Element(             93,      7,     0, "Np"    , "Neptunium"    ,  237.05,                                  0,              0,       1.70,          913,        4175,                                   1.4,             6.19, ElectronConfiguration( util::GetUndefined< size_t>(),                   7,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   4,   1,   6,   1,   0,   2,   0,   0,   0 ),    0.000000000f,      0.501960784f,     1.000000000f, ElementStructureFactor(      36.187400,       23.596400,       15.640200,        4.185500,        0.511929,        3.253960,       15.362200,       97.490800,       13.357300            )))),// Neptunium
    e_Plutonium(     AddEnum(  "Plutonium"    , Element(             94,      7,     0, "Pu"    , "Plutonium"    ,     244,                                  0,              0,       1.70,          913,        3503,                                   1.3,             6.06, ElectronConfiguration( util::GetUndefined< size_t>(),                   8,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   6,   2,   6,   0,   0,   2,   0,   0,   0 ),    0.000000000f,      0.419607843f,     1.000000000f, ElementStructureFactor(      35.510300,       22.578700,       12.776600,        4.921590,        0.498626,        2.966270,       11.948400,       22.750200,       13.211600            )))),// Plutonium
    e_Americium(     AddEnum(  "Americium"    , Element(             95,      7,     0, "Am"    , "Americium"    ,     243,                                  0,              0,       1.70,         1267,        2880,                                   1.3,             5.99, ElectronConfiguration( util::GetUndefined< size_t>(),                   9,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   7,   2,   6,   0,   0,   2,   0,   0,   0 ),    0.329411765f,      0.360784314f,     0.949019608f, ElementStructureFactor(      36.670600,       24.099200,       17.341500,        3.493310,        0.483629,        3.206470,       14.313600,      102.273000,       13.359200            )))),// Americium
    e_Curium(        AddEnum(  "Curium"       , Element(             96,      7,     0, "Cm"    , "Curium"       ,     247,                                  0,              0,       1.70,         1340,        3383,                                   1.3,             6.02, ElectronConfiguration( util::GetUndefined< size_t>(),                  10,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   7,   1,   6,   1,   0,   2,   0,   0,   0 ),    0.470588235f,      0.360784314f,     0.890196078f, ElementStructureFactor(      36.648800,       24.409600,       17.399000,        4.216650,        0.465154,        3.089970,       13.434600,       88.483400,       13.288700            )))),// Curium
    e_Berkelium(     AddEnum(  "Berkelium"    , Element(             97,      7,     0, "Bk"    , "Berkelium"    ,     247,                                  0,              0,       1.70,         1259,           0,                                   1.3,             6.23, ElectronConfiguration( util::GetUndefined< size_t>(),                  11,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,   9,   2,   6,   0,   0,   2,   0,   0,   0 ),    0.541176471f,      0.309803922f,     0.890196078f, ElementStructureFactor(      36.788100,       24.773600,       17.891900,        4.232840,        0.451018,        3.046190,       12.894600,       86.003000,       13.275400            )))),// Berkelium
    e_Californium(   AddEnum(  "Californium"  , Element(             98,      7,     0, "Cf"    , "Californium"  ,     251,                                  0,              0,       1.70,         1173,           0,                                   1.3,              6.3, ElectronConfiguration( util::GetUndefined< size_t>(),                  12,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  10,   2,   6,   0,   0,   2,   0,   0,   0 ),    0.631372549f,      0.211764706f,     0.831372549f, ElementStructureFactor(      36.918500,       25.199500,       18.331700,        4.243910,        0.437533,        3.007750,       12.404400,       83.788100,       13.267400            )))),// Californium
    e_Einsteinium(   AddEnum(  "Einsteinium"  , Element(             99,      7,     0, "Es"    , "Einsteinium"  ,     252,                                  0,              0,       1.70,         1133,           0,                                   1.3,             6.42, ElectronConfiguration( util::GetUndefined< size_t>(),                  13,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  11,   2,   6,   0,   0,   2,   0,   0,   0 ),    0.701960784f,      0.121568627f,     0.831372549f, ElementStructureFactor(                                                                                                                                                                   )))),// Einsteinium
    e_Fermium(       AddEnum(  "Fermium"      , Element(            100,      7,     0, "Fm"    , "Fermium"      ,     257,                                  0,              0,       1.70,            0,           0,                                   1.3,              6.5, ElectronConfiguration( util::GetUndefined< size_t>(),                  14,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  12,   2,   6,   0,   0,   2,   0,   0,   0 ),    0.701960784f,      0.121568627f,     0.729411765f, ElementStructureFactor(                                                                                                                                                                   )))),// Fermium
    e_Mendelevium(   AddEnum(  "Mendelevium"  , Element(            101,      7,     0, "Md"    , "Mendelevium"  ,     258,                                  0,              0,       1.70,            0,           0,                                   1.3,             6.58, ElectronConfiguration( util::GetUndefined< size_t>(),                  15,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  13,   2,   6,   0,   0,   2,   0,   0,   0 ),    0.701960784f,      0.050980392f,     0.650980392f, ElementStructureFactor(                                                                                                                                                                   )))),// Mendelevium
    e_Nobelium(      AddEnum(  "Nobelium"     , Element(            102,      7,     0, "No"    , "Nobelium"     ,     259,                                  0,              0,       1.70,            0,           0,                                   1.3,             6.65, ElectronConfiguration( util::GetUndefined< size_t>(),                  16,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   0,   0,   2,   0,   0,   0 ),    0.741176471f,      0.050980392f,     0.529411765f, ElementStructureFactor(                                                                                                                                                                   )))),// Nobelium
    e_Lawrencium(    AddEnum(  "Lawrencium"   , Element(            103,      7,     0, "Lr"    , "Lawrencium"   ,     260,                                  0,              0,       1.70,            0,           0,                                     0,                0, ElectronConfiguration( util::GetUndefined< size_t>(),                  17,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   1,   0,   2,   0,   0,   0 ),    0.780392157f,      0.000000000f,     0.400000000f, ElementStructureFactor(                                                                                                                                                                   )))),// Lawrencium
    e_Rutherfordium( AddEnum(  "Rutherfordium", Element(            104,      7,     0, "Rf"    , "Rutherfordium",     261,                                  0,              0,       1.70,            0,           0,                                     0,                0, ElectronConfiguration( util::GetUndefined< size_t>(),                  18,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   2,   0,   2,   0,   0,   0 ),    0.800000000f,      0.000000000f,     0.349019608f, ElementStructureFactor(                                                                                                                                                                   )))),// Rutherfordiu
    e_Dubnium(       AddEnum(  "Dubnium"      , Element(            105,      7,     0, "Db"    , "Dubnium"      ,     262,                                  0,              0,       1.70,            0,           0,                                     0,                0, ElectronConfiguration( util::GetUndefined< size_t>(),                  19,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   3,   0,   2,   0,   0,   0 ),    0.819607843f,      0.000000000f,     0.309803922f, ElementStructureFactor(                                                                                                                                                                   )))),// Dubnium
    Seaborgium(    AddEnum(  "Seaborgium"   , Element(            106,      7,     0, "Sg"    , "Seaborgium"   ,     263,                                  0,              0,       1.70,            0,           0,                                     0,                0, ElectronConfiguration( util::GetUndefined< size_t>(),                  20,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   4,   0,   2,   0,   0,   0 ),    0.850980392f,      0.000000000f,     0.270588235f, ElementStructureFactor(                                                                                                                                                                   )))),// Seaborgium
    e_Bohrium(       AddEnum(  "Bohrium"      , Element(            107,      7,     0, "Bh"    , "Bohrium"      ,     262,                                  0,              0,       1.70,            0,           0,                                     0,                0, ElectronConfiguration( util::GetUndefined< size_t>(),                  21,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   5,   0,   2,   0,   0,   0 ),    0.878431373f,      0.000000000f,     0.219607843f, ElementStructureFactor(                                                                                                                                                                   )))),// Bohrium
    e_Hassium(       AddEnum(  "Hassium"      , Element(            108,      7,     0, "Hs"    , "Hassium"      ,     265,                                  0,              0,       1.70,            0,           0,                                     0,                0, ElectronConfiguration( util::GetUndefined< size_t>(),                  22,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   6,   0,   2,   0,   0,   0 ),    0.901960784f,      0.000000000f,     0.180392157f, ElementStructureFactor(                                                                                                                                                                   )))),// Hassium
    e_Meitnerium(    AddEnum(  "Meitnerium"   , Element(            109,      7,     0, "Mt"    , "Meitnerium"   ,     266,                                  0,              0,       1.70,            0,           0,                                     0,                0, ElectronConfiguration( util::GetUndefined< size_t>(),                  23,   2,   0,   0,   0,   2,   6,   0,   0,   2,   6,  10,   0,   2,   6,  10,  14,   2,   6,  10,  14,   2,   6,   7,   0,   2,   0,   0,   0 ),    0.921568627f,      0.000000000f,     0.149019608f, ElementStructureFactor(                                                                                                                                                                   )))) // Meitnerium
*/
}

utility::vector1< core::chemical::gasteiger::GasteigerAtomTypeDataCOP >
initialize_atoms ( core::chemical::ElementSetCOP ele_set) {
	using namespace core::chemical::gasteiger;
	utility::vector1< core::chemical::gasteiger::GasteigerAtomTypeDataCOP > out;

	out.push_back( new
        GasteigerAtomTypeData
        (
        "H_S",
          "H",               // Element type
          GasteigerAtomTypeData::Unhybridized,     // Hybridization
          0,                                          // # hybridized sigma orbitals in bonds
          0,                                          // # non-binding hybrid orbitals
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),    // pi orbitals in bonds
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),        // non-binding atomic orbitals
          13.60,
          0.75,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          0.387,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Li_S",
          "Li",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          5.39,
          0.82,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Li_P",
          "Li",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx1" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          3.54,
          0.56,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "BSP",
          "Be",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx01" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          9.92,
          3.18,
          5.96,
          0.11,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Be_PP",
          "Be",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx12" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          6.11,
          0.76,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Be_DiDi",
          "Be",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          8.58,
          0.99,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Be_DiPi",
          "Be",
          GasteigerAtomTypeData::SP,
          1,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx2" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          8.02,
          0.92,
          6.04,
          0.43,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Be_TrTr",
          "Be",
          GasteigerAtomTypeData::SP2,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          7.61,
          0.59,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Be_TrPi",
          "Be",
          GasteigerAtomTypeData::SP2,
          1,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          7.38,
          0.63,
          6.06,
          0.54,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Be_TeTe",
          "Be",
          GasteigerAtomTypeData::SP3,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          7.18,
          0.51,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_SPP",
          "B",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx012" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          14.91,
          5.70,
          8.42,
          0.32,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_PPP",
          "B",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx123" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          8.40,
          3.46,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_DiDiPi",
          "B",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          12.55,
          2.12,
          8.23,
          0.44,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_DiPiPi",
          "B",
          GasteigerAtomTypeData::SP,
          1,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          11.66,
          2.56,
          8.41,
          1.89,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_TrTrTr",
          "B",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          11.29,
          1.38,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_TrTrPi",
          "B",
          GasteigerAtomTypeData::SP2,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          10.97,
          1.87,
          8.33,
          1.42,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_TeTeTe",
          "B",
          GasteigerAtomTypeData::SP3,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          10.43,
          1.53,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "C_SPPP",
          "C",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0123" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          21.01,
          8.91,
          11.27,
          0.34,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "C_DiDiPiPi",
          "C",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          17.42,
          3.34,
          11.19,
          0.1,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          1.283,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "C_TrTrTrPi",
          "C",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          15.62,
          1.95,
          11.16,
          0.03,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          1.352, // = 1.352 if connected to >0 H's, 1.896 otherwise.  Assume for now that at least one bond is to an H
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "C_TeTeTeTe",
          "C",
          GasteigerAtomTypeData::SP3,
          4,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          14.61,
          1.34,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          1.061,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_S2PPP",
          "N",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx123" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          13.94,
          0.84,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_SP2PP",
          "N",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx023" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx1" ),
          26.92,
          14.05,
          14.42,
          2.54,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_Di2DiPiPi",
          "N",
          GasteigerAtomTypeData::SP,
          1,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          23.91,
          7.45,
          14.18,
          1.66,
          37.024,
          17.254,
          0.956,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_DiDiPi2Pi",
          "N",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx1" ),
          22.10,
          6.84,
          14.11,
          2.14,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          1.012, // N_Di2DiPiPi * N_TrTrTrPi2 / N_Tr2TrTrPi
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_Tr2TrTrPi",
          "N",
          GasteigerAtomTypeData::SP2,
          2,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          20.60,
          5.14,
          14.12,
          1.78,
          34.645,
          15.107,
          1.030,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_TrTrTrPi2",
          "N",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx1" ),
          19.72,
          4.92,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          1.090,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_Te2TeTeTe",
          "N",
          GasteigerAtomTypeData::SP3,
          3,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          18.93,
          4.15,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          33.313,
          14.153,
          0.964,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_S2P2PP",
          "O",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx12" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx03" ),
          17.28,
          2.01,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_SP2P2P",
          "O",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx01" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          36.07,
          18.44,
          18.53,
          3.40,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_Di2Di2PiPi",
          "O",
          GasteigerAtomTypeData::SP,
          0,
          2,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          17.28,
          2.01,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_Di2DiPi2Pi",
          "O",
          GasteigerAtomTypeData::SP,
          1,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx2" ),
          30.17,
          10.23,
          17.91,
          2.71,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          0.569, // similar to Tr2Tr2TrPi,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_DiDiPi2Pi2",
          "O",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          28.71,
          9.51,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          0.637, // similar to Te2Te2TeTe,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_Tr2Tr2TrPi",
          "O",
          GasteigerAtomTypeData::SP2,
          1,
          2,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          26.65,
          7.49,
          17.70,
          2.47,
          42.534,
          20.154,
          0.569,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_Tr2TrTrPi2",
          "O",
          GasteigerAtomTypeData::SP2,
          2,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          26.14,
          7.32,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          0.274,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_Te2Te2TeTe",
          "O",
          GasteigerAtomTypeData::SP3,
          2,
          2,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          24.39,
          6.11,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          40.358,
          18.708,
          0.637,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "F_S2P2P2P",
          "F",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx012" ),
          20.86,
          3.50,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          26.368,
          13.378,
          0.296,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "F_SP2P2P2",
          "F",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx123" ),
          38.24,
          24.37,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          26.368,
          13.378,
          0.296,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Na_S",
          "Na",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          5.14,
          0.47,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Na_P",
          "Na",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx1" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          3.04,
          0.09,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Mg_SP",
          "Mg",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx01" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          8.95,
          2.80,
          4.52,
          0.06,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Mg_PP",
          "Mg",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx12" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          5.65,
          0.01,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Mg_DiDi",
          "Mg",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          7.10,
          1.08,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Mg_DiPi",
          "Mg",
          GasteigerAtomTypeData::SP,
          1,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          7.30,
          0.78,
          5.09,
          0.03,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Mg_TrTr",
          "Mg",
          GasteigerAtomTypeData::SP2,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          6.54,
          0.52,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Mg_TrPi",
          "Mg",
          GasteigerAtomTypeData::SP2,
          1,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          6.75,
          0.38,
          5.27,
          0.02,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Mg_TeTe",
          "Mg",
          GasteigerAtomTypeData::SP3,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          6.28,
          0.32,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_SPP",
          "Al",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx012" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          12.27,
          4.92,
          6.47,
          1.37,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_PPP",
          "Al",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx123" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          6.50,
          4.89,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_DiDiPi",
          "Al",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          9.91,
          2.61,
          6.36,
          1.45,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_DiPiPi",
          "Al",
          GasteigerAtomTypeData::SP,
          1,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          9.39,
          3.66,
          6.49,
          3.13,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_TrTrTr",
          "Al",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          8.83,
          2.11,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_TrTrPi",
          "Al",
          GasteigerAtomTypeData::SP2,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          8.65,
          2.94,
          6.43,
          2.58,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_TeTeTe",
          "Al",
          GasteigerAtomTypeData::SP3,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          8.17,
          2.58,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Si_SPPP",
          "Si",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0123" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          17.31,
          6.94,
          9.19,
          2.82,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Si_DiDiPiPi",
          "Si",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          14.06,
          4.07,
          9.18,
          2.20,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Si_TrTrTrPi",
          "Si",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          12.61,
          3.20,
          9.17,
          2.00,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Si_TeTeTeTe",
          "Si",
          GasteigerAtomTypeData::SP3,
          4,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          11.82,
          2.78,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "P_S2PPP",
          "P",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx123" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          10.73,
          1.42,
          31.172,
          18.612,
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "P_SP2PP",
          "P",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx023" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx1" ),
          20.20,
          8.48,
          12.49,
          1.98,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "P_Di2DiPiPi",
          "P",
          GasteigerAtomTypeData::SP,
          1,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          17.53,
          4.95,
          11.61,
          1.68,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          1.525, // P_Te2TeTeTe * N_Di2DiPiPi/N_Te2TeTeTe,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "P_DiDiPi2Pi",
          "P",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx1" ),
          16.78,
          4.77,
          11.89,
          2.02,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "P_Tr2TrTrPi",
          "P",
          GasteigerAtomTypeData::SP2,
          2,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          15.59,
          3.74,
          11.64,
          1.80,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          1.643, // P_Te2TeTeTe * N_Tr2TrTrPi/N_Te2TeTeTe,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "P_TrTrTrPi2",
          "P",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx1" ),
          15.18,
          3.76,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          1.739, // P_Te2TeTeTe * N_TrTrTrPi2/N_Te2TeTeTe,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "P_Te2TeTeTe",
          "P",
          GasteigerAtomTypeData::SP3,
          3,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          14.57,
          3.24,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          24.041,
          12.095,
          1.538, // Miller 1990 contains contains a typo; the lewis-diagram is clearly of P_Te2TeTeTe, but they call it P_TeTeTeTe
                // Presumably they meant Te2 since TeTeTeTe would be charged,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_S2P2PP",
          "S",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx12" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx03" ),
          12.39,
          2.38,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          22.977,
          11.053,
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_SP2P2P",
          "S",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx01" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          20.08,
          11.54,
          13.32,
          3.50,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_Di2Di2PiPi",
          "S",
          GasteigerAtomTypeData::SP,
          0,
          2,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          12.39,
          2.38,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_Di2DiPi2Pi",
          "S",
          GasteigerAtomTypeData::SP,
          1,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx2" ),
          17.78,
          6.96,
          12.86,
          2.94,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          3.729, // similar to Tr2Tr2TrPi,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_DiDiPi2Pi2",
          "S",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          17.42,
          6.80,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_Tr2Tr2TrPi",
          "S",
          GasteigerAtomTypeData::SP2,
          1,
          2,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          16.33,
          5.43,
          12.70,
          2.76,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          3.729,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_Tr2TrTrPi2",
          "S",
          GasteigerAtomTypeData::SP2,
          2,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          16.27,
          5.49,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          2.700,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_Te2Te2TeTe",
          "S",
          GasteigerAtomTypeData::SP3,
          2,
          2,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          15.50,
          4.77,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          27.728,
          13.638,
          3.000,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Cl_S2P2P2P",
          "Cl",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx012" ),
          15.03,
          3.73,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          26.368,
          13.378,
          2.315,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Cl_SP2P2P2",
          "Cl",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx123" ),
          24.02,
          14.45,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          26.368,
          13.378,
          2.315,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "K_S",
          "K",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          4.341,
          1.95,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "K_P",
          "K",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx1" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          2.7,
          1.195,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
	ele_set
        )
      );
    // cations
      out.push_back( new
        GasteigerAtomTypeData
        (
        "BS",
          "Be",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          18.21,
          9.32,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Be_P",
          "Be",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx1" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          14.25,
          5.32,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_SP",
          "B",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx01" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          25.40,
          14.05,
          19.40,
          7.38,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_PP",
          "B",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx12" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          18.91,
          7.37,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_DiDi",
          "B",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          23.48,
          9.64,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_DiPi",
          "B",
          GasteigerAtomTypeData::SP,
          1,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx2" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          22.16,
          8.94,
          19.16,
          7.37,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_TrTr",
          "B",
          GasteigerAtomTypeData::SP2,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          21.72,
          8.33,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_TrPi",
          "B",
          GasteigerAtomTypeData::SP2,
          1,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          21.08,
          8.02,
          19.08,
          7.37,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_TeTe",
          "B",
          GasteigerAtomTypeData::SP3,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          20.93,
          7.88,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "C_SPP",
          "C",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx012" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          33.03,
          19.42,
          23.93,
          9.91,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "C_PPP",
          "C",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx123" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          23.29,
          11.65,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "C_DiDiPi",
          "C",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx2" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          29.85,
          13.29,
          23.86,
          9.83,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "C_DiPiPi",
          "C",
          GasteigerAtomTypeData::SP,
          1,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          28.16,
          12.96,
          23.61,
          10.78,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "C_TrTrTr",
          "C",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          28.14,
          11.83,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "C_TrTrPi",
          "C",
          GasteigerAtomTypeData::SP2,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          27.36,
          11.91,
          23.68,
          10.45,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "C_TeTeTe",
          "C",
          GasteigerAtomTypeData::SP3,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          26.71,
          11.37,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_SPPP",
          "N",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0123" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          41.84,
          25.59,
          28.69,
          12.48,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_DiDiPiPi",
          "N",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          37.00,
          17.24,
          28.70,
          12.06,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_TrTrTrPi",
          "N",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          34.62,
          15.09,
          28.71,
          11.96,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_TeTeTeTe",
          "N",
          GasteigerAtomTypeData::SP3,
          4,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          33.29,
          14.14,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_S2PPP",
          "O",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx123" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          34.15,
          14.61,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_SP2PP",
          "O",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx012" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          51.41,
          32.29,
          34.22,
          15.86,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_Di2DiPiPi",
          "O",
          GasteigerAtomTypeData::SP,
          1,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          46.80,
          23.45,
          34.19,
          15.24,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_DiDiPi2Pi",
          "O",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx2" ),
          44.56,
          22.34,
          33.95,
          15.53,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_Tr2TrTrPi",
          "O",
          GasteigerAtomTypeData::SP2,
          2,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          42.49,
          20.15,
          34.08,
          15.30,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_TrTrTrPi2",
          "O",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          41.39,
          19.64,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_Te2TeTeTe",
          "O",
          GasteigerAtomTypeData::SP3,
          3,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          40.31,
          18.70,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Mg_S",
          "Mg",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          15.03,
          7.64,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Mg_P",
          "Mg",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx1" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          10.60,
          4.67,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_SP",
          "Al",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx01" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          20.15,
          11.32,
          13.48,
          5.99,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_PP",
          "Al",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx12" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          14.34,
          6.03,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_DiDi",
          "Al",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          17.47,
          8.00,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_DiPi",
          "Al",
          GasteigerAtomTypeData::SP,
          1,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx2" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          17.25,
          7.59,
          13.92,
          6.00,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_TrTr",
          "Al",
          GasteigerAtomTypeData::SP2,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          16.28,
          7.01,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_TrPi",
          "Al",
          GasteigerAtomTypeData::SP2,
          1,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          16.28,
          6.74,
          14.06,
          5.92,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Al_TeTe",
          "Al",
          GasteigerAtomTypeData::SP3,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          15.75,
          6.64,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Si_SPP",
          "Si",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx012" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          24.68,
          14.93,
          16.56,
          8.61,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Si_PPP",
          "Si",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx123" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          16.56,
          11.42,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Si_DiDiPi",
          "Si",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx2" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          21.43,
          10.95,
          16.50,
          8.60,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Si_DiPiPi",
          "Si",
          GasteigerAtomTypeData::SP,
          1,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          20.62,
          11.56,
          16.55,
          10.02,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Si_TrTrTr",
          "Si",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          19.96,
          9.99,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Si_TrTrPi",
          "Si",
          GasteigerAtomTypeData::SP2,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          19.62,
          10.57,
          16.53,
          9.54,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Si_TeTeTe",
          "Si",
          GasteigerAtomTypeData::SP3,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          18.97,
          10.08,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "P_SPPP",
          "P",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0123" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          31.24,
          18.61,
          20.72,
          11.55,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "P_DiDiPiPi",
          "P",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          27.01,
          14.05,
          20.69,
          10.96,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "P_TrTrTrPi",
          "P",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          25.14,
          12.72,
          20.68,
          10.76,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "P_TeTeTeTe",
          "P",
          GasteigerAtomTypeData::SP3,
          4,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          24.10,
          12.09,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_S2PPP",
          "S",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx123" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          22.91,
          11.05,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_SP2PP",
          "S",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx012" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          35.18,
          21.13,
          24.49,
          11.98,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_Di2DiPiPi",
          "S",
          GasteigerAtomTypeData::SP,
          1,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          31.57,
          16.09,
          23.70,
          11.51,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_DiDiPi2Pi",
          "S",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx2" ),
          30.61,
          15.78,
          24.00,
          11.92,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_Tr2TrTrPi",
          "S",
          GasteigerAtomTypeData::SP2,
          2,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          28.99,
          14.38,
          23.74,
          11.65,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_TrTrTrPi2",
          "S",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          28.51,
          14.33,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_Te2TeTeTe",
          "S",
          GasteigerAtomTypeData::SP3,
          3,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          27.65,
          13.64,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "O_Te2Te2Te2Te",
          "O",
          GasteigerAtomTypeData::SP3,
          1,
          3,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          6.11,
          0.0,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "P_TeTeTeTePi",
          "P",
          GasteigerAtomTypeData::SP3,
          4,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          17.704,
          5.694,
          5.385,
          -0.015,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          1.523, // P_Te2TeTeTe * N_TrTrTrPi2/N_Tr2TrTrPi * N_Te2TeTeTe / N_Tr2TrTrPi,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_TeTeTeTePiPi",
          "S",
          GasteigerAtomTypeData::SP3,
          4,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          20.59,
          6.69,
          5.39,
          -2.85,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Br_SP2P2P2",
          "Br",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx123" ),
          22.081,
          14.315,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          26.368,
          13.378,
          3.013,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Br_S2P2P2P",
          "Br",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx012" ),
          13.108,
          3.516,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          26.368,
          13.378,
          3.013,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "I_SP2P2P2",
          "I",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx123" ),
          18.01,
          13.23,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          26.368,
          13.378,
          5.415,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "I_S2P2P2P",
          "I",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx012" ),
          12.677,
          3.375,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          26.368,
          13.378,
          5.415,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Se_Te2Te2TeTe",
          "Se",
          GasteigerAtomTypeData::SP3,
          2,
          2,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          20.908,
          10.469,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "S_Te2TeTeTePi",
          "S",
          GasteigerAtomTypeData::SP3,
          3,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          18.136,
          5.708,
          2.283,
          -4.393,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_TrTrTrPi2Pi",
          "N",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx2" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Sn_TeTeTeTe",
          "Sn",
          GasteigerAtomTypeData::SP3,
          4,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          10.4,
          5.39,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Ge_TeTeTeTe",
          "Ge",
          GasteigerAtomTypeData::SP3,
          4,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          11.48,
          4.66,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_TeTeTeTe",
          "B",
          GasteigerAtomTypeData::SP3,
          4,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          1.53,
          0.0,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "B_TrTrTrPi",
          "B",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          1.87,
          0.0,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Cl_S2P2P2P2",
          "Cl",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0123" ),
          14.45,
          0.0,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Se_Di2DiPi2Pi",
          "Se",
          GasteigerAtomTypeData::SP,
          1,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx2" ),
          17.29,
          6.44,
          13.06,
          2.28,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Te_Te2Te2TeTe",
          "Te",
          GasteigerAtomTypeData::SP3,
          2,
          2,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          15.11,
          4.20,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "I_S2P2P2P2",
          "I",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx0123" ),
          13.38,
          0.0,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          26.368,
          13.378,
          5.573,
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "As_Te2TeTeTe",
          "As",
          GasteigerAtomTypeData::SP3,
          3,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          12.80,
          3.81,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_TrTrTrPiPi",
          "N",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "P_TrTrTrPiPi",
          "P",
          GasteigerAtomTypeData::SP2,
          3,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_TeTeTeTePi",
          "N",
          GasteigerAtomTypeData::SP3,
          4,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_DiDiPi2Pi2",
          "N",
          GasteigerAtomTypeData::SP,
          2,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx23" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_Di2DiPi2Pi",
          "N",
          GasteigerAtomTypeData::SP,
          1,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx2" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_Tr2TrTrPi2",
          "N",
          GasteigerAtomTypeData::SP2,
          2,
          1,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx3" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "N_Te2Te2TeTe",
          "N",
          GasteigerAtomTypeData::SP3,
          2,
          2,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "H_",
          "H",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          13.60 * 2.0,
          13.60,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Li_",
          "Li",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          5.39 * 2.0,
          5.39,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "Na_",
          "Na",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          5.14 * 2.0,
          5.14,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
        ele_set
      ));
      out.push_back( new
        GasteigerAtomTypeData
        (
        "K_",
          "K",
          GasteigerAtomTypeData::Unhybridized,
          0,
          0,
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          parse_enum_set< core::chemical::gasteiger::GasteigerAtomTypeData::AtomicOrbitalTypes > ( "xx" ),
          4.341 * 2.0,
          4.341,
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
          numeric::get_undefined_real(),
	ele_set
        )
      );

    return out;

/* Copy paste from BCL

      H_S
      (
        AddEnum
        (
          "H_S",
          AtomTypeData
          (
            "Li",               // Element type
            GetHybridOrbitalTypes().Unhybridized,     // Hybridization
            0,                                          // # hybridized sigma orbitals in bonds
            0,                                          // # non-binding hybrid orbitals
            storage::Set< AtomicOrbitalTypesEnum>( S),    // pi orbitals in bonds
            storage::Set< AtomicOrbitalTypesEnum>(),        // non-binding atomic orbitals
            13.60,
            0.75,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            0.387
          )
        )
      ),
      Li_S
      (
        AddEnum
        (
          "Li_S",
          AtomTypeData
          (
            GetElementTypes().e_Lithium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( S),
            storage::Set< AtomicOrbitalTypesEnum>(),
            5.39,
            0.82,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Li_P
      (
        AddEnum
        (
          "Li_P",
          AtomTypeData
          (
            GetElementTypes().e_Lithium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            3.54,
            0.56,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      BSP
      (
        AddEnum
        (
          "BSP",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            9.92,
            3.18,
            5.96,
            0.11,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_PP
      (
        AddEnum
        (
          "Be_PP",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            6.11,
            0.76,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_DiDi
      (
        AddEnum
        (
          "Be_DiDi",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.58,
            0.99,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_DiPi
      (
        AddEnum
        (
          "Be_DiPi",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.02,
            0.92,
            6.04,
            0.43,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_TrTr
      (
        AddEnum
        (
          "Be_TrTr",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            7.61,
            0.59,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_TrPi
      (
        AddEnum
        (
          "Be_TrPi",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().SP2,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            7.38,
            0.63,
            6.06,
            0.54,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_TeTe
      (
        AddEnum
        (
          "Be_TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().SP3,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            7.18,
            0.51,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_SPP
      (
        AddEnum
        (
          "B_SPP",
          AtomTypeData
          (
            "C",
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            14.91,
            5.70,
            8.42,
            0.32,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_PPP
      (
        AddEnum
        (
          "B_PPP",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.40,
            3.46,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_DiDiPi
      (
        AddEnum
        (
          "B_DiDiPi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            12.55,
            2.12,
            8.23,
            0.44,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_DiPiPi
      (
        AddEnum
        (
          "B_DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            11.66,
            2.56,
            8.41,
            1.89,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TrTrTr
      (
        AddEnum
        (
          "B_TrTrTr",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            11.29,
            1.38,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TrTrPi
      (
        AddEnum
        (
          "B_TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            10.97,
            1.87,
            8.33,
            1.42,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TeTeTe
      (
        AddEnum
        (
          "B_TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().SP3,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            10.43,
            1.53,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_SPPP
      (
        AddEnum
        (
          "C_SPPP",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            21.01,
            8.91,
            11.27,
            0.34,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_DiDiPiPi
      (
        AddEnum
        (
          "C_DiDiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            17.42,
            3.34,
            11.19,
            0.1,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.283
          )
        )
      ),
      C_TrTrTrPi
      (
        AddEnum
        (
          "C_TrTrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            15.62,
            1.95,
            11.16,
            0.03,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.352 // = 1.352 if connected to >0 H's, 1.896 otherwise.  Assume for now that at least one bond is to an H
          )
        )
      ),
      C_TeTeTeTe
      (
        AddEnum
        (
          "C_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            14.61,
            1.34,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.061
          )
        )
      ),
      N_S2PPP
      (
        AddEnum
        (
          "N_S2PPP",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( S),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            13.94,
            0.84,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_SP2PP
      (
        AddEnum
        (
          "N_SP2PP",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            26.92,
            14.05,
            14.42,
            2.54,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_Di2DiPiPi
      (
        AddEnum
        (
          "N_Di2DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            23.91,
            7.45,
            14.18,
            1.66,
            37.024,
            17.254,
            0.956
          )
        )
      ),
      N_DiDiPi2Pi
      (
        AddEnum
        (
          "N_DiDiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            22.10,
            6.84,
            14.11,
            2.14,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.012 // N_Di2DiPiPi * N_TrTrTrPi2 / N_Tr2TrTrPi
          )
        )
      ),
      N_Tr2TrTrPi
      (
        AddEnum
        (
          "N_Tr2TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            20.60,
            5.14,
            14.12,
            1.78,
            34.645,
            15.107,
            1.030
          )
        )
      ),
      N_TrTrTrPi2
      (
        AddEnum
        (
          "N_TrTrTrPi2",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            19.72,
            4.92,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.090
          )
        )
      ),
      N_Te2TeTeTe
      (
        AddEnum
        (
          "N_Te2TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP3,
            3,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            18.93,
            4.15,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            33.313,
            14.153,
            0.964
          )
        )
      ),
      O_S2P2PP
      (
        AddEnum
        (
          "O_S2P2PP",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Pz),
            17.28,
            2.01,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_SP2P2P
      (
        AddEnum
        (
          "O_SP2P2P",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            36.07,
            18.44,
            18.53,
            3.40,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_Di2Di2PiPi
      (
        AddEnum
        (
          "O_Di2Di2PiPi",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().SP,
            0,
            2,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            17.28,
            2.01,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_Di2DiPi2Pi
      (
        AddEnum
        (
          "O_Di2DiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            30.17,
            10.23,
            17.91,
            2.71,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            0.569 // similar to Tr2Tr2TrPi
          )
        )
      ),
      O_DiDiPi2Pi2
      (
        AddEnum
        (
          "O_DiDiPi2Pi2",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            28.71,
            9.51,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            0.637 // similar to Te2Te2TeTe
          )
        )
      ),
      O_Tr2Tr2TrPi
      (
        AddEnum
        (
          "O_Tr2Tr2TrPi",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().SP2,
            1,
            2,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            26.65,
            7.49,
            17.70,
            2.47,
            42.534,
            20.154,
            0.569
          )
        )
      ),
      O_Tr2TrTrPi2
      (
        AddEnum
        (
          "O_Tr2TrTrPi2",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            26.14,
            7.32,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            0.274
          )
        )
      ),
      O_Te2Te2TeTe
      (
        AddEnum
        (
          "O_Te2Te2TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().SP3,
            2,
            2,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            24.39,
            6.11,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            40.358,
            18.708,
            0.637
          )
        )
      ),
      F_S2P2P2P
      (
        AddEnum
        (
          "F_S2P2P2P",
          AtomTypeData
          (
            GetElementTypes().e_Fluorine,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py),
            20.86,
            3.50,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            0.296
          )
        )
      ),
      F_SP2P2P2
      (
        AddEnum
        (
          "F_SP2P2P2",
          AtomTypeData
          (
            GetElementTypes().e_Fluorine,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( S),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            38.24,
            24.37,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            0.296
          )
        )
      ),
      Na_S
      (
        AddEnum
        (
          "Na_S",
          AtomTypeData
          (
            GetElementTypes().Sodium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( S),
            storage::Set< AtomicOrbitalTypesEnum>(),
            5.14,
            0.47,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Na_P
      (
        AddEnum
        (
          "Na_P",
          AtomTypeData
          (
            GetElementTypes().Sodium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            3.04,
            0.09,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_SP
      (
        AddEnum
        (
          "Mg_SP",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.95,
            2.80,
            4.52,
            0.06,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_PP
      (
        AddEnum
        (
          "Mg_PP",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            5.65,
            0.01,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_DiDi
      (
        AddEnum
        (
          "Mg_DiDi",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            7.10,
            1.08,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_DiPi
      (
        AddEnum
        (
          "Mg_DiPi",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            7.30,
            0.78,
            5.09,
            0.03,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_TrTr
      (
        AddEnum
        (
          "Mg_TrTr",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            6.54,
            0.52,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_TrPi
      (
        AddEnum
        (
          "Mg_TrPi",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().SP2,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            6.75,
            0.38,
            5.27,
            0.02,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_TeTe
      (
        AddEnum
        (
          "Mg_TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().SP3,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            6.28,
            0.32,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_SPP
      (
        AddEnum
        (
          "Al_SPP",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            12.27,
            4.92,
            6.47,
            1.37,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_PPP
      (
        AddEnum
        (
          "Al_PPP",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            6.50,
            4.89,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_DiDiPi
      (
        AddEnum
        (
          "Al_DiDiPi",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            9.91,
            2.61,
            6.36,
            1.45,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_DiPiPi
      (
        AddEnum
        (
          "Al_DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            9.39,
            3.66,
            6.49,
            3.13,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_TrTrTr
      (
        AddEnum
        (
          "Al_TrTrTr",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.83,
            2.11,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_TrTrPi
      (
        AddEnum
        (
          "Al_TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.65,
            2.94,
            6.43,
            2.58,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_TeTeTe
      (
        AddEnum
        (
          "Al_TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().SP3,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            8.17,
            2.58,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_SPPP
      (
        AddEnum
        (
          "Si_SPPP",
          AtomTypeData
          (
            GetElementTypes().Silicon,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            17.31,
            6.94,
            9.19,
            2.82,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_DiDiPiPi
      (
        AddEnum
        (
          "Si_DiDiPiPi",
          AtomTypeData
          (
            GetElementTypes().Silicon,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            14.06,
            4.07,
            9.18,
            2.20,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_TrTrTrPi
      (
        AddEnum
        (
          "Si_TrTrTrPi",
          AtomTypeData
          (
            GetElementTypes().Silicon,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            12.61,
            3.20,
            9.17,
            2.00,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_TeTeTeTe
      (
        AddEnum
        (
          "Si_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().Silicon,
            GetHybridOrbitalTypes().SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            11.82,
            2.78,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_S2PPP
      (
        AddEnum
        (
          "P_S2PPP",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( S),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            10.73,
            1.42,
            31.172,
            18.612,
            util::GetUndefined< double>()
          )
        )
      ),
      P_SP2PP
      (
        AddEnum
        (
          "P_SP2PP",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            20.20,
            8.48,
            12.49,
            1.98,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_Di2DiPiPi
      (
        AddEnum
        (
          "P_Di2DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            17.53,
            4.95,
            11.61,
            1.68,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.525 // P_Te2TeTeTe * N_Di2DiPiPi/N_Te2TeTeTe
          )
        )
      ),
      P_DiDiPi2Pi
      (
        AddEnum
        (
          "P_DiDiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            16.78,
            4.77,
            11.89,
            2.02,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_Tr2TrTrPi
      (
        AddEnum
        (
          "P_Tr2TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            15.59,
            3.74,
            11.64,
            1.80,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.643 // P_Te2TeTeTe * N_Tr2TrTrPi/N_Te2TeTeTe
          )
        )
      ),
      P_TrTrTrPi2
      (
        AddEnum
        (
          "P_TrTrTrPi2",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            15.18,
            3.76,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.739 // P_Te2TeTeTe * N_TrTrTrPi2/N_Te2TeTeTe
          )
        )
      ),
      P_Te2TeTeTe
      (
        AddEnum
        (
          "P_Te2TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().SP3,
            3,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            14.57,
            3.24,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            24.041,
            12.095,
            1.538 // Miller 1990 contains contains a typo; the lewis-diagram is clearly of P_Te2TeTeTe, but they call it P_TeTeTeTe
                  // Presumably they meant Te2 since TeTeTeTe would be charged
          )
        )
      ),
      S_S2P2PP
      (
        AddEnum
        (
          "S_S2P2PP",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Pz),
            12.39,
            2.38,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            22.977,
            11.053,
            util::GetUndefined< double>()
          )
        )
      ),
      S_SP2P2P
      (
        AddEnum
        (
          "S_SP2P2P",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            20.08,
            11.54,
            13.32,
            3.50,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Di2Di2PiPi
      (
        AddEnum
        (
          "S_Di2Di2PiPi",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().SP,
            0,
            2,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            12.39,
            2.38,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Di2DiPi2Pi
      (
        AddEnum
        (
          "S_Di2DiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            17.78,
            6.96,
            12.86,
            2.94,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            3.729 // similar to Tr2Tr2TrPi
          )
        )
      ),
      S_DiDiPi2Pi2
      (
        AddEnum
        (
          "S_DiDiPi2Pi2",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            17.42,
            6.80,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Tr2Tr2TrPi
      (
        AddEnum
        (
          "S_Tr2Tr2TrPi",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().SP2,
            1,
            2,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            16.33,
            5.43,
            12.70,
            2.76,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            3.729
          )
        )
      ),
      S_Tr2TrTrPi2
      (
        AddEnum
        (
          "S_Tr2TrTrPi2",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            16.27,
            5.49,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            2.700
          )
        )
      ),
      S_Te2Te2TeTe
      (
        AddEnum
        (
          "S_Te2Te2TeTe",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().SP3,
            2,
            2,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            15.50,
            4.77,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            27.728,
            13.638,
            3.000
          )
        )
      ),
      Cl_S2P2P2P
      (
        AddEnum
        (
          "Cl_S2P2P2P",
          AtomTypeData
          (
            GetElementTypes().e_Chlorine,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py),
            15.03,
            3.73,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            2.315
          )
        )
      ),
      Cl_SP2P2P2
      (
        AddEnum
        (
          "Cl_SP2P2P2",
          AtomTypeData
          (
            GetElementTypes().e_Chlorine,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( S),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            24.02,
            14.45,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            2.315
          )
        )
      ),
      K_S
      (
        AddEnum
        (
          "K_S",
          AtomTypeData
          (
            GetElementTypes().e_Potassium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( S),
            storage::Set< AtomicOrbitalTypesEnum>(),
            4.341,
            1.95,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      K_P
      (
        AddEnum
        (
          "K_P",
          AtomTypeData
          (
            GetElementTypes().e_Potassium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            2.7,
            1.195,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      // cations
      BS
      (
        AddEnum
        (
          "BS",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( S),
            storage::Set< AtomicOrbitalTypesEnum>(),
            18.21,
            9.32,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Be_P
      (
        AddEnum
        (
          "Be_P",
          AtomTypeData
          (
            GetElementTypes().e_Beryllium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            14.25,
            5.32,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_SP
      (
        AddEnum
        (
          "B_SP",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            25.40,
            14.05,
            19.40,
            7.38,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_PP
      (
        AddEnum
        (
          "B_PP",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            18.91,
            7.37,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_DiDi
      (
        AddEnum
        (
          "B_DiDi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            23.48,
            9.64,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_DiPi
      (
        AddEnum
        (
          "B_DiPi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            22.16,
            8.94,
            19.16,
            7.37,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TrTr
      (
        AddEnum
        (
          "B_TrTr",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            21.72,
            8.33,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TrPi
      (
        AddEnum
        (
          "B_TrPi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().SP2,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            21.08,
            8.02,
            19.08,
            7.37,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TeTe
      (
        AddEnum
        (
          "B_TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().SP3,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            20.93,
            7.88,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_SPP
      (
        AddEnum
        (
          "C_SPP",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            33.03,
            19.42,
            23.93,
            9.91,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_PPP
      (
        AddEnum
        (
          "C_PPP",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            23.29,
            11.65,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_DiDiPi
      (
        AddEnum
        (
          "C_DiDiPi",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            29.85,
            13.29,
            23.86,
            9.83,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_DiPiPi
      (
        AddEnum
        (
          "C_DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            28.16,
            12.96,
            23.61,
            10.78,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_TrTrTr
      (
        AddEnum
        (
          "C_TrTrTr",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            28.14,
            11.83,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_TrTrPi
      (
        AddEnum
        (
          "C_TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            27.36,
            11.91,
            23.68,
            10.45,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      C_TeTeTe
      (
        AddEnum
        (
          "C_TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Carbon,
            GetHybridOrbitalTypes().SP3,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            26.71,
            11.37,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_SPPP
      (
        AddEnum
        (
          "N_SPPP",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            41.84,
            25.59,
            28.69,
            12.48,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_DiDiPiPi
      (
        AddEnum
        (
          "N_DiDiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            37.00,
            17.24,
            28.70,
            12.06,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_TrTrTrPi
      (
        AddEnum
        (
          "N_TrTrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            34.62,
            15.09,
            28.71,
            11.96,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_TeTeTeTe
      (
        AddEnum
        (
          "N_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            33.29,
            14.14,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_S2PPP
      (
        AddEnum
        (
          "O_S2PPP",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( S),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            34.15,
            14.61,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_SP2PP
      (
        AddEnum
        (
          "O_SP2PP",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            51.41,
            32.29,
            34.22,
            15.86,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_Di2DiPiPi
      (
        AddEnum
        (
          "O_Di2DiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            46.80,
            23.45,
            34.19,
            15.24,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_DiDiPi2Pi
      (
        AddEnum
        (
          "O_DiDiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            44.56,
            22.34,
            33.95,
            15.53,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_Tr2TrTrPi
      (
        AddEnum
        (
          "O_Tr2TrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            42.49,
            20.15,
            34.08,
            15.30,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_TrTrTrPi2
      (
        AddEnum
        (
          "O_TrTrTrPi2",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            41.39,
            19.64,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_Te2TeTeTe
      (
        AddEnum
        (
          "O_Te2TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().SP3,
            3,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            40.31,
            18.70,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_S
      (
        AddEnum
        (
          "Mg_S",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( S),
            storage::Set< AtomicOrbitalTypesEnum>(),
            15.03,
            7.64,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Mg_P
      (
        AddEnum
        (
          "Mg_P",
          AtomTypeData
          (
            GetElementTypes().e_Magnesium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            10.60,
            4.67,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_SP
      (
        AddEnum
        (
          "Al_SP",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px),
            storage::Set< AtomicOrbitalTypesEnum>(),
            20.15,
            11.32,
            13.48,
            5.99,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_PP
      (
        AddEnum
        (
          "Al_PP",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            14.34,
            6.03,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_DiDi
      (
        AddEnum
        (
          "Al_DiDi",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            17.47,
            8.00,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_DiPi
      (
        AddEnum
        (
          "Al_DiPi",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            17.25,
            7.59,
            13.92,
            6.00,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_TrTr
      (
        AddEnum
        (
          "Al_TrTr",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            16.28,
            7.01,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_TrPi
      (
        AddEnum
        (
          "Al_TrPi",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().SP2,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            16.28,
            6.74,
            14.06,
            5.92,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Al_TeTe
      (
        AddEnum
        (
          "Al_TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Aluminum,
            GetHybridOrbitalTypes().SP3,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            15.75,
            6.64,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_SPP
      (
        AddEnum
        (
          "Si_SPP",
          AtomTypeData
          (
            GetElementTypes().Silicon,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            24.68,
            14.93,
            16.56,
            8.61,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_PPP
      (
        AddEnum
        (
          "Si_PPP",
          AtomTypeData
          (
            GetElementTypes().Silicon,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            16.56,
            11.42,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_DiDiPi
      (
        AddEnum
        (
          "Si_DiDiPi",
          AtomTypeData
          (
            GetElementTypes().Silicon,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            storage::Set< AtomicOrbitalTypesEnum>(),
            21.43,
            10.95,
            16.50,
            8.60,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_DiPiPi
      (
        AddEnum
        (
          "Si_DiPiPi",
          AtomTypeData
          (
            GetElementTypes().Silicon,
            GetHybridOrbitalTypes().SP,
            1,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            20.62,
            11.56,
            16.55,
            10.02,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_TrTrTr
      (
        AddEnum
        (
          "Si_TrTrTr",
          AtomTypeData
          (
            GetElementTypes().Silicon,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            19.96,
            9.99,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_TrTrPi
      (
        AddEnum
        (
          "Si_TrTrPi",
          AtomTypeData
          (
            GetElementTypes().Silicon,
            GetHybridOrbitalTypes().SP2,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            19.62,
            10.57,
            16.53,
            9.54,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Si_TeTeTe
      (
        AddEnum
        (
          "Si_TeTeTe",
          AtomTypeData
          (
            GetElementTypes().Silicon,
            GetHybridOrbitalTypes().SP3,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            18.97,
            10.08,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_SPPP
      (
        AddEnum
        (
          "P_SPPP",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            31.24,
            18.61,
            20.72,
            11.55,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_DiDiPiPi
      (
        AddEnum
        (
          "P_DiDiPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            27.01,
            14.05,
            20.69,
            10.96,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_TrTrTrPi
      (
        AddEnum
        (
          "P_TrTrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            25.14,
            12.72,
            20.68,
            10.76,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_TeTeTeTe
      (
        AddEnum
        (
          "P_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            24.10,
            12.09,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_S2PPP
      (
        AddEnum
        (
          "S_S2PPP",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( S),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            22.91,
            11.05,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_SP2PP
      (
        AddEnum
        (
          "S_SP2PP",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            35.18,
            21.13,
            24.49,
            11.98,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Di2DiPiPi
      (
        AddEnum
        (
          "S_Di2DiPiPi",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            31.57,
            16.09,
            23.70,
            11.51,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_DiDiPi2Pi
      (
        AddEnum
        (
          "S_DiDiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            30.61,
            15.78,
            24.00,
            11.92,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Tr2TrTrPi
      (
        AddEnum
        (
          "S_Tr2TrTrPi",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            28.99,
            14.38,
            23.74,
            11.65,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_TrTrTrPi2
      (
        AddEnum
        (
          "S_TrTrTrPi2",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            28.51,
            14.33,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Te2TeTeTe
      (
        AddEnum
        (
          "S_Te2TeTeTe",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().SP3,
            3,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            27.65,
            13.64,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      O_Te2Te2Te2Te
      (
        AddEnum
        (
          "O_Te2Te2Te2Te",
          AtomTypeData
          (
            GetElementTypes().e_Oxygen,
            GetHybridOrbitalTypes().SP3,
            1,
            3,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            6.11,
            0.0,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_TeTeTeTePi
      (
        AddEnum
        (
          "P_TeTeTeTePi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            17.704,
            5.694,
            5.385,
            -0.015,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            1.523 // P_Te2TeTeTe * N_TrTrTrPi2/N_Tr2TrTrPi * N_Te2TeTeTe / N_Tr2TrTrPi
          )
        )
      ),
      S_TeTeTeTePiPi
      (
        AddEnum
        (
          "S_TeTeTeTePiPi",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            20.59,
            6.69,
            5.39,
            -2.85,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Br_SP2P2P2
      (
        AddEnum
        (
          "Br_SP2P2P2",
          AtomTypeData
          (
            GetElementTypes().e_Bromine,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( S),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            22.081,
            14.315,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            3.013
          )
        )
      ),
      Br_S2P2P2P
      (
        AddEnum
        (
          "Br_S2P2P2P",
          AtomTypeData
          (
            GetElementTypes().e_Bromine,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py),
            13.108,
            3.516,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            3.013
          )
        )
      ),
      I_SP2P2P2
      (
        AddEnum
        (
          "I_SP2P2P2",
          AtomTypeData
          (
            GetElementTypes().e_Iodine,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( S),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Px, e_Py, e_Pz),
            18.01,
            13.23,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            5.415
          )
        )
      ),
      I_S2P2P2P
      (
        AddEnum
        (
          "I_S2P2P2P",
          AtomTypeData
          (
            GetElementTypes().e_Iodine,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py),
            12.677,
            3.375,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            5.415
          )
        )
      ),
      Se_Te2Te2TeTe
      (
        AddEnum
        (
          "Se_Te2Te2TeTe",
          AtomTypeData
          (
            GetElementTypes().Selenium,
            GetHybridOrbitalTypes().SP3,
            2,
            2,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            20.908,
            10.469,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      S_Te2TeTeTePi
      (
        AddEnum
        (
          "S_Te2TeTeTePi",
          AtomTypeData
          (
            GetElementTypes().Sulfur,
            GetHybridOrbitalTypes().SP3,
            3,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            18.136,
            5.708,
            2.283,
            -4.393,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_TrTrTrPi2Pi
      (
        AddEnum
        (
          "N_TrTrTrPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Sn_TeTeTeTe
      (
        AddEnum
        (
          "Sn_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Tin,
            GetHybridOrbitalTypes().SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            10.4,
            5.39,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Ge_TeTeTeTe
      (
        AddEnum
        (
          "Ge_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Germanium,
            GetHybridOrbitalTypes().SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            11.48,
            4.66,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TeTeTeTe
      (
        AddEnum
        (
          "B_TeTeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            1.53,
            0.0,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      B_TrTrTrPi
      (
        AddEnum
        (
          "B_TrTrTrPi",
          AtomTypeData
          (
            GetElementTypes().e_Boron,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            1.87,
            0.0,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Cl_S2P2P2P2
      (
        AddEnum
        (
          "Cl_S2P2P2P2",
          AtomTypeData
          (
            GetElementTypes().e_Chlorine,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py, e_Pz),
            14.45,
            0.0,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Se_Di2DiPi2Pi
      (
        AddEnum
        (
          "Se_Di2DiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().Selenium,
            GetHybridOrbitalTypes().SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            17.29,
            6.44,
            13.06,
            2.28,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Te_Te2Te2TeTe
      (
        AddEnum
        (
          "Te_Te2Te2TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Tellurium,
            GetHybridOrbitalTypes().SP3,
            2,
            2,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            15.11,
            4.20,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      I_S2P2P2P2
      (
        AddEnum
        (
          "I_S2P2P2P2",
          AtomTypeData
          (
            GetElementTypes().e_Iodine,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>::Create( S, e_Px, e_Py, e_Pz),
            13.38,
            0.0,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            26.368,
            13.378,
            5.573
          )
        )
      ),
      As_Te2TeTeTe
      (
        AddEnum
        (
          "As_Te2TeTeTe",
          AtomTypeData
          (
            GetElementTypes().e_Arsenic,
            GetHybridOrbitalTypes().SP3,
            3,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            12.80,
            3.81,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_TrTrTrPiPi
      (
        AddEnum
        (
          "N_TrTrTrPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      P_TrTrTrPiPi
      (
        AddEnum
        (
          "P_TrTrTrPiPi",
          AtomTypeData
          (
            GetElementTypes().e_Phosphorus,
            GetHybridOrbitalTypes().SP2,
            3,
            0,
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_TeTeTeTePi
      (
        AddEnum
        (
          "N_TeTeTeTePi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP3,
            4,
            0,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_DiDiPi2Pi2
      (
        AddEnum
        (
          "N_DiDiPi2Pi2",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP,
            2,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>::Create( e_Py, e_Pz),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_Di2DiPi2Pi
      (
        AddEnum
        (
          "N_Di2DiPi2Pi",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP,
            1,
            1,
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            storage::Set< AtomicOrbitalTypesEnum>( e_Py),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_Tr2TrTrPi2
      (
        AddEnum
        (
          "N_Tr2TrTrPi2",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP2,
            2,
            1,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>( e_Pz),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      N_Te2Te2TeTe
      (
        AddEnum
        (
          "N_Te2Te2TeTe",
          AtomTypeData
          (
            GetElementTypes().e_Nitrogen,
            GetHybridOrbitalTypes().SP3,
            2,
            2,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      H_
      (
        AddEnum
        (
          "H_",
          AtomTypeData
          (
            GetElementTypes().e_Hydrogen,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            13.60 * 2.0,
            13.60,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Li_
      (
        AddEnum
        (
          "Li_",
          AtomTypeData
          (
            GetElementTypes().e_Lithium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            5.39 * 2.0,
            5.39,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      Na_
      (
        AddEnum
        (
          "Na_",
          AtomTypeData
          (
            GetElementTypes().Sodium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            5.14 * 2.0,
            5.14,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      ),
      K_
      (
        AddEnum
        (
          "K_",
          AtomTypeData
          (
            GetElementTypes().e_Potassium,
            GetHybridOrbitalTypes().Unhybridized,
            0,
            0,
            storage::Set< AtomicOrbitalTypesEnum>(),
            storage::Set< AtomicOrbitalTypesEnum>(),
            4.341 * 2.0,
            4.341,
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>(),
            util::GetUndefined< double>()
          )
        )
      )
 */
}

int
main( int argc, char * argv [] )
{

try {

	devel::init(argc, argv);

	utility::vector1< core::chemical::ElementCOP > elements( initialize_elements() );

	std::ofstream eout;
	eout.open( "element_data.txt" ); // Hardcoded - deal with it

	for ( core::Size ii(1); ii <= elements.size(); ++ii ) {
		eout << *(elements[ii]);
	}
	eout << std::endl;

	eout.close();

	// Now read it back in
	std::ifstream ein;
	ein.open( "element_data.txt" );
	core::Size ecount( 0 );
	while ( ein.good() ) {
		std::streampos start( ein.tellg() );
		std::string tag;
		ein >> tag;
		if ( ! ein.good() ) { // Expended all the values
			break;
		}
		ein.seekg( start );
		core::chemical::Element etd;
		ein >> etd;
		//std::cout << "read data for " << etd.GetChemicalSymbol() << " -- " << etd.GetChemicalName() << std::endl;
		//std::cerr << etd;
		++ecount;
	}
	std::cout << ">>> read " << ecount << " elements versus " << elements.size() << " total " << std::endl << std::endl;

	// Atoms now.

	core::chemical::ElementSetOP ele_set = new core::chemical::ElementSet;
	ele_set->read_file( "element_data.txt" );

	utility::vector1< core::chemical::gasteiger::GasteigerAtomTypeDataCOP > atoms( initialize_atoms( ele_set ) );

	std::ofstream aout;
	aout.open( "atom_type_data.txt" ); // Hardcoded - deal with it

	for ( core::Size ii(1); ii <= atoms.size(); ++ii ) {
		aout << *(atoms[ii]);
	}

	aout.close();

	// Now read it back in

	std::ifstream ain;
	ain.open( "atom_type_data.txt" );
	core::Size acount( 0 );
	while ( ain.good() ) {
		std::streampos start( ain.tellg() );
		std::string tag;
		ain >> tag;
		if ( ! ain.good() ) { // Expended all the values
			break;
		}
		ain.seekg( start );
		core::chemical::gasteiger::GasteigerAtomTypeData atd;
		atd.read( ain, ele_set );
		std::cout << "read data for " << atd.get_name() << std::endl;
		std::cerr << atd;
		++acount;
	}
	std::cout << ">>> read " << acount << " atom types versus " << atoms.size() << " total " << std::endl << std::endl;
	// read back in.

} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
	return -1;
}

}

