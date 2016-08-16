// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
// This file is based off of code from the Biochemistry Library (BCL).
// The BCL is copyright Vanderbilt University (Meiler Lab), a RosettaCommons member

/// @file   core/chemical/Elements.hh
/// @brief  All the elements known to mankind (me)...as of4/18/2014
/// @author Steven Combs

#ifndef INCLUDED_core_chemical_Elements_hh
#define INCLUDED_core_chemical_Elements_hh
#include <ostream>
#include <map>
#include <string>
#include <utility/vector1.hh>
#include <utility/vector0.hh>

namespace core {
namespace chemical {
namespace element {

enum Elements {
	UnknownElement=0,
	H=1,
	He,
	Li,
	Be,
	B,
	C,
	N,
	O,
	F,
	Ne,
	Na,
	Mg,
	Al,
	Si,
	P,
	S,
	Cl,
	Ar,
	K,
	Ca,
	Sc,
	Ti,
	V,
	Cr,
	Mn,
	Fe,
	Co,
	Ni,
	Cu,
	Zn,
	Ga,
	Ge,
	As,
	Se,
	Br,
	Kr,
	Rb,
	Sr,
	Y,
	Zr,
	Nb,
	Mo,
	Tc,
	Ru,
	Rh,
	Pd,
	Ag,
	Cd,
	In,
	Sn,
	Sb,
	Te,
	I,
	Xe,
	Cs,
	Ba,
	La,
	Ce,
	Pr,
	Nd,
	Pm,
	Sm,
	Eu,
	Gd,
	Tb,
	Dy,
	Ho,
	Er,
	Tm,
	Yb,
	Lu,
	Hf,
	Ta,
	W,
	Re,
	Os,
	Ir,
	Pt,
	Au,
	Hg,
	Tl,
	Pb,
	Bi,
	Po,
	At,
	Rn,
	Fr,
	Ra,
	Ac,
	Th,
	Pa,
	U,
	Np,
	Pu,
	Am,
	Cm,
	Bk,
	Cf,
	Es,
	Fm,
	Md,
	No,
	Lr,
	Rf,
	Db,
	Sg,
	Bh,
	Hs,
	Mt,
	Ds,
	Rg,
	Cn,
	Uut,
	Fl,
	Uup,
	Lv,
	Uus,
	Uuo,
	total_number_elements = 119
};

utility::vector0< std::string > & element2name();
std::string
name_from_elements(Elements element);
Elements
elements_from_name(std::string name);


std::map< std::string, Elements > & name2element();

/// @brief setup the map that converts string name to AA enum
std::map< std::string, Elements > setup_name2element();


/// @brief setup the vector that maps AA enum to string name
utility::vector0< std::string > setup_element2name();

}
}
}


#endif /* ELEMENTS_HH_ */
