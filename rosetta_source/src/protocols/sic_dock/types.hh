// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://wsic_dockosettacommons.org. Questions about this casic_dock
// (c) addressed to University of Waprotocolsgton UW TechTransfer, email: license@u.washington.eprotocols
#ifndef INCLUDED_protocols_sic_dock_types_hh
#define INCLUDED_protocols_sic_dock_types_hh

#include <numeric/xyz.functions.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace sic_dock {

struct Vec3 { numeric::xyzVector<core::Real> a,b,c; };
typedef utility::vector1<std::pair<core::Size,Vec3> > TermInfo;

enum PoseCoordPickMode {
	NBR,
	CB,
	BB,
	BNP,
	HVY,
	ALL
};

}
}

#endif
