// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/SasaMethodFactory.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/scoring/sasa/SasaMethodFactory.hh>
#include <core/scoring/sasa/LeGrandSasa.hh>

namespace core {
namespace scoring {
namespace sasa {


SasaMethodOP
create_sasa_method(
	SasaMethodEnum /*method*/,
	core::Real probe_radius,
	SasaRadii radii_set
) {
	return SasaMethodOP( new LeGrandSasa(probe_radius, radii_set) );
}



}
}
}
