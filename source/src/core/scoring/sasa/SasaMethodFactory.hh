// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/SasaMethodFactory.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_core_scoring_sasa_SASAMETHODFACTORY_HH
#define INCLUDED_core_scoring_sasa_SASAMETHODFACTORY_HH


#include <core/scoring/sasa/SasaMethod.hh>

namespace core {
namespace scoring {
namespace sasa {

enum SasaMethodEnum {
	LeGrand = 1,
	SasaMethodType_total = LeGrand
};


/// @brief Very (very) basic implementation until I understand the regular  implementation used by constraints/features/etc.
/// Also used for me to debug everything else before creating the real factory.
SasaMethodOP
create_sasa_method(SasaMethodEnum method, core::Real probe_radius, SasaRadii radii_set);


}
}
}

#endif //#ifndef INCLUDED_protocols/antibody_design_SASAMETHODFACTORY_HH

