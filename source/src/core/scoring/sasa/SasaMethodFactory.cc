// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody_design/SasaMethodFactory.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/scoring/sasa/SasaMethodFactory.hh>
#include <core/scoring/sasa/LeGrandSasa.hh>

namespace core {
namespace scoring {
namespace sasa {
			

SasaMethodOP
create_sasa_method(SasaMethodEnum method, core::Real probe_radius, SasaRadii radii_set){
	
//	switch(method){
//		case LeGrand:
//			return new LeGrandSasa(probe_radius, radii_set);
//	}
	
	
return new LeGrandSasa(probe_radius, radii_set);
}
	
	
	
}
}
}
