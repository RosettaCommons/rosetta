// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/design/AntibodyDesignEnum.hh
/// @brief Antibody Design enumerators
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_ANTIBODYDESIGNENUM_hh
#define INCLUDED_protocols_antibody_design_ANTIBODYDESIGNENUM_hh

namespace protocols {
namespace antibody {
namespace design {

enum MinTypeEnum{
	
	relax = 1,
	centroid_relax,
	minimize,
	minimize_cartesian,
	dualspace,
	repack,
	backrub_protocol,
	no_min,
	MinTypeEnum_total = no_min
			
};


//Can be moved generally somewhere else
enum SeqDesignStrategyEnum {
	seq_design_profiles = 1,
	seq_design_conservative,
	seq_design_profile_sets,
	seq_design_profile_sets_combined,
	seq_design_basic,
	seq_design_none,
	
	SeqDesignStrategyEnum_total = seq_design_none
};

enum AntibodyDesignProtocolEnum {
	generalized_monte_carlo = 1,
	deterministic_graft,
	
	DesignProtocolEnum_total = deterministic_graft
};


}
}
}


#endif	//#ifndef INCLUDED_protocols/antibody_design_ANTIBODYDESIGNENUM_HH

