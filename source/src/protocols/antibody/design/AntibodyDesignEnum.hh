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

#ifndef INCLUDED_protocols_antibody_design_ANTIBODYDESIGNENUM_HH
#define INCLUDED_protocols_antibody_design_ANTIBODYDESIGNENUM_HH

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
	no_min,
	MinTypeEnum_total = no_min
			
};

enum DesignTypeEnum{
	
	relaxed_design = 1,
	docked_design,
	fixbb,
	flxbb,
	DesignTypeEnum_total = flxbb
			
};

//Can be moved generally somewhere else
enum SeqDesignStrategyEnum {
	seq_design_profiles = 1,
	seq_design_conservative,
	seq_design_basic,
	
	SeqDesignStrategyEnum_total = seq_design_basic
};


}
}
}


#endif	//#ifndef INCLUDED_protocols/antibody_design_ANTIBODYDESIGNENUM_HH

