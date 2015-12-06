// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// Rosetta Headers
#include <protocols/farna/fragments/RNA_Fragments.hh>
#include <protocols/toolbox/AtomLevelDomainMap.fwd.hh>


#include <core/types.hh>

#include <utility/vector1.hh>
#include <iostream>


namespace protocols {
namespace farna {
namespace fragments {

RNA_Fragments::RNA_Fragments(){}
RNA_Fragments::~RNA_Fragments(){}

// empty shell.
void
RNA_Fragments::apply_random_fragment(
	core::pose::Pose & /*pose*/,
	core::Size const /*position*/,
	core::Size const /*size*/,
	core::Size const /*type = MATCH_YR*/,
	toolbox::AtomLevelDomainMapCOP  /*atom_level_domain_map*/ ) const {
	std::cout << "Should not be in here! " << std::endl;
}

bool
RNA_Fragments::is_fullatom(){ return false; }

} //fragments
} //farna
} //protocols

