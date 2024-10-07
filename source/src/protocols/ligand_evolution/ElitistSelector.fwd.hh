// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_evolution/ElitistSelector.fwd.hh
/// @brief  Forward declaration of the %ElitistSelector class
/// @author Paul Eisenhuth (eisenhuth451@gmail.com)

#ifndef INCLUDED_protocols_ligand_evolution_ElitistSelector_FWD_HH
#define INCLUDED_protocols_ligand_evolution_ElitistSelector_FWD_HH

// utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols{
namespace ligand_evolution{

    class ElitistSelector;

    typedef utility::pointer::shared_ptr< ElitistSelector > ElitistSelectorOP;
    typedef utility::pointer::shared_ptr< ElitistSelector const > ElitistSelectorCOP;

}
}

#endif
