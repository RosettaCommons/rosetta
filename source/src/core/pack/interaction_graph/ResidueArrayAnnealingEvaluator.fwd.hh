// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/ResidueArrayAnnealingEvaluator.hh
/// @brief  Annealable method interrface for score types evaluated over explicit list of residues.
/// @author Alex Ford (fordas@uw.edu)

#ifndef INCLUDED_core_pack_interaction_graph_ResidueArrayAnnealingEvaluator_fwd_hh
#define INCLUDED_core_pack_interaction_graph_ResidueArrayAnnealingEvaluator_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


namespace core {
namespace pack {
namespace interaction_graph {

class ResidueArrayAnnealingEvaluator;

typedef utility::pointer::shared_ptr< ResidueArrayAnnealingEvaluator > ResidueArrayAnnealingEvaluatorOP;

}
}
}

#endif
