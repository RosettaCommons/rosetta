// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/interaction_graph/NPDHBondInteractionGraph.cc
/// @brief  Non-templated helper functions for the NPDHBondInteractionGraph
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#include <core/pack/interaction_graph/NPDHBondInteractionGraph.hh>

//Rosetta Headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueType.hh>

#include <core/conformation/Residue.hh>

#include <utility/graph/DisjointSets.hh>

#include <core/pack/interaction_graph/AdditionalBackgroundNodesInteractionGraph.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/interaction_graph/LinearMemoryInteractionGraph.hh>
#include <core/pack/interaction_graph/RotamerDots.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>

#include <core/pose/Pose.hh>

//#include <core/scoring/sasa.hh>
#include <core/scoring/sasa/util.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/pack/interaction_graph/SurfacePotential.hh>

#include <basic/Tracer.hh>

//Utility Headers
#include <utility/vector1.hh>
#include <utility/exit.hh>

namespace core {
namespace pack {
namespace interaction_graph {


void
compute_alt_weights_for_npd_hbonds(
	conformation::Residue const & res,
	utility::vector1< utility::vector1< NPDHBondOP > > const & atom_hbonds,
	utility::vector1< Real > & tmp_energies,
	utility::vector1< Real > & tmp_weights
)
{
	for ( Size ii = 1; ii <= res.natoms(); ++ii ) {
		if ( atom_hbonds[ ii ].size() == 0 ) continue;
		tmp_energies.resize( 0 ); tmp_energies.reserve( atom_hbonds[ ii ].size() );
		tmp_weights.resize( atom_hbonds[ ii ].size() );
		std::transform( atom_hbonds[ ii ].begin(), atom_hbonds[ ii ].end(),
			std::back_inserter( tmp_energies ), []( NPDHBondOP const & hb ) { return hb->energy_ * hb->sfxn_wt_; } );

		// invoke the function from NPDHBondSet to compute the weights...
		scoring::hbonds::weights_for_hbonds( res, ii, tmp_energies, tmp_weights );

		// ... and then distribute the weights returned to the alt_ values in the hbond objects.
		// This allows us to reuse a hydrogen bond whose energy has not changed, but whose weight might have.
		if ( res.atom_is_polar_hydrogen( ii ) ) {
			for ( Size jj = 1; jj <= atom_hbonds[ ii ].size(); ++jj ) {
				atom_hbonds[ ii ][ jj ]->don_wt_alt_ = tmp_weights[ jj ];
			}
		} else {
			debug_assert( res.heavyatom_is_an_acceptor( ii ) );
			for ( Size jj = 1; jj <= atom_hbonds[ ii ].size(); ++jj ) {
				atom_hbonds[ ii ][ jj ]->acc_wt_alt_ = tmp_weights[ jj ];
			}
		}
	}

}


} //end namespace
} //end namespace
} //end namespace

