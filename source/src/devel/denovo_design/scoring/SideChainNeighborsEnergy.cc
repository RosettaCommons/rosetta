// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SideChainNeighborsEnergy.cc
/// @author Tom Linsky

#include <devel/denovo_design/scoring/SideChainNeighborsEnergy.hh>
#include <devel/denovo_design/scoring/SideChainNeighborsEnergyCreator.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>

#include <core/pose/Pose.hh>

#include <basic/prof.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <ObjexxFCL/string.functions.hh>

#include <basic/Tracer.hh>
#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR("core.scoring.methods.SideChainNeighborsEnergy");

namespace core {
namespace scoring {
namespace methods {

/// @details This must return a fresh instance of the SideChainNeighborsEnergy class,
/// never an instance already in use
methods::EnergyMethodOP
SideChainNeighborsEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const & ) const
{
	return methods::EnergyMethodOP( new SideChainNeighborsEnergy() );
}

ScoreTypes
SideChainNeighborsEnergyCreator::score_types_for_method() const
{
	ScoreTypes sts;
	sts.push_back(sidechain_neighbors);
	return sts;
}

void
SideChainNeighborsEnergy::setup_for_scoring(
	pose::Pose & pose,
	ScoreFunction const & ) const
{
	pose.update_residue_neighbors();
}

/// clone
EnergyMethodOP
SideChainNeighborsEnergy::clone() const
{
	return EnergyMethodOP( new SideChainNeighborsEnergy() );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

SideChainNeighborsEnergy::SideChainNeighborsEnergy() :
	ContextDependentOneBodyEnergy( methods::EnergyMethodCreatorOP( new SideChainNeighborsEnergyCreator() ) )
{
}

void
SideChainNeighborsEnergy::residue_energy(
	core::conformation::Residue const & rsd,
	core::pose::Pose const & pose,
	EnergyMap & emap ) const
{
	TwelveANeighborGraph const & graph ( pose.energies().twelveA_neighbor_graph() );

	Real my_neighbors = 0.0;
	Vector my_sc_coordinates, my_bb_coordinates;
	if ( rsd.name3() == "GLY" ) {
		if ( rsd.type().residue_type_set()->name() == "fa_standard" ) {
			my_sc_coordinates = rsd.atom(rsd.atom_index("2HA")).xyz();
		} else if ( rsd.type().residue_type_set()->name() == "centroid" ) {
			my_sc_coordinates = rsd.atom(rsd.atom_index("CEN")).xyz();
		} else {
			throw utility::excn::EXCN_BadInput( "Unknown residue type set for gly residue: " + rsd.type().residue_type_set()->name() );
		}
		my_bb_coordinates = rsd.atom(rsd.atom_index("CA")).xyz() ;
	} else {
		core::Size const sc_atom_index = rsd.first_sidechain_atom();
		core::Size const parent_atom_index = rsd.icoor( sc_atom_index ).stub_atom1().atomno();
		// if main chain/sidechain atoms are the same, just ignore this residue.
		// the most likely cause is that this is a ligand residue
		if ( sc_atom_index == parent_atom_index ) {
			return;
		}
		my_sc_coordinates = rsd.atom( sc_atom_index ).xyz();
		my_bb_coordinates = rsd.atom( parent_atom_index ).xyz();
	}

	// compute normalized CA -> CB vector
	Vector const my_sc_vector = (my_sc_coordinates - my_bb_coordinates).normalize();

	// iterate across neighbors within 12 angstroms, add sum of vectors
	for ( graph::Graph::EdgeListConstIter
			ir  = graph.get_node( rsd.seqpos() )->const_edge_list_begin();
			ir != graph.get_node( rsd.seqpos() )->const_edge_list_end(); ++ir ) {
		Size const j = (*ir)->get_other_ind( rsd.seqpos() );
		core::conformation::Residue const & rsd_j = pose.residue(j);

		Vector other_bb_coordinates;
		if ( rsd_j.name3() == "GLY" ) {
			//my_sc_coordinates = rsd_j.atom(pose.residue(j).atom_index("2HA")).xyz();
			other_bb_coordinates = rsd_j.atom(pose.residue(j).atom_index("CA")).xyz();
		} else {
			//my_sc_coordinates = rsd_j.atom(pose.residue(j).first_sidechain_atom()).xyz() ;
			core::Size const parent_atom_index = rsd_j.icoor( rsd_j.first_sidechain_atom() ).stub_atom1().atomno();
			other_bb_coordinates = rsd_j.atom( parent_atom_index ).xyz();
		}

		Vector new_sc_vector = other_bb_coordinates - my_sc_coordinates;
		Real const distance = 1.0 / ( 1.0 + exp(new_sc_vector.length() - 9.0 ) );
		Real angle = my_sc_vector.dot( new_sc_vector.normalize() ) + 0.5;
		if ( angle < 0 ) {
			angle = 0.0;
		}
		my_neighbors += ( distance * angle * angle );
	}

	Real const score = -my_neighbors / 2.25;
	emap[ sidechain_neighbors ] += score;
	TR.Debug << rsd.seqpos() << " " << score << std::endl;
}

core::Size
SideChainNeighborsEnergy::version() const
{
	return 1;
}

void
SideChainNeighborsEnergy::indicate_required_context_graphs( utility::vector1< bool > & context_graphs_required ) const
{
	context_graphs_required[ twelve_A_neighbor_graph ] = true;
}

} // methods
} // scoring
} // core
