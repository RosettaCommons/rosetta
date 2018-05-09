// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/core/conformation/AtomGraph.cc
/// @author Sam DeLuca

/// @detail right now there isnt a canonical mapping from nodes back to atom objects.

#include <core/conformation/AtomGraph.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>

#include <core/conformation/AtomGraphData.hh>
#include <utility/graph/UpperEdgeGraph.hh>


namespace core {
namespace conformation {

void
atom_graph_from_conformation(
	Conformation const & conformation,
	AtomGraphOP atom_graph
)
{
	//TODO: doing this because I don't want to pull in an entire pose :( gotta be a better way
	platform::Size num_atoms = 0;
	for ( platform::Size resid = 1; resid <= conformation.size(); ++resid ) {
		num_atoms += conformation.residue_type(resid).natoms();
	}
	num_atoms++;
	atom_graph->set_num_vertices(num_atoms);
	platform::Size index_id = 1;

	for ( platform::Size resid=1 ; resid <= conformation.size(); ++resid ) {
		core::conformation::Residue current_res(conformation.residue(resid));
		for ( platform::Size atomno=1; atomno <= current_res.natoms(); ++atomno ) {
			AtomGraphVertexData current_vertex = atom_graph->get_vertex(index_id).data();
			current_vertex.xyz() = current_res.xyz(atomno);
			current_vertex.atom_name() = current_res.atom_name(atomno);
			core::Real lj_radius = current_res.atom_type(atomno).lj_radius();
			current_vertex.atom_radius_squared() = lj_radius*lj_radius;


			num_atoms++;
		}
	}

}


platform::Size
annotated_atom_graph_from_conformation(
	Conformation const & conformation,
	AtomGraphOP atom_graph,
	PointPosition const & additional_point

)
{
	//TODO: doing this because I don't want to pull in an entire pose :( gotta be a better way
	platform::Size num_atoms = 0;
	for ( platform::Size resid = 1; resid <= conformation.size(); ++resid ) {
		num_atoms += conformation.residue_type(resid).natoms();
	}
	num_atoms++;
	atom_graph->set_num_vertices(num_atoms);
	platform::Size index_id = 1;
	for ( platform::Size resid=1 ; resid <= conformation.size(); ++resid ) {
		core::conformation::Residue current_res(conformation.residue(resid));
		for ( platform::Size atomno=1; atomno <= current_res.natoms(); ++atomno ) {
			AtomGraphVertexData current_vertex = atom_graph->get_vertex(index_id).data();
			current_vertex.xyz() = current_res.xyz(atomno);
			current_vertex.atom_name() = current_res.atom_name(atomno);
			current_vertex.residue_id() = resid;
			//num_atoms++;
		}
	}
	//add the additional point
	atom_graph->get_vertex(num_atoms).data().xyz() = additional_point;

	return num_atoms; //this is the ID of the additional point
}


}
}
