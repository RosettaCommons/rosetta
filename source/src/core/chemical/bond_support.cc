// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/chemical/ResidueSupport.hh
/// @brief support functions for class Bond; functions that
/// should not be included as part of the class.
/// @author Steven Combs
#include <core/chemical/bond_support.hh>
#include <core/chemical/ResidueSupport.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/Bond.hh>
#include <utility/graph/RingDetection.hh>

namespace core {
namespace chemical {
void find_bonds_in_rings(ResidueType & res){
	//first, we assign all the bonds in the residue to having no rings
	EIter edge_begin, edge_end;
	boost::tie(edge_begin, edge_end) = boost::edges( res.graph() );
	for(EIter edge_iter = edge_begin; edge_iter != edge_end; ++edge_iter){
		Bond & bond = res.bond(*edge_iter);
		//for now, nothing is in the ring. Defaulted to not known when constructed, after ring detection, either
		//in a ring or not in a ring
		bond.ringness( BondNotInRing );
	}

	//first we get the light weight residue graph
	LightWeightResidueGraph lwrg = convert_residuetype_to_light_graph(res);
	//now get the property maps
	boost::property_map<LightWeightResidueGraph, boost::vertex_name_t>::type lwrg_vd_to_VD = boost::get(boost::vertex_name, lwrg);
	//boost::property_map<LightWeightResidueGraph, boost::edge_name_t>::type lwrg_ed_to_ED = boost::get(boost::edge_name, lwrg);


	utility::graph::RingDetection ring_detect(lwrg); //initialize the ring detector. Automatically assigns rings
	utility::vector1<utility::vector1<lwrg_VD> > rings = ring_detect.GetRings(); //these are the path of the rings

	//iterate through the rings, then assign the bonds for ringness
	for(core::Size i=1; i<= rings.size(); ++i){
		utility::vector1< ED > just_the_edges;
		for(core::Size j=1; j<  rings[i].size(); ++j){ //not less than, see explanation below
			VD source =  lwrg_vd_to_VD[ rings[i][j] ];
			//next set of code is a little convulted. The ring code returns all the vertex, but we need
			//the edges to assign (bond), not the vertex. (atom, at least this point in time).
			//Therefore, we have to move one past the vector to get the edge.
			VD target = lwrg_vd_to_VD[ rings[i][ j+1] ];
			ED edge;
			bool edge_exists;
			boost::tie(edge, edge_exists) = boost::edge(source, target, res.graph());
			if(edge_exists){ //if there is an edge, mark it as being a ring
				just_the_edges.push_back(edge); //get the edge
				Bond & bond = res.bond(edge);
				bond.ringness(BondInRing);
			}
		}
	}
}

utility::vector1<VD> get_connecting_atoms(ResidueType const & res, ED const & edge) {
	return get_connecting_atoms(res.graph(), edge);
}

utility::vector1<VD> get_connecting_atoms(ResidueGraph const & graph, ED const & edge){
	utility::vector1<VD> connecting_atoms;
	connecting_atoms.push_back( boost::source( edge, graph ) );
	connecting_atoms.push_back( boost::target( edge, graph ) );
	return connecting_atoms;
}


ED get_bond(ResidueType const & res, VD const & source, VD const & target){
	bool bond_there(false);
	ED edge;
	boost::tie(edge,bond_there) = boost::edge(source, target, res.graph() );
	assert(bond_there);
	return edge;
}

}
}
