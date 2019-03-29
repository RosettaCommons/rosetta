// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/ga_dock/util.cc
///
/// @brief
/// @author Hahnbeom Park and Frank DiMaio

#include <protocols/ligand_docking/GALigandDock/util.hh>

#include <core/chemical/AA.hh>
#include <core/scoring/Energies.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
//#include <utility/graph/Graph.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

using namespace core;

static basic::Tracer TR( "protocols.ligand_docking.GALigandDock" );

core::Size
count_neighbors_on_coord( core::pose::Pose const &pose,
	core::Vector const &xyz1,
	std::string const atomname,
	core::Real const dcut )
{
	core::Size neighbor_counts( 0 );
	core::Real const dcut2( dcut*dcut );

	for ( core::Size resno = 1; resno <= pose.size(); ++resno ) {
		if ( pose.residue(resno).is_virtual_residue() ) continue;

		core::Vector xyz2;
		if ( atomname == "nbr" ) {
			xyz2 = pose.residue(resno).nbr_atom_xyz();
		} else {
			xyz2 = pose.residue(resno).has(atomname)? pose.residue(resno).xyz(atomname)
				: pose.residue(resno).nbr_atom_xyz();
		}
		core::Real d2 = xyz1.distance_squared(xyz2);
		if ( d2 < dcut2 ) neighbor_counts++;
	}

	return neighbor_counts;
}

utility::vector1< core::Size >
count_neighbors( core::pose::Pose const &pose,
	std::string const atomname,
	core::Real const dcut )
{
	utility::vector1< core::Size > neighbor_counts( pose.size(), 0 );

	core::scoring::TwelveANeighborGraph const & graph = pose.energies().twelveA_neighbor_graph();
	core::Real const dcut2( dcut*dcut );

	for ( core::Size resno = 1; resno <= pose.size(); ++resno ) {
		if ( pose.residue(resno).is_virtual_residue() ) continue;
		core::Vector xyz1;
		if ( atomname == "nbr" ) {
			xyz1 = pose.residue(resno).nbr_atom_xyz();
		} else {
			xyz1 = pose.residue(resno).has(atomname)? pose.residue(resno).xyz(atomname) : pose.residue(resno).nbr_atom_xyz();
		}

		for ( utility::graph::Graph::EdgeListConstIter
				ir  = graph.get_node( resno )->const_edge_list_begin(),
				ire = graph.get_node( resno )->const_edge_list_end();
				ir != ire; ++ir ) {
			Size const j( (*ir)->get_other_ind( resno ) );

			// supporting residue level only yet...
			core::Vector xyz2;
			if ( atomname == "nbr" ) {
				xyz2 = pose.residue(j).nbr_atom_xyz();
			} else {
				xyz2 = pose.residue(j).has(atomname)? pose.residue(j).xyz(atomname) : pose.residue(j).nbr_atom_xyz();
			}
			core::Real d2 = xyz1.distance_squared(xyz2);
			if ( d2 < dcut2 ) neighbor_counts[resno]++;
		}
	}

	return neighbor_counts;
}

utility::vector1< core::Size >
get_atomic_contacting_sidechains( core::pose::Pose const &pose,
	core::Size const ligid,
	core::Real const atomic_distance_cut )
{
	utility::vector1< core::Size > contact_scs;

	core::Real const D2_COARSE( 225.0 );
	core::Real const D2_FINE( atomic_distance_cut*atomic_distance_cut );

	core::Vector const &ligcom = pose.residue(ligid).nbr_atom_xyz();
	utility::vector1< bool > is_close( pose.size(), false );
	for ( core::Size ires = 1; ires <= pose.size(); ++ires ) {
		if ( !pose.residue(ires).is_protein() ) continue;
		if ( pose.residue(ires).aa() == core::chemical::aa_gly ||
				pose.residue(ires).aa() == core::chemical::aa_ala ||
				pose.residue(ires).aa() == core::chemical::aa_pro ) continue;

		core::Vector const & resnbr = pose.residue(ires).nbr_atom_xyz();
		core::Real d2 = ligcom.distance_squared( resnbr );
		if ( d2 < D2_COARSE ) is_close[ires] = true;
	}

	for ( core::Size ires = 1; ires <= pose.size(); ++ires ) {
		if ( !is_close[ires] || ires == ligid ) continue;

		bool i_is_contacting = false;
		for ( core::Size iatm = 1; iatm <= pose.residue(ires).nheavyatoms(); ++iatm ) {
			if ( pose.residue(ires).atom_is_backbone(iatm) ) continue;
			core::Vector const &xyz_i = pose.residue(ires).xyz(iatm);

			for ( core::Size jatm = 1; jatm <= pose.residue(ligid).nheavyatoms(); ++jatm ) {
				core::Vector const &xyz_j = pose.residue(ligid).xyz(jatm);

				core::Real d2 = xyz_i.distance_squared( xyz_j );
				if ( d2 < D2_FINE ) {
					i_is_contacting = true;
					break;
				}
			}
			if ( i_is_contacting ) break;
		} // iatm

		if ( i_is_contacting ) {
			contact_scs.push_back( ires );
			/*
			for( core::Size ichi = 1; ichi <= pose.residue(ires).nchi(); ++ichi ){
			if( !pose.residue(ires).type().is_proton_chi(ichi)  ){
			chidefs.push_back( std::make_pair( ires, ichi ) );
			std::cout << "chidef: " << chidefs.size() << " " << ires << " " << ichi <<std::endl;
			}
			}
			*/
		}
	} // ires

	return contact_scs;
}

} // ga_dock
} // ligand_docking
} // protocols


