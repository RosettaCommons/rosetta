// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/chemical/residue_support.hh
/// @brief support functions for class residue; functions that
/// should not be included as part of the class.
/// @author Phil Bradley
/// @author Rocco Moretti (rmorettiase@gmail.com)
/// @author Steven Combs

// Package Headers
#include <core/chemical/residue_support.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/ChemicalManager.hh>

// Project Headers
#include <utility/graph/Graph.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/Bond.hh>
#include <utility/graph/RingDetection.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/graph/BFS_prune.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>

#include <algorithm>

namespace core {
namespace chemical {

static basic::Tracer TR( "core.chemical.residue_support" );

ObjexxFCL::FArray2D_int
get_residue_path_distances( ResidueType const & res )
{
	//return res.get_residue_path_distances();
	using namespace utility::graph;
	Graph g;

	g.set_num_nodes( res.natoms() );
	for ( uint ii = 1; ii <= res.natoms(); ++ii ) {
		AtomIndices const ii_bonded = res.nbrs( ii  );
		for ( Size jj = 1; jj <= ii_bonded.size(); ++jj ) {
			if ( ii_bonded[ jj ] > ii ) {
				g.add_edge( ii, ii_bonded[ jj ] );
			}
		}
	}
	return g.all_pairs_shortest_paths();
}

LightWeightResidueGraph convert_residuetype_to_light_graph( ResidueType const & res ){
	//this is a const reference because vertices change when  you copy the graph
	const  core::chemical::ResidueGraph & full_residue_graph = res.graph(); //get the boost graph structure from residuetype

	LightWeightResidueGraph lwrg;

	//set up the mapping between VD and ED from ResidueGraph for LightWeightResidueGraph
	boost::property_map<LightWeightResidueGraph, boost::vertex_name_t>::type lwrg_vd_to_VD = boost::get(boost::vertex_name, lwrg);
	boost::property_map<LightWeightResidueGraph, boost::edge_name_t>::type lwrg_ed_to_ED = boost::get(boost::edge_name, lwrg);

	//map between the LightWeightResidueGraph vertex and the ResidueGraph vertex. Used when adding edges
	std::map<core::chemical::VD, lwrg_VD> map_of_vertex;
	//first, add all the vertex to light weight residue graph, and map the property maps to point at the ResidueGraph
	for ( core::chemical::VIterPair vp = boost::vertices(full_residue_graph); vp.first != vp.second; ++vp.first ) {
		VD const & full_residue_graph_vd = *vp.first;
		lwrg_VD lwrg_vd = boost::add_vertex(lwrg);
		lwrg_vd_to_VD[lwrg_vd] = full_residue_graph_vd; //set property maps
		map_of_vertex[full_residue_graph_vd] = lwrg_vd; //set mapping for the edges
	}

	//now we add the edges between the vertex
	for ( EIterPair ep = boost::edges(full_residue_graph); ep.first != ep.second; ++ep.first ) {
		VD source = boost::source(*ep.first, full_residue_graph);
		VD target = boost::target(*ep.first, full_residue_graph);
		lwrg_VD source_lwrg_VD = map_of_vertex[source];
		lwrg_VD target_lwrg_VD = map_of_vertex[target];
		lwrg_ED e_added;
		bool added;
		boost::tie(e_added, added) = boost::add_edge(source_lwrg_VD, target_lwrg_VD,lwrg );
		if ( added ) { //only add bonds once!
			lwrg_ed_to_ED[e_added] = *ep.first;
		}
	}
	debug_assert(boost::num_vertices(lwrg) == full_residue_graph.num_vertices()); //fail if the number of vertex are not the same
	debug_assert(boost::num_edges(lwrg) == full_residue_graph.num_edges()); //fail if the number of edges are not the same

	//boost::property_map<LightWeightResidueGraph, boost::vertex_name_t>::type lwrg_vd_to_VD = boost::get(boost::vertex_name, lwrg);
	//boost::property_map<LightWeightResidueGraph, boost::edge_name_t>::type lwrg_ed_to_ED = boost::get(boost::edge_name, lwrg);

	utility::graph::RingDetection<LightWeightResidueGraph> ring_detect(lwrg); //initialize the ring detector. Automatically assigns rings
	utility::vector1<utility::vector1<lwrg_VD> > rings = ring_detect.GetRings(); //these are the path of the rings

	//our light weight graph has been made. Now return it!
	return lwrg;
}

/// @brief Rename atoms in the residue type such that their names are unique.
/// If preserve is true, only rename those which have no names or who have
/// name conflicts. (Both of the conflicting atoms will be renamed.)
void
rename_atoms( ResidueType & res, bool preserve/*=true*/ ) {
	std::map< std::string, core::Size > name_counts;
	ResidueGraph const & graph( res.graph() );
	VIter iter, iter_end;
	for ( boost::tie( iter, iter_end ) = boost::vertices( graph ); iter != iter_end; ++iter ) {
		Atom const & atom( graph[*iter] );
		if ( preserve && atom.name().size() != 0 ) {
			name_counts[ atom.name() ]++; // with map, builtins are zero-initialized according to the C++ spec.
		}
	}

	for ( boost::tie( iter, iter_end ) = boost::vertices( graph ); iter != iter_end; ++iter ) {
		Atom const & atom( graph[*iter] );
		if ( preserve && name_counts[ atom.name() ] == 1 ) continue;

		//Find the first unoccupied name Xnnn type string.
		// Skipping values which were multiply represented in the input is deliberate
		// There's no fair way to choose which one is the "real" one.
		debug_assert( atom.element_type() );
		std::string name;
		core::Size ii(0);
		do {
			++ii;
			name = ObjexxFCL::uppercased( atom.element_type()->get_chemical_symbol() ) + utility::to_string( ii );
			//Align name, preferring to keep start in the second position
			if ( name.size() == 2 ) {
				name = ' ' + name + ' ';
			} else if ( name.size() == 3 ) {
				name = ' ' + name;
			} //Can't be 1, and if it's 4 or greater, leave as is.
		} while ( name_counts.find( name ) != name_counts.end() );
		// Assign new value and mark it used.
		if ( TR.Trace.visible() ) {
			TR.Trace << "Renaming atom from '"<< res.atom(*iter).name() << "' to '" << name << "'" << std::endl;
		}
		res.atom(*iter).name( name );
		res.remap_pdb_atom_names( true ); // We've renamed an atom, so we need to be flexible with PDB loading.
		name_counts[ name ]++;
	}
}

/// @brief Utility visitor for find_nbr_dist
/// Will only traverse those atoms in the "rigid" portion of graph around the starting atom.
/// "Rigid" includes direct neighbors and atoms connected by non-rotatable bonds
/// e.g. all rings, all double/triple bonds, methyl groups, single atoms, etc.
class RigidDistanceVisitor: public utility::graph::null_bfs_prune_visitor {
	using Matrix = utility::vector1<utility::vector1<core::Real> >;

public:
	RigidDistanceVisitor( Matrix & distances, ResidueType const & restype, VD start ) :
		distances_(distances),
		restype_(restype),
		start_(start),
		start_pos_( restype.atom(start).ideal_xyz() )
	{
		start_index_ = restype.atom_index( start );
	}

	template <class ED, class ResidueGraphType>
	bool examine_edge( ED edge, ResidueGraphType & graph) {
		VD source( boost::source( edge, graph ) ), target( boost::target( edge, graph ) );
		//std::string start( restype_.atom_name( start_ ) );
		//std::string source( restype_.atom_name( boost::source( edge, graph ) ) );
		//std::string target( restype_.atom_name( boost::target( edge, graph ) ) );
		if ( source == start_ ) {
			// Directly connected atoms are part of the rigid units, for distance purposes,
			// but we shouldn't propagate across them unless the fit the other criteria
			// (see TODO below, though)
			core::Size index( restype_.atom_index( target ) );
			distances_[start_index_][ index ] = restype_.atom( target ).ideal_xyz().distance( start_pos_ );
		}
		// Follow across non-rotatable bonds. (Ring and double bonds)
		Bond const & bond( restype_.bond( edge ) );
		if ( bond.ringness() == BondInRing ) {
			return false; // Follow rings
		}
		if ( bond.order() == DoubleBondOrder || bond.order() == TripleBondOrder ) {
			return false; // Follow double and triple bonds
		}
		if ( is_nub( source ) || is_nub( target ) ) {
			//Need to test both for symmetry - if either is a nub the bond isn't rotatable and should be followed.
			return false;
		}
		return true;
		// TODO: The atom *directly* across from the rotatable bond is ridgidly connected.
		// Do we want to include that in the unit? (molfile_to_params doesn't seem to).
	}

	template <class VD, class ResidueGraphType>
	bool examine_vertex( VD vertex, ResidueGraphType & ) {
		// If we're examining the vertex, it's in the rigid unit -- add distances to the matrix
		core::Size index( restype_.atom_index( vertex ) );
		distances_[start_index_][index] = restype_.atom( vertex ).ideal_xyz().distance( start_pos_ );
		// The reverse will be taken care of in a later iteration.
		return false; // Always examine out edges.
	}

	bool is_nub( VD atom ) {
		// Additional non-rotatable bonds include single bonds which are connected to atoms
		core::chemical::AdjacentIter iter, iter_end;
		core::Size nheavy(0), nhydro(0);
		for ( boost::tie( iter, iter_end ) = restype_.bonded_neighbor_iterators(atom); iter != iter_end; ++iter ) {
			if ( restype_.atom(*iter).element_type()->element() == core::chemical::element::H ) {
				++nhydro;
			} else {
				++nheavy;
			}
		}
		if ( nheavy >= 2 ) return false; // Multiply-heavy bonded atoms aren't non-rotatable stubs.
		// Non-considered bonds are all hydrogens (bonds to hydrogens use this function)
		if ( nhydro != 1 ) return true; // Multiple hydrogens aren't considered rotameric.
		if ( restype_.atom(atom).element_type()->element() == core::chemical::element::C ) return true; // Carbon with only a single heavy bond is never rotameric
		// Proton rotamer
		return false;
	}


private:
	Matrix & distances_;
	ResidueType const & restype_;
	VD start_;
	core::Size start_index_;
	Vector const & start_pos_;
};

/// @brief Calculate the rigid matrix - assume that distances has been initialized to some really large value, and is square
void calculate_rigid_matrix( ResidueType const & res, utility::vector1< utility::vector1< core::Real > > & distances ) {
	// Set up bonded and rigid distances
	VIter iter, iter_end;
	for ( boost::tie( iter, iter_end ) = res.atom_iterators(); iter != iter_end; ++iter ) {
		// The visitor takes the distance matrix as a reference and fills it out.
		RigidDistanceVisitor vis( distances, res, *iter );
		utility::graph::breadth_first_search_prune( res.graph(), *iter, vis );
	}
	// The Floyd-Warshall algorithm. We do this in-line instead of with boost
	// because some of the distance weights don't correspond to edge weights.
	// (Besides, it's easy enough.)
	for ( core::Size kk(1); kk <= res.natoms(); ++kk ) {
		for ( core::Size jj(1); jj <= res.natoms(); ++jj ) {
			for ( core::Size ii(1); ii <= res.natoms(); ++ii ) {
				core::Real new_dist( distances[ii][kk] + distances[kk][jj] );
				if ( new_dist < distances[ii][jj] ) {
					distances[ii][jj] = new_dist;
				}
			}
		}
	}

	if ( TR.Trace.visible() ) {
		// Print out the full distance matrix for debugging purposes.
		TR.Trace << std::setprecision(4);
		TR.Trace << "Atom distance matrix for " << res.name() << std::endl;
		TR.Trace << "    " << '\t';
		for ( core::Size xx(1); xx <= res.natoms(); ++xx ) {
			TR.Trace << res.atom_name(xx) << '\t';
		}
		TR.Trace << std::endl;
		for (  core::Size yy(1); yy <= res.natoms(); ++yy ) {
			TR.Trace << res.atom_name(yy) << '\t';
			for (  core::Size xx(1); xx <= res.natoms(); ++xx ) {
				TR.Trace << distances[yy][xx] << '\t';
			}
			TR.Trace << std::endl;
		}
	}
}

/// @brief Find the neighbor distance to the given neighbor atom.
/// If nbr_atom is null_vertex, give the smallest neighbor distance,
/// and set nbr_atom to the atom for that distance.
/// @details The neighbor distance here is adjusted for rotatable bonds -
/// It should be at least as large as the maximum neighbor distance
/// in any torsional rotamer
/// If the neighbor atom is not provided, the atom chosen will be a
/// multiply-bonded heavy atom.
///
/// Assumes:
///   * All atoms and bond are present
///   * All ideal_xyz coordinates have been set
///   * All elements have been set
///  * All ring bonds have been annotated
core::Real
find_nbr_dist( ResidueType const & res, VD & nbr_atom ) {
	if ( res.natoms() == 0 ) {
		utility_exit_with_message("Cannot find neighbor atom distance for empty residue type.");
	}
	core::Real maxdist = 1e9; // Hopefully sufficiently large.
	utility::vector1< utility::vector1< core::Real > > distances(res.natoms(), utility::vector1< core::Real >( res.natoms(), maxdist ) );
	calculate_rigid_matrix( res, distances );

	// TODO: Although we throw out hydrogens as potential neighbor atoms,
	// I believe the meaning of neighbor atoms is heavy-atom distances
	// We could get smaller values by removing hydrogens from distance considerations.

	utility::vector1< core::Real > maxdists;
	for ( core::Size ii(1); ii<= res.natoms(); ++ii ) {
		maxdists.push_back( *(std::max_element( distances[ii].begin(), distances[ii].end() )) );
	}

	if ( nbr_atom != ResidueType::null_vertex ) {
		// We have an atom in mind - we just need the distance.
		return maxdists[ res.atom_index( nbr_atom ) ];
	}

	//maxdist initialized above
	for ( core::Size jj(1); jj <= res.natoms(); ++jj ) {
		VD atom_vd( res.atom_vertex( jj ) );
		if ( maxdists[jj] < maxdist &&
				res.atom(jj).element_type()->element() != core::chemical::element::H ) { // not hydrogen
			//Use graph directly, as other neighbor annotations may not be fully baked
			core::Size bonded_heavy(0);
			AdjacentIter itr, itr_end;
			for ( boost::tie(itr, itr_end) = res.bonded_neighbor_iterators(atom_vd); itr != itr_end; ++itr ) {
				if ( res.atom( *itr ).element_type()->element() != core::chemical::element::H ) {
					bonded_heavy += 1;
				}
			}
			if ( bonded_heavy >= 2 ) {
				maxdist = maxdists[jj];
				nbr_atom = atom_vd;
			}
		}
	}

	if ( nbr_atom == ResidueType::null_vertex ) { // No suitable neighbor -- just pick atom 1
		TR.Warning << "No suitable neighbor atom found for " << res.name() << " -- picking first atom (" << res.atom_name(1) << ") instead." << std::endl;
		nbr_atom = res.atom_vertex( 1 );
		maxdist = maxdists[1];
	}
	return maxdist;
}

/// @brief Apply molfile_to_params style partial charges to the ResidueType.
/// @details These partial charges are based off of the Rosetta atom type,
/// adjusted such that the net partial charge is equal to the net formal charge.
///
/// These charges are almost certainly dodgy. If you have any other source of
/// partial charges that are at all reasonable, you probably want to consider those instead.
///
/// Assumes:
///   * All atoms and bond are present.
///   * All atom types have been set.
///   * Formal charges (if any) have been set.

void
rosetta_recharge_fullatom( ResidueType & res ) {
	ResidueGraph const & graph( res.graph() );
	AtomTypeSet const & ats( res.atom_type_set() );
	if ( ! ats.has_extra_parameter( "CHARGE" ) ) {
		TR.Warning << "Atom Type Set " << ats.name() << " is missing charging information - skipping recharging." << std::endl;
		return;
	}
	int charge_extra_param_index( ats.extra_parameter_index("CHARGE") );
	core::Real desired_net(0), current_net(0);
	core::Size natoms(0);
	VIter iter, iter_end;
	for ( boost::tie( iter, iter_end ) = boost::vertices( graph ); iter != iter_end; ++iter ) {
		Atom & atm( res.atom( *iter ) );
		desired_net += atm.formal_charge();
		core::Real charge(0);
		if ( atm.atom_type_index() != 0 ) {
			charge = ats[ atm.atom_type_index() ].extra_parameter( charge_extra_param_index );
		}
		atm.charge( charge );
		TR.Debug << "Residue " << res.name() << ": Charging atom " << atm.name() << " type " << ats[ atm.atom_type_index() ].name() << " at " << atm.charge() << std::endl;
		current_net += charge;
		if ( charge != 0 ) {
			natoms += 1;
		}
	}

	if ( current_net < (desired_net - 0.001) || current_net > (desired_net + 0.001) ) {
		core::Real correction( (desired_net-current_net)/natoms );
		for ( boost::tie( iter, iter_end ) = boost::vertices( graph ); iter != iter_end; ++iter ) {
			Atom & atm( res.atom( *iter ) );
			if ( atm.charge() != 0 ) {
				atm.charge( atm.charge() + correction );
				TR.Debug << "Residue " << res.name() << ": Adjusting charge on atom " << atm.name() << " type " << ats[ atm.atom_type_index() ].name() << " to " << atm.charge() << std::endl;
			}
		}
	}
}

/// @brief Make a (new) centroid version of the (fullatom) ResidueType passed in.
ResidueTypeOP
make_centroid( ResidueType const & res ) {

	AtomTypeSet const & old_ats( res.atom_type_set() );

	if ( old_ats.mode() == core::chemical::CENTROID_t ) {
		TR.Warning << "Residue " << res.name() << " to convert to centroid is already centroid!" << std::endl;
		return nullptr; // Don't bother - it's an error.
	} else if ( old_ats.mode() != core::chemical::FULL_ATOM_t ) {
		TR.Warning << "Residue " << res.name() << " to convert to centroid is not in full atom mode!" << std::endl;
	}

	ResidueTypeOP centroid( res.clone() );

	// Atom type translation from molfile_to_params.py
	// RM: There's been some additions to the centroid atom type set since
	// those might be a better idea ... or not.
	std::map< std::string, std::string > type_translation = {
		{"Hpol","HNbb"},
		{"Ntrp","Nbb" },
		{"NH2O","Nbb" },
		{"Nlys","Nbb" },
		{"Narg","Nbb" },
		{"Npro","Nbb" },
		{"Nhis","OCbb"}, // No, this is not a typo (apparently) - it's an unprotonated hbond acceptor
		{"OH"  ,"OCbb"},
		{"OW"  ,"OCbb"},
		{"ONH2","OCbb"},
		{"OOC" ,"OCbb"},
		{"Oaro","OCbb"},
		{"Hapo",""    }, // Don't translate.
		{"Haro",""    }, // Don't translate.
		{"aroC","CAbb"},
		{"CH1" ,"CAbb"},
		{"CH2" ,"CAbb"},
		{"CH3" ,"CAbb"},
		{"CNH2","CAbb"},
		{"COO" ,"CAbb"}};

	AtomTypeSetCOP centroid_ats_ptr( ChemicalManager::get_instance()->atom_type_set( core::chemical::CENTROID )  );
	AtomTypeSet const & centroid_ats( *centroid_ats_ptr );

	utility::vector1< std::string > new_types;
	for ( core::Size ii(1); ii <= centroid->natoms(); ++ii ) {
		core::Size old_index = res.atom(ii).atom_type_index();
		std::string const & old_string( old_ats[ old_index ].name() );
		std::string new_string;
		if ( type_translation.count(old_string) ) {
			new_string = type_translation[ old_string ];
		} else if ( centroid_ats.has_atom( old_string ) ) {
			// We have a simple as-is translation
			new_string = old_string;
		} else {
			TR.Warning << "Atom type '" << old_string << "' on atom '" << res.atom_name(ii)
				<< "' from residue type '" <<  res.name() <<  "' does not have a centroid mode equivalent: assuming CAbb." << std::endl;
			// This is how the molfile_to_params.py script does the translation for unrecognized atoms.
			new_string = "CAbb";
		}

		new_types.push_back( new_string );
	}

	centroid->set_atom_type_set( centroid_ats_ptr ); // Will zero out unknown types.

	utility::vector1< std::string > to_delete; // Don't modify ResidueType while iterating.
	for ( core::Size ii(1); ii <= centroid->natoms(); ++ii ) {
		std::string const & new_string( new_types[ii] );

		if ( new_string.empty() ) {
			to_delete.push_back( centroid->atom_name( ii ) );
		} else {
			debug_assert( centroid_ats.has_atom( new_string ) );
			// This should reset things properly, as we've set the atom type set to centroid previously.
			centroid->set_atom_type( centroid->atom_name( ii ), new_string );
		}
	}

	// Now delete the defered atoms
	for ( core::Size jj(1); jj <= to_delete.size(); ++jj ) {
		centroid->delete_atom( to_delete[ jj ] );
	}

	// We need to update the internal coordinates (to get rid of the missing hydrogens)
	// The neighbor settings should be fine - they're based on heavy atoms
	// Rotatable bonds should be fine - we shouldn't have rotatable bonds to apolar hydrogens
	// TODO: Charges might be a little funky, but I don't think we use them in centroid mode.
	if ( to_delete.size() ) {
		// Currently assign_internal_coordinates can't work with polymeric types.
		if ( res.n_possible_residue_connections() != 0 ) {
			TR.Warning << "Cannot automatically convert " << res.name() << " to centroid, as it is polymeric/has connections." << std::endl;
			return nullptr;
		}

		centroid->assign_internal_coordinates();
	}

	centroid->finalize();
	return centroid;
}

} // chemical
} // core
