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

#include <core/chemical/MutableResidueType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/ChemicalManager.hh>

// Project Headers
#include <utility/graph/Graph.hh>

#include <core/chemical/AtomType.hh>
#include <core/chemical/Atom.hh>
#include <core/chemical/AtomProperties.hh>
#include <core/chemical/Bond.hh>
#include <core/chemical/Elements.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/MutableResidueConnection.hh>
#include <utility/graph/RingDetection.hh>
#include <core/chemical/gasteiger/GasteigerAtomTypeData.hh>

#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <boost/unordered_map.hpp>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/graph/BFS_prune.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>

#include <boost/graph/breadth_first_search.hpp>

#include <algorithm>
#include <boost/unordered_map.hpp>

namespace core {
namespace chemical {

static basic::Tracer TR( "core.chemical.residue_support" );

ObjexxFCL::FArray2D_int
get_residue_path_distances( ResidueType const & res )
{
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

utility::vector1< VD >
mainchain_path( MutableResidueType const & res ) {
	if ( res.lower_connect_id() == 0 || res.upper_connect_id() == 0 ) {
		return utility::vector1< VD >{}; // Empty vector
	}

	return shortest_path( res, res.lower_connect_atom(), res.upper_connect_atom() );
}

utility::vector1< VD >
shortest_path( MutableResidueType const & res, VD start, VD end ) {
	debug_assert( has( res.graph(), start ) );
	debug_assert( has( res.graph(), end ) );

	typedef typename std::map<VD,VD> Predecessors;
	typedef typename boost::associative_property_map< Predecessors > PredecessorsMap;

	Predecessors predecessors;
	PredecessorsMap predecessors_map( predecessors );

	// During the BFS, record each nodes predecessor
	auto visitor = boost::make_bfs_visitor( boost::record_predecessors(predecessors_map,boost::on_tree_edge()) );

	// We start at the end, such that we can iterate predicessors, resulting in forward iteration
	boost::breadth_first_search( res.graph(), end, boost::visitor(visitor) );

	if ( predecessors.count( start ) == 0 ) {
		return utility::vector1< VD >{}; // Not found in BFS - no path.
	}

	utility::vector1< VD > path;
	VD current = start;
	while ( current != end ) {
		path.push_back( current );
		debug_assert( predecessors.count( current ) != 0 );
		current = predecessors[ current ];
	}
	path.push_back( end );

	return path;
}

void
annotate_backbone( MutableResidueType & restype ) {
	// First thing - reset all atoms to sidechain.
	for ( VD atm: restype.all_atoms() ) {
		restype.atom(atm).is_backbone( false );
	}

	utility::vector1< VD > to_process;

	// Now start with the seeds on the connects.
	if ( restype.lower_connect_id() != 0 ) {
		restype.atom( restype.lower_connect_atom() ).is_backbone( true );
		to_process.push_back( restype.lower_connect_atom() );
	}
	if ( restype.upper_connect_id() != 0 ) {
		restype.atom( restype.upper_connect_atom() ).is_backbone( true );
		to_process.push_back( restype.upper_connect_atom() );
	}

	if ( to_process.empty() ) { return; } // No starting points - all sidechain

	// Annotate rotatable bonds which we shouldn't cross.
	std::map< VD, utility::vector1< VD > > rot_nbr;
	for ( core::Size ii(1); ii <= restype.nchi(); ++ii ) {
		VDs const & chi_atoms = restype.chi_atom_vds( ii );
		rot_nbr[ chi_atoms[2] ].push_back( chi_atoms[3] ); // Center bond is the rotatable one
		rot_nbr[ chi_atoms[3] ].push_back( chi_atoms[2] );
	}

	// Flood fill atoms
	while ( ! to_process.empty() ) {
		VD atm = to_process.back();
		to_process.pop_back();
		for ( VD nbr: restype.bonded_neighbors(atm) ) {
			if ( restype.atom(nbr).is_backbone() ) { continue; } // Try to halt infinite recursion
			if ( rot_nbr[atm].contains( nbr ) ) { continue; } // Don't cross rotatable bonds
			if ( restype.bond( atm, nbr ).cut_bond() ) { continue; } // Don't cross cut bonds
			// If we made it here, we should be good.
			restype.atom( nbr ).is_backbone( true );
			to_process.push_back( nbr );
		}
	}
}

//////////////////////////////////////////////////////
/// Make all atoms virtual
/// @author Sebastian Rämisch <raemisch@scripps.edu>
void
real_to_virtual( MutableResidueType & restype ) {
	std::string const VIRT = "VIRT";
	for ( VD atm: restype.all_atoms() ) {
		restype.set_atom_type( atm, VIRT );
		restype.atom(atm).charge(0.0);
		restype.atom(atm).is_virtual( true );
	}
	restype.add_property("VIRTUAL_RESIDUE");
}

LightWeightResidueGraph convert_residuetype_to_light_graph( MutableResidueType const & res ){
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
rename_atoms( MutableResidueType & res, bool preserve/*=true*/ ) {
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
		res.rename_atom( *iter, name );
		res.remap_pdb_atom_names( true ); // We've renamed an atom, so we need to be flexible with PDB loading.
		name_counts[ name ]++;
	}
}

core::Real &
VDDistanceMatrix::operator() ( VD a, VD b ) {
	// Accessing matrix_[a] will force creation if it doesn't exist
	if ( ! matrix_[a].count(b) ) {
		matrix_[a][b] = default_;
	}
	return matrix_[a][b];
}

core::Real
VDDistanceMatrix::find_max_over( VD a ) {
	InternalVector const & submap( matrix_[a] );
	core::Real max_element = 0; // Distances are always greater than zero.
	for ( InternalVector::const_iterator iter( submap.begin() ), iter_end(submap.end()); iter != iter_end; ++iter ) {
		if ( iter->second > max_element ) {
			max_element = iter->second;
		}
	}
	return max_element;
}


/// @brief Utility visitor for find_nbr_dist
/// Will only traverse those atoms in the "rigid" portion of graph around the starting atom.
/// "Rigid" includes direct neighbors and atoms connected by non-rotatable bonds
/// e.g. all rings, all double/triple bonds, methyl groups, single atoms, etc.
class RigidDistanceVisitor: public utility::graph::null_bfs_prune_visitor {

public:
	RigidDistanceVisitor( VDDistanceMatrix & distances, MutableResidueType const & restype, VD start ) :
		distances_(distances),
		restype_(restype),
		start_(start),
		start_pos_( restype.atom(start).ideal_xyz() )
	{
		//start_index_ = restype.atom_index( start );
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
			//core::Size index( restype_.atom_index( target ) );
			distances_(start_, target) = restype_.atom( target ).ideal_xyz().distance( start_pos_ );
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
		//core::Size index( restype_.atom_index( vertex ) );
		distances_(start_,vertex) = restype_.atom( vertex ).ideal_xyz().distance( start_pos_ );
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
	VDDistanceMatrix & distances_;
	MutableResidueType const & restype_;
	VD start_;
	//core::Size start_index_;
	Vector const & start_pos_;
};

/// @brief Calculate the rigid matrix - assume that distances has been initialized to some really large value, and is square
void calculate_rigid_matrix( MutableResidueType const & res, VDDistanceMatrix & distances ) {
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
	VIter kk, kk_end, jj, jj_end, ii, ii_end;
	for ( boost::tie(kk, kk_end) = res.atom_iterators(); kk != kk_end; ++kk ) {
		for ( boost::tie(jj, jj_end) = res.atom_iterators(); jj != jj_end; ++jj ) {
			for ( boost::tie(ii, ii_end) = res.atom_iterators(); ii != ii_end; ++ii ) {
				core::Real new_dist( distances(*ii,*kk) + distances(*kk,*jj) );
				if ( new_dist < distances(*ii,*jj) ) {
					distances(*ii,*jj) = new_dist;
				}
			}
		}
	}

	if ( TR.Trace.visible() ) {
		// Print out the full distance matrix for debugging purposes.
		TR.Trace << std::setprecision(4);
		TR.Trace << "Atom distance matrix for " << res.name() << std::endl;
		TR.Trace << "    " << '\t';
		for ( VD xx: res.all_atoms() ) {
			TR.Trace << res.atom_name(xx) << '\t';
		}
		TR.Trace << std::endl;
		for ( VD yy: res.all_atoms() ) {
			TR.Trace << res.atom_name(yy) << '\t';
			for ( VD xx: res.all_atoms() ) {
				TR.Trace << distances(yy,xx) << '\t';
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
find_nbr_dist( MutableResidueType const & res, VD & nbr_atom ) {
	if ( res.natoms() == 0 ) {
		utility_exit_with_message("Cannot find neighbor atom distance for empty residue type.");
	}
	core::Real maxdist = 1e9; // Hopefully sufficiently large.
	VDDistanceMatrix distances( maxdist );
	calculate_rigid_matrix( res, distances );

	// TODO: Although we throw out hydrogens as potential neighbor atoms,
	// I believe the meaning of neighbor atoms is heavy-atom distances
	// We could get smaller values by removing hydrogens from distance considerations.

	boost::unordered_map< VD, core::Real > maxdists;
	VIter ii, ii_end;
	for ( boost::tie(ii, ii_end) = res.atom_iterators(); ii != ii_end; ++ii ) {
		maxdists[*ii] = distances.find_max_over( *ii );
	}

	if ( nbr_atom != MutableResidueType::null_vertex ) {
		if ( maxdists.count( nbr_atom ) ) {
			// We have an atom in mind - we just need the distance.
			if ( res.atom(nbr_atom).element_type() && res.atom(nbr_atom).element_type()->element() == element::H ) {
				TR.Warning << "Specified neighbor atom is a Hydrogen!!." << std::endl;
			}
			return maxdists[ nbr_atom ];
		} else {
			TR.Warning << "Specified neighbor atom is not present in the distance calculation - recomputing." << std::endl;
			nbr_atom = MutableResidueType::null_vertex;
		}
	}

	//maxdist initialized above
	for ( VD atom_vd: res.all_atoms() ) {
		if ( maxdists[atom_vd] < maxdist &&
				res.atom(atom_vd).element_type()->element() != core::chemical::element::H ) { // not hydrogen
			//Use graph directly, as other neighbor annotations may not be fully baked
			core::Size bonded_heavy(0);
			AdjacentIter itr, itr_end;
			for ( boost::tie(itr, itr_end) = res.bonded_neighbor_iterators(atom_vd); itr != itr_end; ++itr ) {
				if ( res.atom( *itr ).element_type()->element() != core::chemical::element::H ) {
					bonded_heavy += 1;
				}
			}
			if ( bonded_heavy >= 2 ) {
				maxdist = maxdists[atom_vd];
				nbr_atom = atom_vd;
			}
		}
	}

	if ( nbr_atom == MutableResidueType::null_vertex ) { // No suitable neighbor -- just pick any atom, preferring a heavy atom
		nbr_atom = maxdists.begin()->first;
		maxdist = maxdists.begin()->second;
		for ( boost::unordered_map< VD, core::Real >::const_iterator iter(maxdists.begin()), iter_end(maxdists.end()); iter != iter_end; ++iter ) {
			if ( res.atom(iter->first).element_type() && res.atom(iter->first).element_type()->element() != element::H ) {
				nbr_atom = iter->first;
				maxdist = iter->second;
				break;
			}
		}
		TR.Warning << "No suitable neighbor atom found for " << res.name() << " -- picking atom (" << res.atom_name(nbr_atom) << ") instead." << std::endl;
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
rosetta_recharge_fullatom( MutableResidueType & res ) {
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

MutableResidueTypeOP
make_centroid( ResidueType const & res ) {
	return make_centroid( MutableResidueType( res ) ); // I'm not 100% on the implicit double copy here, but it's easier
}

MutableResidueTypeOP
make_centroid( MutableResidueType const & res ) {

	AtomTypeSet const & old_ats( res.atom_type_set() );

	if ( old_ats.mode() == core::chemical::CENTROID_t ) {
		TR.Warning << "Residue " << res.name() << " to convert to centroid is already centroid!" << std::endl;
		return nullptr; // Don't bother - it's an error.
	} else if ( old_ats.mode() != core::chemical::FULL_ATOM_t ) {
		TR.Warning << "Residue " << res.name() << " to convert to centroid is not in full atom mode!" << std::endl;
	}

	MutableResidueTypeOP centroid( new MutableResidueType( res ) );

	// Atom type translation from molfile_to_params.py
	// RM: There's been some additions to the centroid atom type set since
	// those might be a better idea ... or not.
	std::map< std::string, std::string > type_translation = {
		{"Hpol","HNbb"},
		{"Ntrp","Nbb" },
		{"NH2O","Nbb" },
		{"Nlys","Nbb" },
		{"Narg","Nbb" },
		{"NtrR","Nbb" },
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
		{"CH0", "CAbb"},
		{"CH1" ,"CAbb"},
		{"CH2" ,"CAbb"},
		{"CH3" ,"CAbb"},
		{"CNH2","CAbb"},
		{"COO" ,"CAbb"}};

	AtomTypeSetCOP centroid_ats_ptr( ChemicalManager::get_instance()->atom_type_set( core::chemical::CENTROID )  );
	AtomTypeSet const & centroid_ats( *centroid_ats_ptr );

	std::map< VD, std::string > new_types;
	for ( VD atm: centroid->all_atoms() ) {
		// We assume that the order is the same between the two
		core::Size old_index = res.atom( res.all_atoms()[ centroid->atom_index(atm) ] ).atom_type_index();
		std::string const & old_string( old_ats[ old_index ].name() );
		std::string new_string;
		if ( type_translation.count(old_string) ) {
			new_string = type_translation[ old_string ];
		} else if ( centroid_ats.has_atom( old_string ) ) {
			// We have a simple as-is translation
			new_string = old_string;
		} else {
			TR.Warning << "Atom type '" << old_string << "' on atom '" << centroid->atom_name(atm)
				<< "' from residue type '" <<  centroid->name() <<  "' does not have a centroid mode equivalent: assuming CAbb." << std::endl;
			// This is how the molfile_to_params.py script does the translation for unrecognized atoms.
			new_string = "CAbb";
		}

		new_types[ atm ] = new_string;
	}

	centroid->update_atom_type_set( centroid_ats_ptr ); // Will zero out unknown types.

	utility::vector1< std::string > to_delete; // Don't modify ResidueType while iterating.
	for ( VD atm: centroid->all_atoms() ) {
		std::string const & new_string( new_types[ atm ] );

		if ( new_string.empty() ) {
			to_delete.push_back( centroid->atom_name( atm ) );
		} else {
			debug_assert( centroid_ats.has_atom( new_string ) );
			// This should reset things properly, as we've set the atom type set to centroid previously.
			centroid->set_atom_type( centroid->atom_name( atm ), new_string );
		}
	}

	// Now delete the defered atoms
	for ( core::Size jj(1); jj <= to_delete.size(); ++jj ) {
		centroid->delete_atom( to_delete[ jj ] );
	}

	if ( ! to_delete.empty() ) {
		// We may have invalidated some of the internal coordinates/chis, due to the deletion of the hydrogens.
		// It's unlikely, but check if that's the case.
		// TODO: Charges might be a little funky, but I don't think we use them in centroid mode.

		// Check if we have any rotatable bonds defined in terms of missing atoms
		for ( core::Size ii(1); ii <= centroid->nchi(); ++ii ) {
			if ( ! centroid->chi_valid( ii ) ) {
				TR.Warning << "Cannot automatically convert " << res.name() << " to centroid, as it has rotatable bond which depend on now-missing atoms." << std::endl;
				return nullptr;
			}
		}

		// Check if we have any ICOOR that still depend on now-missing atoms.
		bool needs_icoor_update = false;
		for ( VD atm: centroid->all_atoms() ) {
			MutableICoorRecordCOP ic = centroid->atom( atm ).icoor();
			if ( ic == nullptr ) {
				needs_icoor_update = true;
				break;
			}
			for ( core::Size stub(1); stub <= 3; ++stub ) {
				if ( ic->stub_type(stub) == ICoordAtomIDType::INTERNAL && ! centroid->has( ic->stub_atom(stub) ) ) {
					needs_icoor_update = true;
					break;
				}
			}
		}

		// Need to check connections, too.
		// They *should* be based off of heavy atoms, but potentially not.
		for ( core::Size cc(1); cc <= centroid->n_possible_residue_connections(); ++cc ) {
			MutableResidueConnection const & connect( centroid->residue_connection(cc) );
			if ( ! centroid->has( connect.vertex() ) ) {
				needs_icoor_update = true;
				break;
			}
			MutableICoorRecord const & icoor( connect.icoor() );
			for ( core::Size stub(1); stub <= 3; ++stub ) {
				if ( icoor.stub_type(stub) == ICoordAtomIDType::INTERNAL && ! centroid->has( icoor.stub_atom(stub) ) ) {
					needs_icoor_update = true;
					break;
				}
			}
		}

		if ( needs_icoor_update ) {
			// Currently assign_internal_coordinates can't work with polymeric types.
			if ( res.n_possible_residue_connections() != 0 ) {
				TR.Warning << "Cannot automatically convert " << res.name() << " to centroid, as it needs ICoor reassignment and it is polymeric/has connections." << std::endl;
				return nullptr;
			}
			// Fix up all the internal icoords.
			TR.Debug << "Reassiging internal coordinates due to missing atoms." << std::endl;
			centroid->assign_internal_coordinates();
		}
	}

	return centroid;
}

bool
residue_type_bases_identical( ResidueTypeBase const & r1, ResidueTypeBase const & r2 ) {
	if ( r1.mode() != r2.mode() ) return false;
	// Assume that we need *exactly* the same type sets to be equal.
	if ( r1.atom_type_set_ptr() != r2.atom_type_set_ptr() ) return false;
	if ( r1.element_set_ptr() != r2.element_set_ptr() ) return false;
	if ( r1.mm_atom_types_ptr() != r2.mm_atom_types_ptr() ) return false;
	if ( r1.gasteiger_atom_typeset() != r2.gasteiger_atom_typeset() ) return false;
	if ( r1.orbital_types_ptr() != r2.orbital_types_ptr() ) return false;

	if ( r1.aa() != r2.aa() ) return false;
	if ( r1.backbone_aa() != r2.backbone_aa() ) return false;
	//if ( r1.rotamer_aa() != r2.rotamer_aa() ) return false;
	if ( r1.na_analogue() != r2.na_analogue() ) return false;
	if ( r1.base_analogue() != r2.base_analogue() ) return false;

	if ( r1.base_name() != r2.base_name() ) return false;
	if ( r1.name() != r2.name() ) return false;
	if ( r1.name1() != r2.name1() ) return false;
	if ( r1.name3() != r2.name3() ) return false;
	if ( r1.interchangeability_group() != r2.interchangeability_group() ) return false;

	if ( r1.get_rama_prepro_mainchain_torsion_potential_name(false) != r2.get_rama_prepro_mainchain_torsion_potential_name(false) ) return false;
	if ( r1.get_rama_prepro_mainchain_torsion_potential_name(true) != r2.get_rama_prepro_mainchain_torsion_potential_name(true) ) return false;
	if ( r1.get_rama_prepro_map_file_name(false) != r2.get_rama_prepro_map_file_name(false) ) return false;
	if ( r1.get_rama_prepro_map_file_name(true) != r2.get_rama_prepro_map_file_name(true) ) return false;

	if ( r1.net_formal_charge() != r2.net_formal_charge() ) return false;
	if ( r1.force_nbr_atom_orient() != r2.force_nbr_atom_orient() ) return false;
	if ( r1.remap_pdb_atom_names() != r2.remap_pdb_atom_names() ) return false;
	if ( r1.is_metapatched() != r2.is_metapatched() ) return false;

	if ( r1.properties() != r2.properties() ) return false;
	// TODO: Probably should check rotamer library specification.
	//if ( *r1.rotamer_library_specification() != *r2.rotamer_library_specification() ) return false;

	if ( r1.get_metal_binding_atoms() != r2.get_metal_binding_atoms() ) return false;
	if ( r1.get_disulfide_atom_name() != r2.get_disulfide_atom_name() ) return false;

	if ( r1.atom_aliases() != r2.atom_aliases() ) return false;
	if ( r1.canonical_atom_aliases() != r2.canonical_atom_aliases() ) return false;

	if ( r1.defined_adducts() != r2.defined_adducts() ) return false;
	if ( r1.n_orbitals() != r2.n_orbitals() ) return false;
	// TODO: Probably should check orbitals themselves, but don't bother.

	return true;
}

bool
residue_types_identical( ResidueType const & r1, ResidueType const & r2 ) {
	if ( ! residue_type_bases_identical( r1, r2 ) ) return false;

	if ( r1.natoms() != r2.natoms() ) return false;
	if ( r1.nbonds() != r2.nbonds() ) return false;
	if ( r1.nheavyatoms() != r2.nheavyatoms() ) return false;
	if ( r1.n_hbond_acceptors() != r2.n_hbond_acceptors() ) return false;
	if ( r1.n_hbond_donors() != r2.n_hbond_donors() ) return false;
	if ( r1.last_backbone_atom() != r2.last_backbone_atom() ) return false;
	if ( r1.first_sidechain_hydrogen() != r2.first_sidechain_hydrogen() ) return false;

	//if ( r1.rotamer_aa() != r2.rotamer_aa() ) return false;

	// atom_name_to_index_ ??
	for ( core::Size ii(1); ii<= r1.natoms(); ++ii ) {
		if ( r1.atom_name(ii) != r2.atom_name(ii) ) return false;
		if ( r1.atom_type_index(ii) != r2.atom_type_index(ii) ) return false;
		//if ( r1.atom_type(ii) != r2.atom_type(ii) ) ) return false;
		if ( r1.element(ii) != r2.element(ii) ) return false;
		if ( r1.mm_atom_type_index(ii) != r2.mm_atom_type_index(ii) ) return false;
		if ( r1.gasteiger_atom_type(ii) != r2.gasteiger_atom_type(ii) ) return false;
		if ( r1.formal_charge(ii) != r2.formal_charge(ii) ) return false;
		if ( r1.atom_charge(ii) != r2.atom_charge(ii) ) return false;
		if ( r1.ideal_xyz(ii) != r2.ideal_xyz(ii) ) return false; // Potentially an issue???

		if ( r1.atom_properties(ii) != r2.atom_properties(ii) ) return false;
		if ( ! compare_atom_icoor( r1.icoor(ii), r2.icoor(ii), true ) ) return false;

		if ( ii <= r1.nheavyatoms() ) {
			if ( r1.attached_H_begin(ii) != r2.attached_H_begin(ii) ) return false;
			if ( r1.attached_H_end(ii) != r2.attached_H_end(ii) ) return false;
		}

		// Shouldn't have ordering issues (should all be in numeric order.)
		if ( r1.bonded_neighbor(ii) != r2.bonded_neighbor(ii) ) return false;
		if ( r1.bonded_neighbor_types(ii) != r2.bonded_neighbor_types(ii) ) return false;
		if ( r1.bonded_neighbor_ringnesses(ii) != r2.bonded_neighbor_ringnesses(ii) ) return false;

		if ( r1.dihedrals_for_atom(ii) != r2.dihedrals_for_atom(ii) ) return false;
		if ( r1.bondangles_for_atom(ii) != r2.bondangles_for_atom(ii) ) return false;

		if ( r1.last_controlling_chi(ii) != r2.last_controlling_chi(ii) ) return false;

		if ( r1.bonded_orbitals(ii) != r2.bonded_orbitals(ii) ) return false;

		if ( r1.heavyatom_has_polar_hydrogens(ii) != r2.heavyatom_has_polar_hydrogens(ii) ) return false;

		if ( r1.atom_base(ii) != r2.atom_base(ii) ) return false;
		if ( r1.abase2(ii) != r2.abase2(ii) ) return false;

		if ( r1.cut_bond_neighbor(ii) != r2.cut_bond_neighbor(ii) ) return false;
		if ( r1.atom_being_shadowed(ii) != r2.atom_being_shadowed(ii) ) return false;

		if ( r1.atom_depends_on_lower_polymeric_connection(ii) != r2.atom_depends_on_lower_polymeric_connection(ii) ) return false;
		if ( r1.atom_depends_on_upper_polymeric_connection(ii) != r2.atom_depends_on_upper_polymeric_connection(ii) ) return false;

	}

	if ( r1.ndihe() != r2.ndihe() ) return false;
	for ( core::Size dd(1); dd <= r1.ndihe(); ++dd ) {
		if ( r1.dihedral(dd) != r2.dihedral(dd) ) return false;
	}
	if ( r1.num_bondangles() != r2.num_bondangles() ) return false;
	for ( core::Size ba(1); ba <= r1.num_bondangles(); ++ba ) {
		if ( r1.bondangle(ba) != r2.bondangle(ba) ) return false;
	}

	if ( r1.atoms_with_orb_index() != r2.atoms_with_orb_index() ) return false;

	if ( r1.Haro_index() != r2.Haro_index() ) return false;
	if ( r1.Hpol_index() != r2.Hpol_index() ) return false;
	if ( r1.accpt_pos() != r2.accpt_pos() ) return false;
	if ( r1.Hpos_polar() != r2.Hpos_polar() ) return false;
	if ( r1.Hpos_apolar() != r2.Hpos_apolar() ) return false;
	if ( r1.accpt_pos_sc() != r2.accpt_pos_sc() ) return false;
	if ( r1.Hpos_polar_sc() != r2.Hpos_polar_sc() ) return false;
	if ( r1.all_bb_atoms() != r2.all_bb_atoms() ) return false;
	if ( r1.all_sc_atoms() != r2.all_sc_atoms() ) return false;

	if ( r1.nchi() != r2.nchi() ) return false;
	if ( r1.chi_atoms() != r2.chi_atoms() ) return false;

	for ( core::Size cc(1); cc <= r1.nchi(); ++cc ) {
		if ( r1.is_proton_chi(cc) != r2.is_proton_chi(cc) ) return false;
		if ( r1.chi_2_proton_chi(cc) != r2.chi_2_proton_chi(cc) ) return false;
		if ( r1.chi_rotamers(cc) != r2.chi_rotamers(cc) ) return false;
		if ( r1.atoms_last_controlled_by_chi(cc) != r2.atoms_last_controlled_by_chi(cc) ) return false;
	}
	for ( core::Size cc(1); cc <= r1.n_proton_chi(); ++cc ) {
		if ( r1.proton_chi_2_chi( cc ) != r2.proton_chi_2_chi( cc ) ) return false;
		if ( r1.proton_chi_samples(cc) != r2.proton_chi_samples(cc) ) return false;
		if ( r1.proton_chi_extra_samples(cc) != r2.proton_chi_extra_samples(cc) ) return false;
	}

	if ( r1.nu_atoms() != r2.nu_atoms() ) return false;
	if ( r1.ring_atoms() != r2.ring_atoms() ) return false;
	if ( r1.n_ring_conformer_sets() != r2.n_ring_conformer_sets() ) return false;
	for ( core::Size rr(1); rr <= r1.n_ring_conformer_sets(); ++rr ) {
		if ( r1.ring_conformer_set(rr) != r2.ring_conformer_set(rr) ) {
			if ( r1.ring_conformer_set(rr) == nullptr ) return false; // If one is nullptr, they both should be
			if ( r2.ring_conformer_set(rr) == nullptr ) return false;
			// The Ring conformer sets are created deterministically from other data - skip for now
			// if ( *r1.ring_conformer_set(rr) != *r2.ring_conformer_set(rr) ) return false;
		}
	}
	if ( r1.lowest_ring_conformers() != r2.lowest_ring_conformers() ) return false;
	if ( r1.low_ring_conformers() != r2.low_ring_conformers() ) return false;

	if ( r1.path_distances() != r2.path_distances() ) return false;

	if ( r1.root_atom() != r2.root_atom() ) return false;
	if ( r1.nbr_atom() != r2.nbr_atom() ) return false;
	if ( r1.nbr_radius() != r2.nbr_radius() ) return false;
	if ( r1.mass() != r2.mass() ) return false;

	if ( r1.n_possible_residue_connections() != r2.n_possible_residue_connections() ) return false;
	for ( core::Size rr(1); rr <= r1.n_possible_residue_connections(); ++ rr ) {
		if ( ! compare_residue_connection( r1.residue_connection(rr), r2.residue_connection(rr), true) ) return false;
		for ( core::Size ii(1); ii<= r1.natoms(); ++ii ) {
			if ( r1.atom_depends_on_connection(ii,rr) != r2.atom_depends_on_connection(ii,rr) ) return false;
		}
	}
	// Bonds within distance of connection-- skip for now
	if ( r1.lower_connect_id() != r2.lower_connect_id() ) return false;
	if ( r1.upper_connect_id() != r2.upper_connect_id() ) return false;
	if ( r1.n_non_polymeric_residue_connections() != r2.n_non_polymeric_residue_connections() ) return false;
	if ( r1.n_polymeric_residue_connections() != r2.n_polymeric_residue_connections() ) return false;

	if ( r1.get_base_type_cop().get() == &r1 ) {
		if ( r2.get_base_type_cop().get() != &r2 ) return false; // Self-base
	} else if ( r1.get_base_type_cop() != r2.get_base_type_cop() ) {
		if ( !residue_types_identical(*r1.get_base_type_cop(),*r2.get_base_type_cop()) ) return false;
	}

	// RNA and Carbohydrate info is created deterministicly from other data - skip for now
	//if ( r1.RNA_info() != r2.RNA_info() ) return false;
	//if ( r1.carbohydrate_info() != r2.carbohydrate_info() || *r1.carbohydrate_info() != *r2.carbohydrate_info() ) return false;

	if ( r1.mainchain_atoms() != r2.mainchain_atoms() ) return false;
	if ( r1.actcoord_atoms() != r2.actcoord_atoms() ) return false;

	if ( r1.has_polymer_dependent_groups() != r2.has_polymer_dependent_groups() ) return false;

	return true;
}

bool
compare_residue_connection( ResidueConnection const & rc1, ResidueConnection const & rc2, bool fuzzy) {
	if ( rc1.atomno() != rc2.atomno() ) return false;
	if ( rc1.index() != rc2.index() ) return false;
	return compare_atom_icoor( rc1.icoor(), rc2.icoor(), fuzzy );
}

bool
compare_atom_icoor( AtomICoor const & aic1, AtomICoor const & aic2, bool fuzzy ) {
	if ( aic1.built_atom() != aic2.built_atom() ) return false;
	if ( aic1.stub_atom1() != aic2.stub_atom1() ) return false;
	if ( aic1.stub_atom2() != aic2.stub_atom2() ) return false;
	if ( aic1.stub_atom3() != aic2.stub_atom3() ) return false;
	if ( !fuzzy ) {
		if ( aic1.d() != aic2.d() ) return false;
		if ( aic1.theta() != aic2.theta() ) return false;
		if ( aic1.phi() != aic2.phi() ) return false;
	} else {
		constexpr core::Real delta = 1e-10; // Enough to accomodate for slight noise from building, not large enough to be an appreciable different.
		if ( aic1.d() > aic2.d() + delta ) return false;
		if ( aic1.d() < aic2.d() - delta ) return false;
		if ( aic1.theta() > aic2.theta() + delta ) return false;
		if ( aic1.theta() < aic2.theta() - delta ) return false;
		if ( aic1.phi() > aic2.phi() + delta ) return false;
		if ( aic1.phi() < aic2.phi() - delta ) return false;
	}
	return true;
}

} // chemical
} // core
