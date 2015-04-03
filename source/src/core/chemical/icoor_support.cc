// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/icoor_support.cc
/// @brief  external support for manipulating the ResidueType's icoor/ideal_xyz
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/icoor_support.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/ResidueConnection.hh>

#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/tree/BondedAtom.hh>
#include <core/kinematics/tree/JumpAtom.hh>

#include <core/id/AtomID.hh>
//#include <core/id/types.hh>

#include <core/types.hh>

#include <basic/Tracer.hh>
#include <numeric/xyzVector.io.hh>

#include <utility/graph/DFS_sort.hh>

//#include <boost/graph/depth_first_search.hpp>
//#include <boost/graph/properties.hpp>

#include <set>
#include <queue>

namespace core {
namespace chemical {

static thread_local basic::Tracer TR( "core.chemical.icoor_support" );

class RerootRestypeVisitor: public boost::default_dfs_visitor {
private:
	core::chemical::VD root_;
	VdTreeatomMap & treeatom_map_;
	core::chemical::ResidueType const & restype_;
public:
	RerootRestypeVisitor( VdTreeatomMap & map, core::chemical::VD root, core::chemical::ResidueType const & restype ):
		root_( root ),
		treeatom_map_( map ),
		restype_( restype )
	{ }

	template <class Vertex, class Graph>
	void initialize_vertex(Vertex u, const Graph& /*g*/) {
		if( treeatom_map_.count(u) == 0 ) {
			treeatom_map_[ u ] = core::kinematics::tree::AtomOP( new core::kinematics::tree::BondedAtom );
		}
		treeatom_map_[ u ]->xyz( restype_.atom( u ).ideal_xyz() );
		// Fill in vertex identity for debugging purposes:
		core::id::AtomID id(restype_.atom_index(u), 1);
		treeatom_map_[ u ]->id( id );
	}
	template <class Vertex, class Graph>
	void start_vertex(Vertex u, const Graph& /*g*/) {
		// The DFS algorithm will restart on disconnected subgraphs --
		// we currently can't handle the disconnected portion of it,
		// so raise an error if that happens.
		// (This means the root needs to be pre-loaded into the treeatom_map
	debug_assert( restype_.has( u ) );
		if( u != root_ ) {
			TR.Error << "ERROR: For ResidueType " << restype_.name() << ", atoms "
					<< restype_.atom_name( root_ ) << " and " << restype_.atom_name( u )
					<< " are not connected via bonds." << std::endl;
			utility_exit_with_message( "Cannot reroot a disconnected ResidueType: ");
		}
	}
	template <class Edge, class Graph>
	void tree_edge(Edge u, const Graph& g) {
		core::chemical::VD parent( boost::source(u,g) );
		core::chemical::VD child( boost::target(u,g) );
	debug_assert( treeatom_map_.count( parent ) );
	debug_assert( treeatom_map_.count( child ) );
		treeatom_map_[ parent ]->append_atom( treeatom_map_[ child ] );
	}
};

/// @brief Edge sorting:
/// Return true if we should prefer edge1 over edge2
/// * Non-cut bonds before cut bonds
/// * Actual bonds before pseudobonds
/// * Bonds to Heavy atoms before light atoms
/// * Bonds to Concrete atoms before virtual atoms
///
/// This doesn't (need to?) quite match the logic in core/conformation/util.cc:setup_atom_links()
class RerootEdgeSorter {
public:
	RerootEdgeSorter(core::chemical::ResidueGraph const & graph, core::chemical::ResidueType const & restype):
		graph_(graph),
		restype_(restype)
	{}

	/// Return true if the first argument goes before the second argument
	bool operator()(core::chemical::ED edge1, core::chemical::ED edge2 ) {
		core::chemical::Bond const & bond1( graph_[edge1] );
		core::chemical::Bond const & bond2( graph_[edge2] );

		/// * Non-cut bonds before cut bonds
		if( bond1.cut_bond() && ! bond2.cut_bond() ) {
			return false;
		} else if ( bond2.cut_bond() && ! bond1.cut_bond() ) {
			return true;
		}

		// * Actual bonds before pseudobonds - TODO better pseudobond testing
		if( bond1.bond_name() == core::chemical::UnknownBond && bond2.bond_name() != core::chemical::UnknownBond) {
			return false;
		} else if ( bond1.bond_name() != core::chemical::UnknownBond && bond2.bond_name() == core::chemical::UnknownBond) {
			return true;
		}
	debug_assert( boost::source(edge1,graph_) == boost::source(edge2,graph_) );
		//core::chemical::VD source( boost::source(edge1,graph_) );
		core::chemical::VD target1( boost::target(edge1,graph_) );
		core::chemical::VD target2( boost::target(edge2,graph_) );
		core::chemical::Atom const & atom1( graph_[ target1 ] );
		core::chemical::Atom const & atom2( graph_[ target2 ] );

/*
		// * Rotatable bonds before non-rotatable bonds.
		// This is slightly manky - should we pre annotate the bonds as rotatable?
		bool b1rot(false), b2rot(false);
		for ( core::Size chino(1); chino <= restype_.nchi(); ++chino ) {
			core::chemical::AtomIndices const & ai( restype_.chi_atoms(chino) );
			if( (ai[2] == restype_.atom_index(source) && ai[3] == restype_.atom_index(target1)) ||
					(ai[3] == restype_.atom_index(source) && ai[2] == restype_.atom_index(target1)) ) {
				b1rot = true;
				TR << "Rotatable: " << restype_.atom_name(source) << " -- " << restype_.atom_name(target1) << std::endl;
			}
			if( (ai[2] == restype_.atom_index(source) && ai[3] == restype_.atom_index(target2)) ||
					(ai[3] == restype_.atom_index(source) && ai[2] == restype_.atom_index(target2)) ) {
				b2rot = true;
				TR << "Rotatable: " << restype_.atom_name(source) << " -- " << restype_.atom_name(target2) << std::endl;
			}
		}
		if( ! b1rot && b2rot ) {
			return false;
		} else if (b1rot && ! b2rot ) {
			return true;
		}
*/

		// * Bonds to Heavy atoms before light atoms
		if( atom1.is_hydrogen() && ! atom2.is_hydrogen() ) {
			return false;
		} else if ( ! atom1.is_hydrogen() && atom2.is_hydrogen() ) {
			return true;
		}
		/// * Bonds to Concrete atoms before virtual atoms
		if( atom1.is_virtual() && ! atom2.is_virtual() ) {
			return false;
		} else if( ! atom1.is_virtual() && atom2.is_virtual() ) {
			return true;
		}
		/// For reproducible ordering, the "smaller" atom name should be first
		return atom1.name() < atom2.name();
	}
private:
	core::chemical::ResidueGraph const & graph_;
	core::chemical::ResidueType const & restype_;
};


/// @brief Reroot the Icoord records of a ResidueType on the given atom
/// We need direct access to the ResidueGraph, so this function can only be called by
/// ResidueType itself
///
/// @details Doing a depth first search here because that's what molfile_to_params.py does:
/// "Protein residues appear to go depth first, so that all chi angles ride on each other."
/// RM: Doing a breadth first search would likely result in a shallower tree, but with possibly
/// different behavior on how ring atom trees are built.
///
/// Note that updating the ICOOR records will also update the atom_base values.
///
/// Assumes:
/// * All bonds and atoms exist in the graph,
/// * The graph is completely connected.
/// * All ideal xyz coordinates are updated.
///
void
reroot_restype( core::chemical::ResidueType & restype, core::chemical::ResidueGraph const & graph, core::chemical::VD root) {
	if( restype.natoms() < 3 ) {
		TR.Warning << "Warning: Cannot re-root residue type with less than three atoms." << std::endl;
		return;
	}
	VdTreeatomMap treeatom_map;
	std::map< core::chemical::VD, boost::default_color_type > colormap;
	boost::associative_property_map< std::map< core::chemical::VD, boost::default_color_type > > color_property_map(colormap);

	// The root of the atom tree needs to be a jump atom,
	// as opposed to the Bonded atom that would be initialized by default in the DFS
	core::kinematics::tree::AtomOP rootatom( new core::kinematics::tree::JumpAtom );
	treeatom_map[ root ] = rootatom;
	TR << "Rooting on atom: " << restype.atom_name(root) << std::endl;

	RerootRestypeVisitor visitor( treeatom_map, root, restype );

	utility::graph::depth_first_search_sort(graph, visitor, color_property_map, root, RerootEdgeSorter( graph, restype ) );

	rootatom->update_internal_coords( true );

	// Reverse the map for easy lookup
	std::map< core::kinematics::tree::AtomCOP, core::chemical::VD > revmap;
	for( VdTreeatomMap::const_iterator itr(treeatom_map.begin()), itr_end(treeatom_map.end()); itr != itr_end; ++itr) {
		revmap[ itr->second ] = itr->first;
	}

	// Now pull the icoord data out of the kinematic atoms

	core::chemical::VIter iter, iter_end;
	for( boost::tie(iter, iter_end) = boost::vertices(graph); iter != iter_end; ++iter ) {
		core::chemical::VD const & atomVD( *iter );
	debug_assert( treeatom_map.count( atomVD ) );

		core::kinematics::tree::AtomCOP atom = treeatom_map[ atomVD ];

		core::kinematics::tree::AtomCOP parent, angle, torsion;
		core::Real phi, theta, d;
		parent = atom->parent();
		if( parent ) {
			// Regular, non-root atom
		debug_assert( utility::pointer::dynamic_pointer_cast< core::kinematics::tree::BondedAtom const >( atom ) );

			// parent = atom->input_stub_atom1();
			angle = atom->input_stub_atom2();
			torsion = atom->input_stub_atom3();

			phi = atom->dof(core::id::PHI);
			theta = atom->dof(core::id::THETA);
			d = atom->dof(core::id::D);
		} else {
			// This is the root atom.
			if( atom->n_children() < 1 ) {
				utility_exit_with_message("ERROR: Rerooting restype for root with no children.");
			}
			// The children use zero based indices for some reason ...
			parent = atom;
			angle = atom->child(0);
			// The following logic is taken from molfile_to_params
			// Define the torsion progressively, if possible, else use an improper through the root atoms' other child.
			if( angle->n_children() > 0 ) {
				torsion = angle->child(0);
			} else {
				if( atom->n_children() < 2 ) {
					utility_exit_with_message("ERROR: Rerooting restype for root with only one child and no grandchildren.");
				}
				torsion = atom->child(1);
			}
			phi = 0;
			theta = 0;
			d = 0;
		}

		core::chemical::VD parentVD( revmap[parent] ), angleVD( revmap[angle] ), torsionVD( revmap[torsion] );
	debug_assert( restype.has( parentVD ) && restype.has( angleVD ) && restype.has( torsionVD ) );

		// Note: set_icoor will automatically setup the atom base for us.
		restype.set_icoor( atomVD, phi, theta, d, parentVD, angleVD, torsionVD );
	}
}

/// @brief Utility function for fill_ideal_xyz_from_icoor() -- does this ICoorAtomID have all the dependancies filled?
bool has_assigned_coords(ICoorAtomID const & stub, std::set< VD > const & assigned, core::chemical::ResidueType const & restype) {
	if( stub.type() == ICoorAtomID::INTERNAL ) {
	debug_assert( restype.has( stub.vertex() ) );
		return assigned.count( stub.vertex() );
	} else {
		// For connections, they have assigned coords if all their dependancies have assigned coords.
		ResidueConnection connection;
		if( stub.type() == ICoorAtomID::POLYMER_LOWER ) {
			connection = restype.lower_connect();
		} else if( stub.type() == ICoorAtomID::POLYMER_UPPER ) {
			connection = restype.upper_connect();
		} else if( stub.type() == ICoorAtomID::CONNECT ) {
			connection = restype.residue_connection( stub.atomno() );
		} else {
			utility_exit_with_message("Unable to assign coordinates for "+restype.name()+" - bad ICOOR specification." );
		}
		AtomICoor const & conicoor( connection.icoor() );
		// TODO: This has a possibility of an infinite loop if you have connection points which mutually depend on each other
		return ( has_assigned_coords( conicoor.stub_atom(1), assigned, restype ) &&
				has_assigned_coords( conicoor.stub_atom(2), assigned, restype ) &&
				has_assigned_coords( conicoor.stub_atom(3), assigned, restype ) );
	}
	return false; // make the compilier happy.
}

/// @details Contains logic originally from read_topology_file()
void
fill_ideal_xyz_from_icoor(
		core::chemical::ResidueType & restype,
		core::chemical::ResidueGraph const & graph) {
	if( restype.natoms() == 0 ) {
		TR.Warning << "fill_ideal_xyz_from_icoor: residue type has no atoms." << std::endl;
		return;
	}
	if( restype.natoms() == 1 ) {
		restype.set_ideal_xyz( *boost::vertices( graph ).first, Vector(0,0,0) );
		return;
	}
	// As we can't guarantee we'll cover the atoms in any reasonable order,
	// keep going through the atoms, assigning those we're able, and deferring
	// those we aren't. To avoid infinite loops, for each cycle through the atoms
	// make sure we've made progress (assigned at least one xyz).
	std::set< VD > assigned;
	std::queue< VD > atom_queue;
	VIterPair allverts( boost::vertices( graph ) );
	for( VIter iter(allverts.first); iter != allverts.second; ++iter ) {
		atom_queue.push( *iter );
		// For debug:
		restype.set_ideal_xyz( *iter, Vector(99.44, 99.44, 99.44) );
	}
	core::Size natoms( atom_queue.size() ); // The number of atoms left to see this cycle.
	bool progress = false; // A variable to keep track of progress, to avoid infinite loops.

	while( ! atom_queue.empty() ) {
		if( natoms == 0 ) {
			if( ! progress ) {
				TR.Error << "Cannot assign ideal coordinates for atoms ";
				while( ! atom_queue.empty() ) {
					VD a( atom_queue.front() );
					TR.Error << restype.atom_name( a ) << " ";
					atom_queue.pop();
				}
				TR.Error << "in residue type " << restype.name() << std::endl;
				TR.Error << "Do you have a circular ICOOR dependency?" << std::endl;
				utility_exit_with_message( "Unable to assign ideal coordinates for residue type "+restype.name() );
			} else {
				// We've done a round, but have made progress.
				natoms = atom_queue.size();
				progress = false;
			}
		}
		VD child_atom = atom_queue.front();
		atom_queue.pop();
		--natoms;
		AtomICoor const & icoor( restype.icoor(child_atom) );

		ICoorAtomID parent_stub( icoor.stub_atom(1) ), angle_stub( icoor.stub_atom(2) ), torsion_stub( icoor.stub_atom(3) );
		core::Real phi( icoor.phi() ), theta( icoor.theta() ), d( icoor.d() );
	//debug_assert( restype.has(parent_atom) && restype.has(angle_atom) );
		if ( child_atom == parent_stub.vertex() ) { // root atom
			restype.set_ideal_xyz( child_atom, Vector(0,0,0) );
			assigned.insert( child_atom );
			progress = true;
		} else if ( child_atom == angle_stub.vertex() ) { // second atom
			restype.set_ideal_xyz( child_atom, Vector(d,0,0) );
			assigned.insert( child_atom );
			progress = true;
		} else {
			if( ! has_assigned_coords( parent_stub, assigned, restype ) || ! has_assigned_coords( angle_stub, assigned, restype ) ) {
				// Defer until dependencies are met.
				atom_queue.push( child_atom );
				continue;
			}
			Vector torsion_xyz;
			if ( child_atom == torsion_stub.vertex() ) { // third atom
				torsion_xyz = Vector( 1.0, 1.0, 0.0 );
			} else {
				if( ! has_assigned_coords( torsion_stub, assigned, restype ) ) {
					// Defer until dependencies are met.
					atom_queue.push( child_atom );
					continue;
				}
				torsion_xyz = torsion_stub.xyz( restype );
			}
			kinematics::Stub const stub( parent_stub.xyz( restype ),
					angle_stub.xyz( restype ),
					torsion_xyz );
			restype.set_ideal_xyz( child_atom , stub.spherical( phi, theta, d ) );
			assigned.insert( child_atom );
			progress = true;
		}
	} //while( ! atom_queue.empty() )
}


} //chemical
} //core

