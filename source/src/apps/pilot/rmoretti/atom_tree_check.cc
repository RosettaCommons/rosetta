// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/rmoretti/atom_tree_check.cc
/// @brief  Check retreeing with core::kinematics::tree
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/AtomICoor.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/residue_io.hh>

#include <core/id/AtomID.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/tree/BondedAtom.hh>
#include <core/kinematics/tree/JumpAtom.hh>
#include <core/id/types.hh>

#include <core/types.hh>

#include <devel/init.hh>

#include <utility/options/FileVectorOption.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <fstream>


#include <utility/graph/DFS_sort.hh>
#include <numeric/xyzVector.io.hh>

#include <boost/graph/depth_first_search.hpp>
#include <boost/graph/properties.hpp>

using namespace basic::options;

static thread_local basic::Tracer TR( "apps.atomtreeprint" );

void tweak_coords( core::kinematics::tree::AtomOP atom) {
	using namespace core::kinematics::tree;
	using namespace core::id;
	if ( atom->parent() ) {
		atom->set_dof(D,-1);
		atom->set_dof(PHI,-1);
		atom->set_dof(THETA,-1);
		atom->xyz( atom->xyz() + 0.0000001 );
	}

	for( core::Size ii(0); ii < atom->n_children(); ++ii ) {
		tweak_coords(atom->child(ii));
	}

}


void
print_tree( core::kinematics::tree::AtomCOP atom, core::pose::Pose const & pose ) {
	using namespace core::kinematics::tree;
	using namespace core::id;

//	TR << " Atom " << pose.residue(atom->id().rsd()).atom_name(atom->id().atomno());
/*
	AtomCOP parent, gparent, ggparent;
	parent = atom->parent();
	if( parent == 0 ) {
		TR << " None";
		gparent = 0;
	} else {
		TR << " " << pose.residue(parent->id().rsd()).atom_name(parent->id().atomno());
		gparent = parent->parent();
	}
	if( gparent == 0 ) {
		TR << " None";
		ggparent = 0;
	} else {
		TR << " " << pose.residue(gparent->id().rsd()).atom_name(gparent->id().atomno());
		ggparent = gparent->parent();
	}
	if( ggparent == 0 ) {
		TR << " None";
	} else {
		TR << " " << pose.residue(ggparent->id().rsd()).atom_name(ggparent->id().atomno());
	}
	core::kinematics::tree::BondedAtomCOP batom;
	batom = dynamic_cast< core::kinematics::tree::BondedAtom const * >( atom() );
	if( batom == 0 ) {
		TR << " (Not bonded)" << std::endl;
	} else {
		TR << std::fixed;
		TR.precision(3);
		TR << " "  << 57.2957795*batom->dof(PHI) << " " << 57.2957795 * batom->dof(THETA) << " " << batom->dof(D) << std::endl;
	}
*/

	AtomCOP stub, parent;
	parent = atom->parent();
	if( parent ) {
	stub = atom;
	if( stub ) { TR << " " << pose.residue(stub->id().rsd()).atom_name(stub->id().atomno()); } else { TR << " --"; }

		TR << std::fixed;
		TR.precision(3);
	TR << " "  << 57.2957795*atom->dof(PHI) << " " << 57.2957795 * atom->dof(THETA) << " " << atom->dof(D);

	stub = atom->input_stub_atom1();
	if( stub ) { TR << " " << pose.residue(stub->id().rsd()).atom_name(stub->id().atomno()); } else { TR << " --"; }
	stub = atom->input_stub_atom2();
	if( stub ) { TR << " " << pose.residue(stub->id().rsd()).atom_name(stub->id().atomno()); } else { TR << " --"; }
	stub = atom->input_stub_atom3();
	if( stub ) { TR << " " << pose.residue(stub->id().rsd()).atom_name(stub->id().atomno()); } else { TR << " --"; }

	TR << std::endl;
	}

	for( core::Size ii(0); ii < atom->n_children(); ++ii ) {
		print_tree(atom->child(ii), pose);
	}
}

void
print_tree( core::kinematics::tree::AtomCOP atom, core::chemical::ResidueType const & restype ) {
	using namespace core::kinematics::tree;
	using namespace core::id;

	AtomCOP stub, parent;
	parent = atom->parent();
	if( parent ) {
		stub = atom;
		if( stub && stub->id().atomno() ) {
			TR << " " << restype.atom_name(stub->id().atomno()); }
		else if ( stub ) {
			TR << " 00";
		} else {
			TR << " --";
		}

		TR << std::fixed;
		TR.precision(2);
		TR << '\t'  << 57.2957795*atom->dof(PHI) << '\t' << 57.2957795 * atom->dof(THETA) << '\t' << atom->dof(D);

		stub = atom->input_stub_atom1();
		if( stub && stub->id().atomno() ) { TR << " " << restype.atom_name(stub->id().atomno()); } else { TR << " ----"; }
		stub = atom->input_stub_atom2();
		if( stub && stub->id().atomno() ) { TR << " " << restype.atom_name(stub->id().atomno()); } else { TR << " ----"; }
		stub = atom->input_stub_atom3();
		if( stub && stub->id().atomno() ) { TR << " " << restype.atom_name(stub->id().atomno()); } else { TR << " ----"; }
		stub = atom->stub_atom3();
		if( stub && stub->id().atomno() ) { TR << " (" << restype.atom_name(stub->id().atomno())<<")"; } else { TR << " (----)"; }

		TR << '\t'  << atom->x() << '\t' << atom->y() << '\t' << atom->z();

		TR << std::endl;
	}

	for( core::Size ii(0); ii < atom->n_children(); ++ii ) {
		print_tree(atom->child(ii), restype);
	}
}


// This is a test function to make sure that all the critical information about
// the atom tree can be correctly encapsulated by copying over just the atom_id
// and xyz coordinate information.
core::kinematics::tree::AtomOP
rebuild( core::kinematics::tree::AtomCOP atom ) {
	using namespace core::kinematics::tree;
	BondedAtomCOP batom = dynamic_cast< BondedAtom const * >( atom() );
	core::kinematics::tree::AtomOP new_atom;
	if( ! batom ) {
		new_atom = new JumpAtom;
	}
	else {
		new_atom = new BondedAtom;
	}

	new_atom->set_weak_ptr_to_self( new_atom() );

	for( core::Size ii(0); ii < atom->n_children(); ++ii ) {
		new_atom->append_atom( rebuild( atom->child(ii) ) );
	}

	new_atom->id( atom->atom_id() );
	new_atom->xyz( atom->xyz() );

	return new_atom;
}


//The core::kinematics::tree::Atom objects take care of the internal/xyz exchange.
//We just need to chain them together appropriately.
typedef std::map< core::chemical::VD, core::kinematics::tree::AtomOP > VdTreeatomMap;


class RerootRestypeVisitor: boost::default_dfs_visitor {
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
			treeatom_map_[ u ] = new core::kinematics::tree::BondedAtom;
			treeatom_map_[ u ]->set_weak_ptr_to_self( treeatom_map_[ u ].get() );
		}
		treeatom_map_[ u ]->xyz( restype_.atom( u ).ideal_xyz() );
		// For debugging purposes:
		core::id::AtomID id(restype_.atom_index(u), 1);
		treeatom_map_[ u ]->id( id );
		//TR << "Initializing " << restype_.atom_name(u) << " to " << treeatom_map_[ u ]->xyz() << std::endl;
	}
	template <class Vertex, class Graph>
	void start_vertex(Vertex u, const Graph& /*g*/) {
		// The DFS algorithm will restart on disconnected subgraphs --
		// we currently can't handle the disconnected portion of it,
		// so raise an error if that happens.
		// (This means the root needs to be pre-loaded into the treeatom_map
		assert( restype_.has( u ) );
		if( u != root_ ) {
			TR.Error << "ERROR: For ResidueType " << restype_.name() << ", atoms "
					<< restype_.atom_name( root_ ) << " and " << restype_.atom_name( u )
					<< " are not connected via bonds." << std::endl;
			utility_exit_with_message( "Cannot reroot a disconnected ResidueType: ");
		}
	}
	template <class Vertex, class Graph>
	void discover_vertex(Vertex /*u*/, const Graph& /*g*/) {
	}
	template <class Edge, class Graph>
	void examine_edge(Edge /*u*/, const Graph& /*g*/) {
	}
	template <class Edge, class Graph>
	void tree_edge(Edge u, const Graph& g) {
		core::chemical::VD parent( boost::source(u,g) );
		core::chemical::VD child( boost::target(u,g) );
		assert( treeatom_map_.count( parent ) );
		assert( treeatom_map_.count( child ) );
		treeatom_map_[ parent ]->append_atom( treeatom_map_[ child ] );
		//TR << "Adding tree edge from " << restype_.atom_name(parent) << " to " << restype_.atom_name(child) << std::endl;
		//TR << "%%% Current tree" << std::endl;
		//print_tree( treeatom_map_[ root_ ], restype_ );
	}
	template <class Edge, class Graph>
	void back_edge(Edge /*u*/, const Graph& /*g*/) {
	}
	template <class Edge, class Graph>
	void forward_or_cross_edge(Edge /*u*/, const Graph& /*g*/) {
	}
	template <class Vertex, class Graph>
	void finish_vertex(Vertex /*u*/, const Graph& /*g*/) {
	}
};

/// @brief Edge sorting:
/// * Actual bonds before pseudobonds
/// * Rotatable bonds before fixed bonds
/// * Bonds to Heavy atoms before light atoms
/// * Bonds to Concrete atoms before virtual atoms
///
/// This doesn't (need to?) quite match the logic in core/conformation/util.cc:setup_atom_links()
class RerootEdgeSorter {
public:
	RerootEdgeSorter(core::chemical::ResidueGraph graph, core::chemical::ResidueType const & restype):
		graph_(graph),
		restype_(restype)
	{}

	/// Return true if the first argument goes before the second argument
	bool operator()(core::chemical::ED edge1, core::chemical::ED edge2 ) {
		core::chemical::Bond const & bond1( graph_[edge1] );
		core::chemical::Bond const & bond2( graph_[edge2] );

		// * Actual bonds before pseudobonds - TODO better pseudobond testing
		if( bond1.bond_name() == core::chemical::UnknownBond && bond2.bond_name() != core::chemical::UnknownBond) {
			return false;
		} else if ( bond1.bond_name() != core::chemical::UnknownBond && bond2.bond_name() == core::chemical::UnknownBond) {
			return true;
		}
		assert( boost::source(edge1,graph_) == boost::source(edge2,graph_) );
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
///@brief
/// Doing a depth first search here because that's what molfile_to_params.py does:
/// "Protein residues appear to go depth first, so that all chi angles ride on each other."
/// RM: Doing a breadth first search would likely result in a shallower tree, but with possibly
/// different behavior on how ring atom trees are built.
///
/// Assumes:
/// * All bonds and atoms exist in the graph,
/// * The graph is completely connected.
/// * All ideal xyz coordinates are updated.
///
void
reroot_restype( core::chemical::ResidueType & restype, core::chemical::VD root) {
	if( restype.natoms() < 3 ) {
		TR.Warning << "Warning: Cannot re-root residue type with less than three atoms." << std::endl;
		return;
	}
	VdTreeatomMap treeatom_map;
	std::map< core::chemical::VD, boost::default_color_type > colormap;
	boost::associative_property_map< std::map< core::chemical::VD, boost::default_color_type > > color_property_map(colormap);

	// The root of the atom tree needs to be a jump atom,
	// as opposed to the Bonded atom that would be initialized by default in the DFS
	core::kinematics::tree::AtomOP rootatom = new core::kinematics::tree::JumpAtom;
	rootatom->set_weak_ptr_to_self( rootatom() );
	treeatom_map[ root ] = rootatom;
	TR << "Rooting on atom: " << restype.atom_name(root) << std::endl;

	RerootRestypeVisitor visitor( treeatom_map, root, restype );

	//TR << "$$$ before search tree" << std::endl;
	//print_tree( rootatom, restype );

	utility::graph::depth_first_search_sort(restype.graph_, visitor, color_property_map, root, RerootEdgeSorter( restype.graph_, restype ) );

	//TR << "$$$ Pre update tree" << std::endl;
	//print_tree( rootatom, restype );

	rootatom->update_internal_coords( true );

	//TR << "$$$ Post update tree" << std::endl;
	//print_tree( rootatom, restype );

	// Reverse the map for easy lookup
	std::map< core::kinematics::tree::AtomCOP, core::chemical::VD > revmap;
	for( VdTreeatomMap::const_iterator itr(treeatom_map.begin()), itr_end(treeatom_map.end()); itr != itr_end; ++itr) {
		revmap[ itr->second ] = itr->first;
	}

	// Now pull the icoord data out of the kinematic atoms

	core::chemical::VIter iter, iter_end;
	for( boost::tie(iter, iter_end) = boost::vertices(restype.graph_); iter != iter_end; ++iter ) {
		core::chemical::VD const & atomVD( *iter );
		assert( treeatom_map.count( atomVD ) );

		core::kinematics::tree::AtomCOP atom = treeatom_map[ atomVD ];

		core::kinematics::tree::AtomCOP parent, angle, torsion;
		core::Real phi, theta, d;
		parent = atom->parent();
		if( parent ) {
			// Regular, non-root atom
			assert( dynamic_cast< core::kinematics::tree::BondedAtom const * >( atom() ) );

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
			// Define the torsion linearly, if possible, else use an improper through the root atoms' other child.
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
		assert( restype.has( parentVD ) && restype.has( angleVD ) && restype.has( torsionVD ) );

		// Note: set_icoor will automatically setup the atom base for us.
		restype.set_icoor( atomVD, phi, theta, d, parentVD, angleVD, torsionVD );
	}
}

bool
comp(core::Real x, core::Real y ) {
	return x > (y - 1e-6) && x < (y + 1e-6);
}

int
main( int argc, char * argv [] )
{

try {

	  using namespace basic::options;
	  using namespace basic::options::OptionKeys;
	  using namespace core::chemical;

	devel::init(argc, argv);

	utility::options::FileVectorOption & fvec
		= basic::options::option[ basic::options::OptionKeys::in::file::extra_res_fa ];

	core::chemical::ChemicalManager * chem_mang = core::chemical::ChemicalManager::get_instance();
	core::chemical::AtomTypeSetCAP atom_types = chem_mang->atom_type_set("fa_standard");
	core::chemical::ElementSetCAP elements = chem_mang->element_set("fa_standard");
	core::chemical::MMAtomTypeSetCAP mm_atom_types = chem_mang->mm_atom_type_set("fa_standard");
	core::chemical::orbitals::OrbitalTypeSetCAP orbital_types = chem_mang->orbital_type_set("fa_standard");

	core::chemical::ResidueTypeSet res_set;

	TR << "Loading residue types " << std::endl;
	// We don't need to load all the residue types - just the extra_res_fa ones.
	// Grab each and go.

	for(core::Size i = 1, e = fvec.size(); i <= e; ++i) {
		utility::file::FileName fname = fvec[i];
		std::string filename = fname.name();

		core::chemical::ResidueTypeCOP restype( read_topology_file(
				filename, atom_types, elements, mm_atom_types, orbital_types, &res_set ) );

		TR << "############## Tree for " << fname.base()<< std::endl;

		core::pose::Pose pose;
		core::conformation::Residue res(*restype, true );
		pose.append_residue_by_jump( res, 0 ); // Anchor doesn't matter, as it is the first residue

		core::kinematics::AtomTree const & atree( pose.atom_tree() );

		core::kinematics::tree::AtomOP root( const_cast< core::kinematics::tree::Atom* >( atree.root()() ) );

		//print_tree( root, pose );

		//TR << "############## Rebuilding " << std::endl;

		core::kinematics::tree::AtomOP rebuilt( rebuild( root) );
		//root->update_xyz_coords();
		//tweak_coords( root );
		rebuilt->update_internal_coords( true );

		//print_tree( rebuilt, pose );

		//TR << "############# Retreeing " << std::endl;

		// Assume that the root is the neighbor atom -- typically the case for ligands
		core::chemical::ResidueTypeOP restype2( new core::chemical::ResidueType(*restype) );
		core::chemical::VD root2( restype2->vd_from_index( restype2->nbr_atom() ) );

		//TR << "### Before tree " << std::endl;
		//print_tree( root, pose );

		reroot_restype( *restype2, root2 );

		core::pose::Pose pose2;
		core::conformation::Residue res2(*restype2, true );
		pose2.append_residue_by_jump( res2, 0 ); // Anchor doesn't matter, as it is the first residue

		core::kinematics::AtomTree const & atree2( pose2.atom_tree() );
		core::kinematics::tree::AtomOP root2a( const_cast< core::kinematics::tree::Atom* >( atree2.root()() ) );

		//print_tree( root2a, pose2 );

		bool mismatch(false);

//		for( core::Size ii(1); ii <= restype->natoms(); ++ii ) {
//			core::chemical::AtomICoor const icoor( restype->icoor(ii) );
//			core::chemical::AtomICoor const icoor2( restype2->icoor( restype2->atom_index(restype->atom_name(ii) ) ) );
//			if( restype->atom_name( icoor.stub_atom(1).vertex() ) == restype2->atom_name( icoor2.stub_atom(1).vertex() ) &&
//					restype->atom_name( icoor.stub_atom(2).vertex() ) == restype2->atom_name( icoor2.stub_atom(2).vertex() ) &&
//					restype->atom_name( icoor.stub_atom(3).vertex() ) == restype2->atom_name( icoor2.stub_atom(3).vertex() )  ) {
//				if( ! comp( icoor.d(), icoor2.d() ) ) {
//					TR.Error << "Mismatched distances for " << restype->atom_name(ii) << std::endl;
//				}
//				if( ! comp(icoor.theta(), icoor2.theta() ) ) {
//					TR.Error << "Mismatched thetas for " << restype->atom_name(ii) << std::endl;
//				}
//				if( ! comp(icoor.phi(), icoor2.phi() ) ) {
//					TR.Error << "Mismatched phis for " << restype->atom_name(ii) << std::endl;
//				}
//			} else {
//				TR.Error << "Mismatched atoms for " << restype->atom_name(ii) << std::endl;
//				TR.Error << "'" << restype->atom_name( icoor.stub_atom(1).vertex() ) << "' '" << restype2->atom_name( icoor2.stub_atom(1).vertex() ) << "' '" <<
//						restype->atom_name( icoor.stub_atom(2).vertex() ) << "' '" << restype2->atom_name( icoor2.stub_atom(2).vertex() ) << "' '" <<
//						restype->atom_name( icoor.stub_atom(3).vertex() ) << "' '" << restype2->atom_name( icoor2.stub_atom(3).vertex() )  << "'" << std::endl;
//
//				mismatch=true;
//			}
//		}

		// Better test - build ideal coordinates from new icoors tables, check for differences.

		restype2->fill_ideal_xyz_from_icoor();

		//Now make sure we're in the same reference frame.
		//Code based on Residue::orient_onto_residue()
		core::Size center, nbr1, nbr2;
		restype2->select_orient_atoms( center, nbr1, nbr2 );

		core::Size  src_center = restype->atom_index( restype->atom_name( center ));
		core::Size  src_nbr1 = restype->atom_index( restype->atom_name( nbr1 ));
		core::Size  src_nbr2 = restype->atom_index( restype->atom_name( nbr2 ));
		using core::kinematics::Stub;
		core::Vector const
			rot_midpoint ( 0.5 * ( restype2->atom(     nbr1 ).ideal_xyz() + restype2->atom(     nbr2 ).ideal_xyz() ) ),
			src_midpoint ( 0.5 * (  restype->atom( src_nbr1 ).ideal_xyz() +  restype->atom( src_nbr2 ).ideal_xyz() ) );

		Stub rot_stub( restype2->atom( center ).ideal_xyz(),
									 rot_midpoint,
									 restype2->atom( nbr1 ).ideal_xyz() );

		Stub src_stub( restype->atom( src_center ).ideal_xyz(),
									 src_midpoint,
									 restype->atom( src_nbr1 ).ideal_xyz() );

		// this could be made faster by getting the composite rotation and translation

		assert( restype->natoms() == restype2->natoms() );
		for ( core::Size i=1; i<= restype2->natoms(); ++i ) {
			std::string atom_name(restype2->atom_name(i));
			assert( restype->has( atom_name ) );
			core::Vector const old_xyz( restype2->atom(i).ideal_xyz() );
			core::Vector const new_xyz( src_stub.local2global( rot_stub.global2local( old_xyz ) ) );
			core::Vector const original_xyz( restype->atom(atom_name).ideal_xyz() );
			if( original_xyz.distance_squared( new_xyz ) > 0.00000001 ) {
				TR << "Atom " << atom_name << " inappropriately placed " << new_xyz << std::endl;
				TR << "     " << "    "    << "              should be " << original_xyz << std::endl;
				TR << "     " << "    "    << "  (without translation) " << old_xyz << std::endl;
			mismatch=true;
			}
			//restype2->atom(i).ideal_xyz( new_xyz );
		}


//		assert( restype->natoms() == restype2->natoms() );
//		for( core::Size ii(1); ii <= restype->natoms(); ++ii ) {
//			std::string atom_name(restype->atom_name(ii));
//			assert( restype2->has( atom_name) );
//			core::Vector p1(restype->atom(atom_name).ideal_xyz()), p2(restype2->atom(atom_name).ideal_xyz());
//			if( p1.distance_squared(p2 ) > 0.00000001 ) {
//				TR << "Atom " << atom_name << " inappropriately placed " << p2 << " should be " << p1 << std::endl;
//				mismatch=true;
//			}
//		}

		if( mismatch ) {
			pose.dump_pdb("original_" + fname.base() + ".pdb");
			pose2.dump_pdb( "modified_" + fname.base() + ".pdb");
			TR << "### Before tree " << std::endl;
			print_tree( root, *restype );
			TR << "### After tree " << std::endl;
			print_tree( root2a, *restype2 );
			TR << "### For Restype Before " << std::endl;
			TR << "\n";
			// Here we assume that the neighbor atom is the ICOOR root.
			core::chemical::pretty_print_atomicoor(TR, restype->icoor( restype->nbr_atom() ), *restype);
			TR << std::endl;
			TR << "### For Restype After" << std::endl;
			TR << "\n";
			core::chemical::pretty_print_atomicoor(TR, restype2->icoor( restype2->atom_index(root2) ), *restype2);
			TR << std::endl;

		}
		TR << "############## End tree for " << fname.base()<< std::endl;
	}


	TR << "Done trees" << std::endl;

} catch ( utility::excn::EXCN_Base const & e ) {
	std::cout << "caught exception " << e.msg() << std::endl;
}

}

