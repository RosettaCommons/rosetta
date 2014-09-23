// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/util.cc
/// @brief  Kinematics utility functions
/// @author Phil Bradley

// Unit headers
#include <core/kinematics/util.hh>
//#include <core/kinematics/tree/IntraResidueStubJumpAtom.hh>
// Package headers
// AUTO-REMOVED #include <core/id/DOF_ID_Mask.hh>
// AUTO-REMOVED #include <core/kinematics/constants.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/JumpAtom.hh>
#include <core/kinematics/tree/BondedAtom.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/types.hh>

// Project headers
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED #include <core/chemical/ResidueConnection.hh>
// AUTO-REMOVED #include <core/chemical/ResConnID.hh>

// AUTO-REMOVED #include <core/id/DOF_ID_Map.hh>
#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>

// Numeric headers
#include <numeric/random/random.hh>

// C++ headers
#include <cassert>

#include <utility/vector1.hh>

#include <ObjexxFCL/format.hh>


namespace core {
namespace kinematics {

using namespace id;
static thread_local basic::Tracer TR( "core.kinematics.util" );


///////////////////////////////////////////////////////////////////////////////
// wrapper for the add_atom recursive function which doesn't know anything
// about Residues
//
/// @details recursively called until all atoms in this reside are added.
/// @note The Atoms that are allocated within this function are owned by the atom_ptr
/// which is a container of AtomOP's (see AtomPointer.fwd.hh)
///
/// @note  Returns a raw pointer for internal tree links but memory management is handled by atom_ptr
///


tree::AtomOP
add_atom(
	int const atomno,
	int const seqpos,
	Links const & links,
	AtomPointer1D & atom_ptr, // owns all the newly allocated atoms
	bool const add_jump_atom
)
{
	using namespace tree;

	// create new atom
	AtomOP const atom_p( add_jump_atom ? static_cast< Atom* >( new JumpAtom()) : static_cast< Atom* >(new BondedAtom()));

	// fill in the atom_ptr data
	assert( atom_ptr[ atomno ] == 0 );
	atom_ptr[ atomno ] = atom_p;

	// set the atom_id information
	atom_p->id( AtomID( atomno, seqpos ) );
	atom_p->parent( 0 );

	utility::vector1< Size > const & nbrs( links[ atomno ] );
	for ( Size i=1, i_end = nbrs.size(); i<= i_end; ++i ) {
		int const nbr( nbrs[i] );
		if ( atom_ptr[ nbr ] != 0 ) continue;
		atom_p->append_atom( add_atom( nbr, seqpos, links, atom_ptr, false ));
	}

	return atom_p();
}

///////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

/// \brief pick a postion in n_res as the cutpoint
///
///  this is done based on probability info stored in cut_bias_sum.
/// This function is used during fold_tree construction.
int
pick_loopy_cutpoint(
	Size const n_res,
	ObjexxFCL::FArray1D_float const & cut_bias_sum
)
{

	float r = numeric::random::rg().uniform() * cut_bias_sum( n_res );

	int cutpoint( 0 );

	for ( Size i = 1; i <= n_res; ++i ) {
		if ( r > cut_bias_sum(i-1) && r <= cut_bias_sum(i) ) {
			cutpoint = i;
		}
	}

	if ( cutpoint == 0 ) {
		TR.Warning << "pick loopy cutpoint = 0! setting = 1" << std::endl;
		cutpoint = 1;
	}

	return cutpoint;
}



///////////////////////////////////////////////////////////////////////////////
/// helper function for setup_backrub_atom_tree
tree::AtomOP
setup_cloned_atom(
	tree::AtomCOP, // old_atom,
	utility::vector1< id::AtomID > const & // exclude
)
{
	utility_exit_with_message("needs to be refactored to meet new tree-building guidelines");
	tree::AtomOP new_atom( new tree::BondedAtom() );
	/*
	new_atom->parent(0);
	new_atom->id ( old_atom->id () );
	new_atom->xyz( old_atom->xyz() );

	// add old_atom's (non-exclude-) children
	for ( Size i=0; i< Size(old_atom->n_children()); ++i ) { // 0-indexing AAARRGGHHH!!!
		AtomID const & child_id( old_atom->child(i)->id() );
		if ( std::find( exclude.begin(), exclude.end(), child_id ) == exclude.end() ) {
			new_atom->append_atom( old_atom->child(i)->clone(0) ); // recursively copies child's children
		}
	}
	*/
	return new_atom;
}


///////////////////////////////////////////////////////////////////////////////
/// in principal, this should be the inverse of the build_tree function at the beginning
///
// void
// fold_tree_from_atom_tree(
// 	conformation::ResidueCAPs const & residues
// 	AtomTree const & atom_tree,
// 	FoldTree & fold_tree
// )
// {

// 	/// get all inter-residue connections
// 	utility::vector1< Edge

tree::AtomOP
setup_backrub_atom_tree(
	utility::vector1< AtomID >, // mainchain, // make our own local copy
	AtomID const &, // downstream_id, // mainchain child of last mainchain atom
	AtomPointer2D const &, // old_atom_pointer,
	utility::vector1< std::pair< Size, Size > > const &, // edges,
	Size const //first_new_pseudo_residue
)
{

	// In the interest of making the gameguys' deadline for refactoring atomtree I'm not thinking about this right now.
	utility_exit_with_message("This code needs to be refactored to use AtomOPs" );
	/*
	using utility::vector1;

	Vector const pseudo_offset( 0,0,0 );
// 	Vector const pseudo_offset( 0.25, 0.25, 0.25 );

	// constraints on the edge list (note that edges are ordered, ie i,j means edge *from* i *to* j (i-->j)
	//
	// -- 1,nbb is the first edge
	//
	// -- for each (i,j) after the first there exists an earlier edge (x,i) (ie, have to have seen first index already)
	//
	// -- no crossing: no pairs (i,j) and (k,l) with i<k<j<l
	//
	// -- for each edge (a,b), a must be downstream of all edges that contain (a,b)
	//    ie if (i,j) contains (a,b) (meaning that the interval [min(i,j),max(i,j)] contains the corresponding ab interval
	//    then there must be a path of edges leading from j to a (j,p1) (p1,p2) (p2,p3) .... (pN,a)
	//    so eg we cant have 1,10 10,7 7,4 1,4 (violated by (1,4) edge: no path from 10 --> 1 and [1,10] contains [1,4])
	//
	// -- no edges (i,j) with |i-j|<=1
	//

	Size const nbb( mainchain.size() );
	mainchain.push_back( downstream_id ); // so that we dont clone this guy
	assert( edges[1].first == 1 && edges[1].second == nbb );

	/// we're going to add pseudo residues, one for each edge in the edges vector
	Size n(0);

	// book-keeping

	// map from each mainchain_atom to the (possibly empty) list of pseudo-residues
	vector1< vector1< Size > > mainchain_to_pseudo( nbb );

	vector1< bool > seen( nbb, false );

	vector1< Atom* > mainchain_atom_pointer( nbb, 0 );
	vector1< Atom* > pseudo_atom_pointer( edges.size(), 0 );

	// create the root of the tree
	AtomOP root( setup_cloned_atom( old_atom_pointer[ mainchain[1] ], mainchain ) );
	mainchain_atom_pointer[ 1 ] = root;
	seen[ 1 ] = true;

	while ( n < edges.size() ) {
		++n;
		std::pair< Size, Size > const & edge( edges[n] );

		Size const a( edge.first );
		Size const b( edge.second );
		assert( seen[a] );

		// what should our anchor atom be?
		// look at all pseudo rsds associated to a
		AtomOP anchor_atom;
		Size anchor(0);
		vector1< Size > const & pseudo_a( mainchain_to_pseudo[a] );
		for ( Size i=1; i<= pseudo_a.size(); ++i ) {
			Size const nn( pseudo_a[i] );
			assert( edges[nn].second == a );
			if ( ( edges[nn].first < b && b < a ) || ( edges[nn].first > b && b > a ) ) {
				anchor = nn; // dont break -- want the smallest segment containing a,b
			}
		}
		if ( anchor ) {
			anchor_atom = pseudo_atom_pointer[ anchor ];
		} else {
			// connect to the authentic a
			anchor_atom = mainchain_atom_pointer[ a ];
		}

		// has b been seen before?
		AtomCOP old_b_atom( old_atom_pointer[ mainchain[b] ] );
		if ( !seen[ b ] ) {
			seen[b] = true;
			// add b and b's non-mainchain children
			Atom* b_atom( setup_cloned_atom( old_b_atom, mainchain ) );
			anchor_atom->insert_atom( b_atom ); // NOTE: insert atom
			mainchain_atom_pointer[ b ] = b_atom;
		}


		// add a new pseudo rsd at b's position
		mainchain_to_pseudo[ b ].push_back( n );
		Size const pseudo_b_seqpos( first_new_pseudo_residue + n - 1 );
		AtomID const pseudo_b_id( AtomID( 1, pseudo_b_seqpos ) );
		// delete the old one
		Atom* old_pseudo_b_atom( old_atom_pointer[ pseudo_b_id ] );
		old_pseudo_b_atom->parent()->delete_atom( old_pseudo_b_atom );
		delete old_pseudo_b_atom;
		// create a new one
		Atom* pseudo_b_atom( new BondedAtom() );
		pseudo_b_atom->id( pseudo_b_id );
		pseudo_b_atom->xyz( old_b_atom->xyz() + pseudo_offset );
		anchor_atom->append_atom( pseudo_b_atom );
		pseudo_atom_pointer[ n ] = pseudo_b_atom;
		std::cout << "new pseudo! " << n << ' ' << a << ' ' << b << ' ' << mainchain[a] << ' ' << mainchain[b] << ' ' <<
			pseudo_b_seqpos << ' ' << old_b_atom->id() << ' ' <<
			old_b_atom->xyz()[0] << ' ' <<
			old_b_atom->xyz()[1] << ' ' <<
			old_b_atom->xyz()[2] << std::endl;

		// check if this edge is terminal, if so add all intervening mainchain atoms and their children
		{
			bool terminal( true );
			for ( Size n2=n+1; n2<= edges.size(); ++n2 ) {
				if ( ( edges[n2].first == b ) &&
						 ( a < b && edges[n2].second < b || a > b && edges[n2].second > b ) ) {
					terminal = false;
				}
			}

			if ( terminal ) {
				int const dir( b<a ? 1 : -1 );
				AtomOP parent_atom( pseudo_b_atom );
				for ( int c=b+dir; c != (int)a; c += dir ) {
					assert( !seen[c] );

					// add c and c's non-mainchain children
					AtomCOP old_c_atom( old_atom_pointer[ mainchain[c] ] );
					Atom* c_atom( setup_cloned_atom( old_c_atom, mainchain ) );
					parent_atom->insert_atom( c_atom ); // at front of list since this is mainchain. may have already added kids
					mainchain_atom_pointer[ c ] = c_atom;
					seen[c] = true;

					parent_atom = c_atom;

				} // walk from b->a adding the mainchain atoms and their children to the tree

			} // if terminal

		} // scope


	} // loop over edges, add one pseudo rsd for each edge

	// confirm that all mainchain atoms have been seen
	for ( Size i=1; i<= nbb; ++i ) {
		assert( seen[i] );
	}
	return root;
	*/
	return 0;

}

/// @brief prints something like this ***1***C***1*********2***C********3****C****2********3*****
void
simple_visualize_fold_tree( FoldTree const & fold_tree, std::ostream& out ) {
	for ( Size pos = 1; pos <= fold_tree.nres(); pos++ ) {
		bool special( false );
		if ( fold_tree.is_jump_point( pos ) ) {
			for ( Size jnr=1; jnr<=fold_tree.num_jump(); jnr++ ) {
				if ( (Size) fold_tree.jump_edge( jnr ).start() == pos || (Size)fold_tree.jump_edge( jnr ).stop() == pos ) {
					if ( special ) out << "/";
					out << jnr;
					special = true;
				}
			}
		}
		if ( fold_tree.is_cutpoint( pos ) ) {
			if ( special ) out << "/";
			out << "C";
			special = true;
		}
		if (!special ) out << "*";
	}
	out << std::endl;
}


/// @brief prints something like this ***1***C***1*********2***C********3****C****2********3*****
///                                   **********xxxxxxxxxxxxx************************************
void
simple_visualize_fold_tree_and_movemap( FoldTree const & fold_tree, MoveMap const& mm, std::ostream& out ) {
	std::string move;
	out << "\n";
	for ( Size pos = 1; pos <= fold_tree.nres(); pos++ ) {
		bool special( false );
		if ( fold_tree.is_jump_point( pos ) ) {
			for ( Size jnr=1; jnr<=fold_tree.num_jump(); jnr++ ) {
				if ( (Size) fold_tree.jump_edge( jnr ).start() == pos || (Size)fold_tree.jump_edge( jnr ).stop() == pos ) {
					if ( special ) {
						out << "/";
						move.push_back( '.' );
					}
					out << jnr;
					move.push_back( mm.get_bb( pos ) ? '*' : 'x' );
					special = true;
				}
			}
		}
		if ( fold_tree.is_cutpoint( pos ) ) {
			if ( special ) {
				out << "/";
				move.push_back( '.' );
			}
			out << "C";
			move.push_back( mm.get_bb( pos ) ? '*' : 'x' );
			special = true;
		}
		if (!special ) {
			out << "*";
			move.push_back( mm.get_bb( pos ) ? '*' : 'x' );
		}
	}
	out << "\n" << move;
	out << std::endl;
}

/// @brief prints something like this ***1***C***1*********2***C********3****C****2********3*****
///                                   **********xxxxxxxxxxxxx************************************
void
simple_visualize_fold_tree_and_movemap_bb_chi( FoldTree const & fold_tree, MoveMap const& mm, std::ostream& out ) {
	std::string move;
	std::string move_chi;
	out << "\n";
	for ( Size pos = 1; pos <= fold_tree.nres(); pos++ ) {
		bool special( false );
		if ( fold_tree.is_jump_point( pos ) ) {
			for ( Size jnr=1; jnr<=fold_tree.num_jump(); jnr++ ) {
				if ( (Size) fold_tree.jump_edge( jnr ).start() == pos || (Size)fold_tree.jump_edge( jnr ).stop() == pos ) {
					if ( special ) {
						out << "/";
						move.push_back( '.' );
					}
					out << jnr;
					move.push_back( mm.get_bb( pos ) ? '*' : 'x' );
					move_chi.push_back( mm.get_chi( pos ) ? '*' : 'x' );
					special = true;
				}
			}
		}
		if ( fold_tree.is_cutpoint( pos ) ) {
			if ( special ) {
				out << "/";
				move.push_back( '.' );
				move.push_back( '.' );
			}
			out << "C";
			move.push_back( mm.get_bb( pos ) ? '*' : 'x' );
			move_chi.push_back( mm.get_chi( pos ) ? '*' : 'x' );
			special = true;
		}
		if (!special ) {
			out << "*";
			move_chi.push_back( mm.get_chi( pos ) ? '*' : 'x' );
			move.push_back( mm.get_bb( pos ) ? '*' : 'x' );
		}
	}
	out << "\n" << move << "\n" << move_chi;
	out << std::endl;
}


///@brief linearizes (or defoliates, if you prefer) a FoldTree.  "default" FoldTrees produced by the PDB reader have all chains (peptide edges) starting from jumps relative to residue 1.  This code modifies the tree to instead have all the jumps be relative to the preceding edge.  It is not tested with ligands and will not work with "functional" jumps.  From A to B:
///A:FOLD_TREE  EDGE 1 78 -1  EDGE 1 79 1   EDGE 79 454 -1  EDGE 1 455 2    EDGE 455 540 -1  EDGE 1 541 3    EDGE 541 697 -1
///B:FOLD_TREE  EDGE 1 78 -1  EDGE 78 79 1  EDGE 79 454 -1  EDGE 454 455 2  EDGE 455 540 -1  EDGE 540 541 3  EDGE 541 697 -1
core::kinematics::FoldTree
linearize_fold_tree( core::kinematics::FoldTree const & tree ) {
	core::kinematics::FoldTree newtree;
	for( core::kinematics::FoldTree::const_iterator it(tree.begin()), end(tree.end()); it != end; ++it){
		//if it is not a jump, we don't modify it
		if( !it->is_jump() ) newtree.add_edge(*it);
		//if it is a jump, we move start() to stop-1.  This is naive but works for the intended case.
		else newtree.add_edge(core::kinematics::Edge(it->stop()-1, it->stop(), it->label()));
	}
	return newtree;
}







////////////////////////// sheffler visualize FT ////////////////////////////////////

void replace_substr(std::string& str, const std::string from, const std::string to){
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
    }
}

std::string operator*(std::string const s, size_t n) {
    std::string r; // empty string
    r.reserve(n * s.size());
    for (size_t i=0; i<n; i++) r += s;
    return r;
}

std::string pad_dash(Size npad, std::string s) {
	// return s;
	int topad = npad-s.size();
	if( topad > 0 && topad%2==0) return std::string("-")*(topad/2  ) + s + std::string("-")*(topad/2);
	if( topad > 0 && topad%2==1) return std::string("-")*(topad/2+1) + s + std::string("-")*(topad/2);
	else return s;
}
std::string pad_dash_left(Size npad, std::string s) {
	// return s;
	int topad = npad-s.size();
	if( topad > 0 ) return s + std::string("-")*topad;
	else return s;
}
std::string pad_dash_right(Size npad, std::string s) {
	// return s;
	int topad = npad-s.size();
	if( topad > 0 ) return std::string("-")*topad + s;
	else return s;
}

struct Node {
	Node(std::string _name, Size _jnum, Size _jumpfrom, Size _jumpto, char _jumpmark=(char)NULL, Size _follows=0)
		: name(_name), jnum(_jnum), jumpfrom(_jumpfrom), jumpto(_jumpto), prefix_len(8), follows(_follows), jumpmark(_jumpmark), parent(NULL) {}
	~Node(){ for(utility::vector1<Node*>::iterator i = children.begin(); i != children.end(); ++i) delete *i; }
	void setparent(Node *p) {
		parent = p;
		parent->children.push_back(this);
	}
	Node* root() { return parent==NULL ? this : parent->root();	}
	std::string str() {
		using ObjexxFCL::string_of;
		using std::string;
		string s;
		for(utility::vector1<Node*>::const_iterator ic = children.begin(); ic != children.end(); ++ic) {
			string mark = ((*ic)->jumpmark==(char)NULL) ? "-" : string("")+(*ic)->jumpmark;
			string vchar = children.size()==1 ? "" : (*ic==children.back() ? " "  : "|");
			string schar = children.size()==1 ? "" : (*ic==children.back() ? "\\" : "|");
			std::string uprsd = ((*ic)->jumpfrom!=0) ? pad_dash_left(3,string_of((*ic)->jumpfrom)) : std::string("");
			std::string dnrsd = ((*ic)->jumpto!=0) ? pad_dash_right(4,">"+string_of((*ic)->jumpto))+":" : std::string(">");
			std::string jstr;
			if( (*ic)->follows!=0 ) { // no mark if follows
				jstr = "j" + string_of((*ic)->jnum) + ( ((*ic)->follows) ? "="+string_of((*ic)->follows) : "" );
			} else {
				jstr = mark + "j" + string_of((*ic)->jnum) + mark;
			}
			std::string prefix = schar + uprsd + "--"+pad_dash(prefix_len,jstr) + "--" + dnrsd;
			string news = (*ic)->str();
 			string pad = string(" ")*(prefix.size()-1);
 			if(children.size()==1) pad += string(" ")*(name.size()+1);
			replace_substr(news,"\n","\n"+vchar+pad);
			if(children.size()==1) s +=      prefix+news;
			else                   s += "\n"+prefix+news;
		}
		s = name + s;
		return s;
	}
	std::string name;
	Size jnum, jumpfrom, jumpto, prefix_len, follows;
	char jumpmark;
	Node *parent; // worried about cycles, no owning_ptr
	utility::vector1<Node*> children;
};

struct TreeVizBuilder {
	utility::vector1<Size> lb_,ub_;
	FoldTree const & ft;

	TreeVizBuilder(core::kinematics::FoldTree const & _ft) : ft(_ft) {
		lb_.resize(ft.nres(),0);
		ub_.resize(ft.nres(),0);
	}

	void get_ft_node_bounds(Size res, Size & out_lb, Size & out_ub) {
		if(lb_[res]==0)	{
			lb_[res]=1, ub_[res]=ft.nres();
			for(int i = 1; i <= ft.num_cutpoint(); ++i) {
				Size c = (Size)ft.cutpoint(i);
				if( c <  res ) lb_[res] = std::max(lb_[res],c+1);
				if( c >= res ) ub_[res] = std::min(ub_[res],c  );
			}
		}
		out_lb = lb_[res];
		out_ub = ub_[res];
	}

	Size get_ft_node_lower_bound(Size res) {
		Size lb,ub;	get_ft_node_bounds(res,lb,ub);
		return lb;
	}

	Size is_single(Size res) {
		Size lb,ub;	get_ft_node_bounds(res,lb,ub);
		return lb==ub;
	}


	Size get_ft_node_subroot(Size res) {
		Size lb,ub;	get_ft_node_bounds(res,lb,ub);
		for(Size i = 1; i <= ft.num_jump(); ++i){
			Size dn = ft.downstream_jump_residue(i);
			if( lb <= dn && dn <= ub ) return dn;
		}
		return lb; // default subroot is first res
	}

	void expand_node_labels_partial_by_contig(std::map<Size,std::string> & node_labels_partial) {
		utility::vector1<Size> tocheck;
		for(std::map<Size,std::string>::iterator i = node_labels_partial.begin(); i != node_labels_partial.end(); ++i) {
			tocheck.push_back(i->first);
		}
		for(utility::vector1<Size>::const_iterator i = tocheck.begin(); i != tocheck.end(); ++i) {
		for(utility::vector1<Size>::const_iterator j = tocheck.begin(); j != tocheck.end(); ++j) {
			if( *i == *j ) continue;
			// TR << "check " << *i << " " << *j << endl;
			if( get_ft_node_lower_bound(*i) == get_ft_node_lower_bound(*j) ) {
				if(node_labels_partial[*i] != node_labels_partial[*j]) {
					utility_exit_with_message("non-matching node_labels_partial requested in same FT contig");
				}
			}
			}
			Size lb,ub;	get_ft_node_bounds(*i,lb,ub);
			for(Size ir = lb; ir <= ub; ++ir) {
				node_labels_partial[ir] = node_labels_partial[*i];
			}
		}
	}

	utility::vector1<std::string>
	get_res_nodenames( std::map<Size,std::string> node_labels_partial ){
		using ObjexxFCL::format::I;
		int npad = 1;
		if( ft.nres() >    9 ) npad = 2;
		if( ft.nres() >   99 ) npad = 3;
		if( ft.nres() >  999 ) npad = 4;
		if( ft.nres() > 9999 ) npad = 5;
		npad = 0;
		utility::vector1<std::string> names(ft.nres(),"ERROR_this_name_was_not_set");
		for(Size i = 1; i <= ft.nres(); ++i) {
			Size lb,ub;	get_ft_node_bounds(i,lb,ub);
			std::string lbl = node_labels_partial.count(lb) ? node_labels_partial[lb] : "Contig";
			// Size rt = get_ft_node_subroot(i);
			// std::string resrange = ((lb==rt)?"":I(npad,lb)+"<-") + I(npad,rt) + ((ub==rt)?"":"->"+I(npad,ub));
			std::string resrange = I(npad,lb) + ((ub==lb)?"":"-"+I(npad,ub));
			names[i] = lbl + "(" + resrange + ")";
		}
		return names;
	}

	Size get_jump_num_to_contig_of_resi( Size resi) {
		Size lb,ub;	get_ft_node_bounds(resi,lb,ub);
		for(Size i=1; i <= ft.num_jump(); ++i) {
			if( lb <= (Size)ft.downstream_jump_residue(i) && (Size)ft.downstream_jump_residue(i) <= ub ) return i;
		}
		return 0;
	}

};

std::string
visualize_fold_tree(
	FoldTree const & ft,
	std::map<Size,std::string> const & node_labels_partial_in,
	std::map<Size,char> const & mark_jump_to_res,
	std::map<Size,Size> const & jump_follows
){
	TreeVizBuilder tvb(ft);

	std::map<Size,std::string> node_labels_partial(node_labels_partial_in);
	tvb.expand_node_labels_partial_by_contig(node_labels_partial);
	utility::vector1<std::string> res_to_nodename = tvb.get_res_nodenames(node_labels_partial);

	// make Node array
	std::map<std::string,Node*> nodemap;
	for(Size ir = 1; ir <= ft.nres(); ++ir) {
		if( nodemap.count(res_to_nodename[ir]) ) continue;
		Size lb,ub;	tvb.get_ft_node_bounds(ir,lb,ub);
		Size jnum = tvb.get_jump_num_to_contig_of_resi(ir);
		Size jumpfrom = (jnum&&!tvb.is_single(ft.  upstream_jump_residue(jnum))) ? ft.  upstream_jump_residue(jnum) : 0;
		Size jumpto   = (jnum&&!tvb.is_single(ft.downstream_jump_residue(jnum))) ? ft.downstream_jump_residue(jnum) : 0;
		char markjump = (char)NULL;
		for(Size jr = lb; jr <= ub; ++jr) if(mark_jump_to_res.count(jr)) markjump = mark_jump_to_res.find(jr)->second;
		Size follows = jump_follows.find(jnum)!=jump_follows.end() ? jump_follows.find(jnum)->second : (Size)0;
		nodemap[res_to_nodename[ir]] = new Node(res_to_nodename[ir],jnum,jumpfrom,jumpto,markjump,follows);
	}
	// set tree topology
	for(Size i = 1; i <= ft.num_jump(); ++i) {
		std::string up = res_to_nodename[ft.  upstream_jump_residue(i)];
		std::string dn = res_to_nodename[ft.downstream_jump_residue(i)];
		if(up==dn) {
			std::cerr << "bad jump " << i << " from " << ft.upstream_jump_residue(i) << " to " << ft.downstream_jump_residue(i) << std::endl;
			utility_exit_with_message("BAD FOLD TREE NODES!!!!!!!");
		}
		nodemap[dn]->setparent(nodemap[up]);
	}
	// sanity check: make sure tree is connected
	for(std::map<std::string,Node*>::const_iterator i = nodemap.begin(); i != nodemap.end(); ++i) {
		std::string rootname0 = nodemap.begin()->second->root()->name;
		if( rootname0 != i->second->root()->name ) utility_exit_with_message("Nodes not connected!!!");
		// std::cerr << "=========================== " << i->second->name << " ===========================" << std::endl;
		// std::cerr << i->second->str() << std::endl;
		// std::cerr << "========================================================================" << std::endl;
	}

	return nodemap.begin()->second->root()->str();
}
std::string
visualize_fold_tree( FoldTree const & ft ) {
	std::map<Size,std::string> empty_labels;
	std::map<Size,char> empty_marks;
	std::map<Size,Size> empty_follows;
	return visualize_fold_tree( ft, empty_labels, empty_marks, empty_follows );
}

std::string
visualize_fold_tree(	FoldTree const & ft, std::map<Size,std::string> const & node_labels_partial ) {
	std::map<Size,char> empty_marks;
	std::map<Size,Size> empty_follows;
	return visualize_fold_tree( ft, node_labels_partial, empty_marks, empty_follows );
}

std::string
visualize_fold_tree(	FoldTree const & ft, std::map<Size,char> const & mark_jump_to_res ) {
	std::map<Size,std::string> empty_labels;
	std::map<Size,Size> empty_follows;
	return visualize_fold_tree( ft, empty_labels, mark_jump_to_res, empty_follows );
}

// std::string show_foldtree(core::conformation::symmetry::SymmetricConformation const & symm_conf, SymmData const & symmdata) {
// 	Size Nreal = symm_conf.Symmetry_Info()->num_total_residues_without_pseudo();
// 	// get optional labels
// 	std::map<Size,std::string> node_labels_partial;
// 	for(std::map<Size,std::string>::const_iterator i = symmdata.get_virtual_num_to_id().begin(); i != symmdata.get_virtual_num_to_id().end(); ++i) {
// 		node_labels_partial[i->first+Nreal] = i->second;
// 	}
// 	// get non-unique names for each res, inefficient but enures coverage
// 	utility::vector1<std::string> res_to_nodename = get_res_nodenames(symm_conf.fold_tree(),node_labels_partial);
// 	// make Node array
// 	std::map<std::string,Node*> nodemap;
// 	for(Size ir = 1; ir <= symm_conf.size(); ++ir) {
// 		std::string dofstr = "";
// 		std::map<Size,SymDof> dofs( symm_conf.Symmetry_Info()->get_dofs() );
// 		for(std::map<Size,SymDof>::const_iterator i = dofs.begin(); i != dofs.end(); ++i) {
// 			if(symm_conf.fold_tree().downstream_jump_residue(i->first)==ir) dofstr = "ISDOF";
// 		}
// 		if( 0==nodemap.count(res_to_nodename[ir]) ) nodemap[res_to_nodename[ir]] = new Node(res_to_nodename[ir],dofstr);
// 	}
// 	// set tree topology
// 	for(int i = 1; i <= symm_conf.fold_tree().num_jump(); ++i) {
// 		std::string up = res_to_nodename[symm_conf.fold_tree().  upstream_jump_residue(i)];
// 		std::string dn = res_to_nodename[symm_conf.fold_tree().downstream_jump_residue(i)];
// 		if(up==dn) utility_exit_with_message("BAD FOLD TREE NODES!!!!!!!");
// 		nodemap[dn]->setparent(nodemap[up]);
// 	}
// 	// sanity check: make sure tree is connected
// 	for(std::map<std::string,Node*>::const_iterator i = nodemap.begin(); i != nodemap.end(); ++i) {
// 		std::string rootname0 = nodemap.begin()->second->root()->name;
// 		if( rootname0 != i->second->root()->name ) utility_exit_with_message("Nodes not connected!!!");
// 	}
// 	return nodemap.begin()->second->root()->str();
// }

//////////////////////////// END sheffler visualize fold tree

///@brief remodel a fold tree to account for a large insertion by adding the size of the insert to upstream positions
///@author Steven Lewis smlewi@gmail.com as a favor for Jared
core::kinematics::FoldTree
remodel_fold_tree_to_account_for_insertion(
	core::kinematics::FoldTree const & input_tree, //return a remodeled version of this tree
	core::Size insert_after, //add insert_size to points after this in primary sequence in the tree
	core::Size insert_size){

	if(input_tree.is_jump_point(insert_after)){
		throw utility::excn::EXCN_Msg_Exception("FoldTree utility remodel_fold_tree_to_account_for_insertion: I do not know how to handle insertion points that are also jump points - does the jump stay where it was or move to the end of the insert?");
	}

	core::kinematics::FoldTree return_tree;

	typedef std::vector< core::kinematics::Edge > EdgeList; //I am not responsible for the std::vector!
	typedef EdgeList::const_iterator ELconst_iterator;

	for( ELconst_iterator it(input_tree.begin()), end(input_tree.end()); it!=end; ++it){
		//get a copy of the old Edge's start/stop, and update them as necessary
		core::Size start(it->start()), stop(it->stop());
		if(start>insert_after) start = start+insert_size;
		if(stop>insert_after)  stop  = stop+insert_size;

		//copy old edge to new edge
		core::kinematics::Edge const new_edge(
			start,
			stop,
			it->label(),
			it->start_atom(),
			it->stop_atom(),
			it->keep_stub_in_residue());

		//put in new fold tree
		return_tree.add_edge(new_edge);
	}

	//return_tree.reorder(input_tree.root()); //I am not convinced this is necessary or useful here but welcome input on the topic
	return return_tree;
}

} // namespace kinematics
} // namespace core
