// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file   core/kinematics/util.cc
/// @brief  Kinematics utility functions
/// @author Phil Bradley
// Unit headers
#include <core/kinematics/util.hh>
// Package headers
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/JumpAtom.hh>
#include <core/kinematics/tree/BondedAtom.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/types.hh>
#include <basic/Tracer.hh>
#include <utility>
#include <utility/excn/Exceptions.hh>
// Numeric headers
#include <numeric/random/random.hh>
// C++ headers
#include <utility/assert.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/format.hh>
namespace core {
namespace kinematics {
using namespace id;
static basic::Tracer TR( "core.kinematics.util" );
///////////////////////////////////////////////////////////////////////////////
// wrapper for the add_atom recursive function which doesn't know anything
// about Residues
//
/// @details recursively called until all atoms in this residue are added.
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
	AtomOP const atom_p( add_jump_atom ? static_cast< Atom* >( new JumpAtom() ) : static_cast< Atom* >( new BondedAtom() ) );
	// fill in the atom_ptr data
	debug_assert( atom_ptr[ atomno ] == nullptr );
	atom_ptr[ atomno ] = atom_p;
	// set the atom_id information
	atom_p->id( AtomID( atomno, seqpos ) );
	atom_p->parent( AtomAP() );
	utility::vector1< Size > const & nbrs( links[ atomno ] );
	for ( Size i = 1, i_end = nbrs.size(); i <= i_end; ++i ) {
		int const nbr( nbrs[i] );
		if ( atom_ptr[ nbr ] != nullptr ) continue;
		atom_p->append_atom( add_atom( nbr, seqpos, links, atom_ptr, false ) );
	}
	return atom_p;
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
		if ( r > cut_bias_sum( i - 1 ) && r <= cut_bias_sum( i ) ) {
			cutpoint = i;
		}
	}
	if ( cutpoint == 0 ) {
		TR.Warning << "pick loopy cutpoint = 0! setting = 1" << std::endl;
		cutpoint = 1;
	}
	return cutpoint;
}

core::Size
jump_which_partitions( FoldTree const & fold_tree, utility::vector1< bool > residues ) {

	if ( residues.size() != fold_tree.nres() ) {
		utility_exit_with_message("Vector passed to jump_which_partitions() must be the same length as the FoldTree.");
	}

	for ( core::Size jj(1); jj <= fold_tree.num_jump(); ++jj ) {
		bool valid_jump( true );

		utility::vector1< bool > partition( fold_tree.partition_by_jump(jj) );
		debug_assert( residues.size() == partition.size() );

		bool parity( residues[ 1 ] == partition[ 1 ] ); // Are we looking for all same ( true ) or different ( false ) values?

		for ( core::Size ii(2); ii <= partition.size(); ++ii ) {
			bool match( residues[ ii ] == partition[ ii ] );
			if ( parity != match ) {
				valid_jump = false;
				break;
			}
		}

		if ( valid_jump ) {
			return jj; // The assumption there's only one such jump is a good one.
		}
	}

	return 0;
}

FoldTree
get_foldtree_which_partitions( FoldTree const & fold_tree, utility::vector1< bool >  residues ) {

	// Is there a better way of doing this with a recursive algorithm?

	debug_assert( fold_tree.nres() == residues.size() );

	FoldTree new_fold_tree( fold_tree.nres() );
	new_fold_tree.clear();

	Size original_root( fold_tree.root() ); // The root of the original fold_tree
	Size other_root(0); // The root for the other side

	bool root_color( residues[ original_root ] );
	core::Size dividing_jump(0); // The jump which connects the two sets

	utility::vector1< bool > jumps_transfered( fold_tree.num_jump(), false ); // Which jumps have been transfered - we'll have at least this many jumps

	// First transfer all of the consistent jumps ( or the first in/out jump )
	for ( core::Size jj(1); jj <= fold_tree.num_jump(); ++jj ) {
		Edge const & edge( fold_tree.jump_edge(jj) );
		if ( residues[ edge.start() ] == residues[ edge.stop() ] ) { // On same side of the set
			TR << "Transfering over " << edge << std::endl;
			new_fold_tree.add_edge( edge );
			jumps_transfered[ jj ] = true;
		} else if ( edge.start() == original_root && dividing_jump == 0 ) {
			// This is a jump starting at the original and going across the divide (the lowest numbered such.)
			TR << "Transfering over as dividing jump " << edge << std::endl;
			debug_assert( residues[ edge.stop() ] != root_color );
			new_fold_tree.add_edge( edge );
			jumps_transfered[ jj ] = true;
			other_root = edge.stop();
			dividing_jump = jj;
		} else {
			// Defer the rest to a second pass
			TR << "Defering Jump processing " << edge << std::endl;
			TR << "    " << edge.start() << "  " << original_root << " " << dividing_jump << std::endl;
		}
	}

	TR << "other root: " << other_root << " dividing jump: " << dividing_jump << std::endl;

	// Okay, now we can transfer the inconsistent jumps
	//
	// Note: The way we handle inconsistent jumps is not necessarily ideal, as an in1->out1->in2->out2->in3
	// style jump pattern will re-root in3 on in1 rather than in2 ... but handling that properly is
	// more complexity than I'm interested in doing at the moment.
	for ( core::Size jj(1); jj <= fold_tree.num_jump(); ++jj ) {
		TR << "Looking at jump " << jj << " again " << jumps_transfered[ jj ] << std::endl;
		if ( jumps_transfered[ jj ] == true ) { continue; }
		Edge const & edge( fold_tree.jump_edge(jj) );
		debug_assert( residues[ edge.start() ] != residues[ edge.stop() ] ); // Should have dealt with all the consistent ones
		TR << "Processing deferred Jump " << edge << std::endl;
		if ( residues[ edge.stop() ] == root_color ) {
			TR << "Rehoming jump to root " << edge << std::endl;
			new_fold_tree.add_edge( original_root, edge.stop(), jj );
			jumps_transfered[ jj ] = true;
		} else {
			if ( other_root == 0 ) {
				// We didn't transfer a jump that's from the root - we can just re-use the current one.
				TR << "Taking different jump to make new non-root root " << edge << std::endl;
				new_fold_tree.add_edge( edge );
				jumps_transfered[ jj ] = true;
				other_root = edge.stop();
				dividing_jump = jj;
			} else {
				TR << "Rehoming jump to non-root root " << edge << std::endl;
				new_fold_tree.add_edge( other_root, edge.stop(), jj );
				jumps_transfered[ jj ] = true;
			}
		}
	}

	TR << "other root: " << other_root << " dividing jump: " << dividing_jump << std::endl;

	debug_assert( std::all_of(jumps_transfered.begin(), jumps_transfered.end(), [](bool b) { return b; }) ); // All true?

	core::Size max_jump( jumps_transfered.size() );

	// Now transfer all non-jump edges
	for ( Edge const & edge: fold_tree ) {
		if ( edge.label() > 0 ) {
			// Jump edge -- we should have already handled this.
			continue;
		} else if ( edge.label() == Edge::CHEMICAL ) {
			// TODO: I (RM) don't know enough about how to sensibly re-do a Chemical edge across the interface
			// So for now, if we encounter this die. (We can handle consistent chemical edges, though.)
			if ( residues[ edge.start() ] != residues[ edge.stop() ] ) {
				// Handle this sanely, at some point.
				utility_exit_with_message("Error: cannot currently reorganize a FoldTree in a way that splits a chemical edge!");
			} else {
				// Consistent chemical edge - just copy over.
				TR << "Transfering over chemical edge" << edge << std::endl;
				new_fold_tree.add_edge( edge );
			}
		} else if ( edge.label() == Edge::PEPTIDE ) {
			// Ideal case - the entire peptide is all the same side;
			bool consistent = true;
			for ( core::Size ii( std::min( edge.start(), edge.stop() ) ), ii_end( std::max( edge.start(), edge.stop() ) );
					ii <= ii_end; ++ii ) {
				if ( residues[ ii ] !=  residues[ edge.start() ] ) {
					consistent = false;
					break;
				}
			}
			if ( consistent ) {
				// consistent edge - just copy over.
				TR << "Transfering over peptide edge " << edge << std::endl;
				new_fold_tree.add_edge( edge );
			} else {
				// Have to make new edges
				int dir = ( edge.start() <= edge.stop() ) ? 1 : -1;
				Size new_edge_start( edge.start() );
				Size off_edge_end( 0 );
				// We've already built the start with the jump here. Figure out the rest of the items
				for ( Size ii( edge.start() + dir ); ii != (edge.stop() + dir); ii += dir ) {
					if ( residues[ ii ] == residues[ ii - dir ] ) { // Same side as previous residue
						if ( ii == edge.stop() && ii != new_edge_start ) {
							// One last peptide edge
							TR << "Making new peptide edge " << edge << std::endl;
							new_fold_tree.add_edge( new_edge_start, ii, Edge::PEPTIDE );
						}
						continue;
					}
					// Transition. Need to make an edge so far, and then a jump to the new residue
					Size new_edge_end( ii - dir );
					if ( new_edge_end != new_edge_start  ) { // Don't make peptide edge for single residue
						TR << "Making new peptide edge " << edge << std::endl;
						new_fold_tree.add_edge( new_edge_start, new_edge_end, Edge::PEPTIDE );
					}
					// Make jump to the current residue
					if ( off_edge_end != 0 ) {
						TR << "Making new mid-peptide jump " << edge << std::endl;
						new_fold_tree.add_edge( off_edge_end, ii, ++max_jump );
					} else {
						// First off-edge of peptide-edge
						if ( residues[ ii ] == root_color ) {
							TR << "Making new mid-peptide jump to root " << edge << std::endl;
							new_fold_tree.add_edge( original_root, ii, ++max_jump );
						} else {
							// Have to catch the possibility that there were no direct jumps to the non-root color
							if ( other_root == 0 ) {
								TR << "Making new mid-peptide jump as dividing jump " << edge << std::endl;
								new_fold_tree.add_edge( original_root, ii, ++max_jump );
								other_root = ii;
								dividing_jump = max_jump;
							} else {
								TR << "Making new mid-peptide jump to non-root root " << edge << std::endl;
								new_fold_tree.add_edge( other_root, ii, ++max_jump );
							}
						}
					}
					new_edge_start = ii;
					off_edge_end = new_edge_end;
				}
			}
			// TODO: XXX
		} else {
			utility_exit_with_message( "Don't know how to handle edge of type " + utility::to_string( edge.label() ) );
		}
	}

	TR << "other root: " << other_root << " dividing jump: " << dividing_jump << std::endl;

	debug_assert( dividing_jump != 0 );

	// We should have all the edges transfered over now: check the new FoldTree for sanity and return.
	new_fold_tree.reorder( original_root ); // Reorder does some sanity checks/fixes
	runtime_assert( new_fold_tree.connected() );
	runtime_assert( new_fold_tree.check_fold_tree() );

	TR << "Converted a FoldTree with " << fold_tree.size() << " edges (" << fold_tree.num_jump() << " jumps) into one with "
		<< new_fold_tree.size() << " edges (" << new_fold_tree.num_jump() << " jumps) where jump " << dividing_jump << " splits the tree as desired." << std::endl;
	return new_fold_tree;
}


/// @brief Return a list of residue numbers which are on the upstream side of the jump.
utility::vector1< core::Size >
residues_upstream_of_jump( FoldTree const & fold_tree, core::Size jump_id ) {

	utility::vector1< core::Size > reslist;

	utility::vector1< bool > partition( fold_tree.partition_by_jump( jump_id ) );

	bool upstream_color( partition[ fold_tree.upstream_jump_residue( jump_id ) ] );

	for ( core::Size ii(1); ii <= partition.size(); ++ii ) {
		if ( partition[ ii ] == upstream_color ) {
			reslist.push_back( ii );
		}
	}

	return reslist;
}

/// @brief Return a list of residue numbers which are on the downstream side of the jump.
utility::vector1< core::Size >
residues_downstream_of_jump( FoldTree const & fold_tree, core::Size jump_id ) {

	utility::vector1< core::Size > reslist;

	utility::vector1< bool > partition( fold_tree.partition_by_jump( jump_id ) );

	bool downstream_color( partition[ fold_tree.downstream_jump_residue( jump_id ) ] );

	for ( core::Size ii(1); ii <= partition.size(); ++ii ) {
		if ( partition[ ii ] == downstream_color ) {
			reslist.push_back( ii );
		}
	}

	return reslist;
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
		if ( !special ) out << "*";
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
		if ( !special ) {
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
		if ( !special ) {
			out << "*";
			move_chi.push_back( mm.get_chi( pos ) ? '*' : 'x' );
			move.push_back( mm.get_bb( pos ) ? '*' : 'x' );
		}
	}
	out << "\n" << move << "\n" << move_chi;
	out << std::endl;
}

/// @brief linearizes (or defoliates, if you prefer) a FoldTree.  "default" FoldTrees produced by the PDB reader have all chains (peptide edges) starting from jumps relative to residue 1.  This code modifies the tree to instead have all the jumps be relative to the preceding edge.  It is not tested with ligands and will not work with "functional" jumps.  From A to B:
///A:FOLD_TREE  EDGE 1 78 -1  EDGE 1 79 1   EDGE 79 454 -1  EDGE 1 455 2    EDGE 455 540 -1  EDGE 1 541 3    EDGE 541 697 -1
///B:FOLD_TREE  EDGE 1 78 -1  EDGE 78 79 1  EDGE 79 454 -1  EDGE 454 455 2  EDGE 455 540 -1  EDGE 540 541 3  EDGE 541 697 -1
core::kinematics::FoldTree
linearize_fold_tree( core::kinematics::FoldTree const & tree ) {
	core::kinematics::FoldTree newtree;
	for ( auto const & it : tree ) {
		//if it is not a jump, we don't modify it
		if ( !it.is_jump() ) newtree.add_edge(it);
		//if it is a jump, we move start() to stop-1.  This is naive but works for the intended case.
		else newtree.add_edge(core::kinematics::Edge(it.stop()-1, it.stop(), it.label()));
	}
	return newtree;
}

////////////////////////// sheffler visualize FT ////////////////////////////////////
void
replace_substr( std::string& str, const std::string & from, const std::string & to ) {
	size_t start_pos = 0;
	while ( ( start_pos = str.find( from, start_pos ) ) != std::string::npos ) {
		str.replace( start_pos, from.length(), to );
		start_pos += to.length(); // In case 'to' contains 'from', like replacing 'x' with 'yx'
	}
}

std::string
operator*( std::string const & s, size_t n ) {
	std::string r; // empty string
	r.reserve( n * s.size() );
	for ( size_t i = 0; i < n; i++ ) r += s;
	return r;
}

std::string
pad_dash( Size npad, std::string s ) {
	// return s;
	int topad = npad-s.size();
	if ( topad > 0 && topad % 2 == 0 ) return std::string( "-" ) * ( topad / 2    ) + s + std::string( "-" ) * ( topad / 2 );
	if ( topad > 0 && topad % 2 == 1 ) return std::string( "-" ) * ( topad / 2 + 1) + s + std::string( "-" ) * ( topad / 2 );
	else return s;
}

std::string
pad_dash_left ( Size npad, std::string s ) {
	// return s;
	int topad = npad - s.size();
	if ( topad > 0 ) return s + std::string( "-" ) * topad;
	else return s;
}

std::string
pad_dash_right ( Size npad, std::string s ) {
	// return s;
	int topad = npad - s.size();
	if ( topad > 0 ) return std::string( "-" ) * topad + s;
	else return s;
}

struct Node;
using NodeOP = utility::pointer::shared_ptr<Node>;
using NodeAP = utility::pointer::weak_ptr<Node>;

struct Node {
	Node( std::string const & _name, Size _jnum, Size _jumpfrom, Size _jumpto, char _jumpmark = ( char )NULL, Size _follows = 0 )
	: name( _name ), jnum( _jnum ), jumpfrom( _jumpfrom ), jumpto( _jumpto ), prefix_len( 8 ), follows( _follows ), jumpmark( _jumpmark ), parent() {}

	~Node() = default;

	void
	setparent( NodeOP p ) {
		parent = p;
		parent.lock()->children.push_back( this );
	}

	Node &
	root() {
		return parent.lock() == nullptr ? *this : parent.lock()->root();
	}

	std::string
	str() {
		using ObjexxFCL::string_of;
		using std::string;
		string s;
		for ( utility::vector1< Node* >::const_iterator ic = children.begin(); ic != children.end(); ++ic ) {
			string mark = ( ( *ic )->jumpmark == ( char ) NULL ) ? "-" : string( "" ) + ( *ic )->jumpmark;
			string vchar = children.size() == 1 ? "" : ( *ic == children.back() ? " "  : "|" );
			string schar = children.size() == 1 ? "" : ( *ic == children.back() ? "\\" : "|" );
			std::string uprsd = ( ( *ic )->jumpfrom != 0 ) ? pad_dash_left( 3, string_of( ( *ic )->jumpfrom ) ) : std::string( "" );
			std::string dnrsd = ( ( *ic )->jumpto   != 0 ) ? pad_dash_right( 4, ">" + string_of( ( *ic )->jumpto ) ) + ":" : std::string( ">" );
			std::string jstr;
			if ( ( *ic )->follows != 0 ) { // no mark if follows
				jstr = "j" + string_of( ( *ic )->jnum ) + ( ( ( *ic )->follows ) ? "=" + string_of( ( *ic )->follows ) : "" );
			} else {
				jstr = mark + "j" + string_of( ( *ic )->jnum ) + mark;
			}
			std::string prefix = schar + uprsd + "--" + pad_dash( prefix_len, jstr ) + "--" + dnrsd;
			string news = ( *ic )->str();
			string pad = string( " " ) * ( prefix.size() - 1 );
			if ( children.size() == 1 ) pad += string( " " ) * ( name.size() + 1 );
			replace_substr( news, "\n", "\n" + vchar + pad );
			if ( children.size() == 1 ) s +=        prefix + news;
			else                        s += "\n" + prefix + news;
		}
		s = name + s;
		return s;
	}
	std::string name;
	Size jnum, jumpfrom, jumpto, prefix_len, follows;
	char jumpmark;
private:
	NodeAP parent;
	utility::vector1< Node * > children;
};


struct TreeVizBuilder {
	utility::vector1< Size > lb_, ub_;
	FoldTree const & ft;

	explicit TreeVizBuilder( core::kinematics::FoldTree const & _ft ) : ft( _ft ) {
		lb_.resize( ft.nres(), 0 );
		ub_.resize( ft.nres(), 0 );
	}

	void
	get_ft_node_bounds( Size res, Size & out_lb, Size & out_ub ) {
		if ( lb_[ res ] == 0 ) {
			lb_[ res ] = 1;
			ub_[ res ] = ft.nres();
			for ( Size i = 1; i <= ft.num_cutpoint(); ++i ) {
				auto c = ( Size ) ft.cutpoint( i );
				if ( c <  res ) lb_[ res ] = std::max( lb_[ res ], c + 1 );
				if ( c >= res ) ub_[ res ] = std::min( ub_[ res ], c     );
			}
		}
		out_lb = lb_[ res ];
		out_ub = ub_[ res ];
	}

	Size
	get_ft_node_lower_bound( Size res ) {
		Size lb, ub;
		get_ft_node_bounds( res, lb, ub );
		return lb;
	}

	Size
	is_single( Size res ) {
		Size lb, ub;
		get_ft_node_bounds( res, lb, ub );
		return lb == ub;
	}

	Size
	get_ft_node_subroot( Size res ) {
		Size lb, ub;
		get_ft_node_bounds( res, lb, ub);
		for ( Size i = 1; i <= ft.num_jump(); ++i ) {
			Size dn = ft.downstream_jump_residue( i );
			if ( lb <= dn && dn <= ub ) return dn;
		}
		return lb; // default subroot is first res
	}

	void
	expand_node_labels_partial_by_contig( std::map< Size, std::string > & node_labels_partial ) {
		utility::vector1< Size > tocheck;
		for ( auto & i : node_labels_partial ) {
			tocheck.push_back( i.first );
		}
		for ( utility::vector1< Size >::const_iterator i = tocheck.begin(); i != tocheck.end(); ++i ) {
			for ( utility::vector1< Size >::const_iterator j = tocheck.begin(); j != tocheck.end(); ++j ) {
				if ( *i == *j ) continue;
				// TR << "check " << *i << " " << *j << endl;
				if ( get_ft_node_lower_bound( *i ) == get_ft_node_lower_bound( *j ) ) {
					if ( node_labels_partial[ *i ] != node_labels_partial[ *j ] ) {
						utility_exit_with_message( "non-matching node_labels_partial requested in same FT contig" );
					}
				}
			}
			Size lb, ub;
			get_ft_node_bounds( *i, lb, ub );
			for ( Size ir = lb; ir <= ub; ++ir ) {
				node_labels_partial[ ir ] = node_labels_partial[ *i ];
			}
		}
	}

	utility::vector1<std::string>
	get_res_nodenames( std::map< Size, std::string > node_labels_partial ) {
		using ObjexxFCL::format::I;
		//int npad = 1;
		//if ( ft.nres() >    9 ) npad = 2;
		//if ( ft.nres() >   99 ) npad = 3;
		//if ( ft.nres() >  999 ) npad = 4;
		//if ( ft.nres() > 9999 ) npad = 5;
		int npad = 0;
		utility::vector1<std::string> names(ft.nres(),"ERROR_this_name_was_not_set");
		for ( Size i = 1; i <= ft.nres(); ++i ) {
			Size lb, ub;
			get_ft_node_bounds( i, lb, ub );
			std::string lbl = node_labels_partial.count( lb ) ? node_labels_partial[ lb ] : "Contig";

			std::string resrange = I( npad, lb ) + ( ( ub == lb ) ? "" : "-" + I( npad, ub ) );
			names[ i ] = lbl + "(" + resrange + ")";
		}
		return names;
	}

	Size
	get_jump_num_to_contig_of_resi( Size resi ) {
		Size lb, ub;
		get_ft_node_bounds( resi, lb, ub );
		for ( Size i = 1; i <= ft.num_jump(); ++i ) {
			if ( lb <= ( Size )ft.downstream_jump_residue( i ) && ( Size )ft.downstream_jump_residue( i ) <= ub ) return i;
		}
		return 0;
	}
};

std::string
visualize_fold_tree(
	FoldTree const & ft,
	std::map< Size, std::string > const & node_labels_partial_in,
	std::map< Size, char > const & mark_jump_to_res,
	std::map< Size, Size > const & jump_follows
){
	TreeVizBuilder tvb( ft );
	std::map< Size, std::string > node_labels_partial( node_labels_partial_in );
	tvb.expand_node_labels_partial_by_contig( node_labels_partial );
	utility::vector1<std::string> res_to_nodename = tvb.get_res_nodenames( node_labels_partial );
	// make Node array -- OPs here so when the map goes out of scope we don't leak memory
	std::map< std::string, NodeOP > nodemap;
	for ( Size ir = 1; ir <= ft.nres(); ++ir ) {
		if ( nodemap.count(res_to_nodename[ ir ] ) ) continue;
		Size lb, ub;
		tvb.get_ft_node_bounds( ir, lb, ub );
		Size jnum = tvb.get_jump_num_to_contig_of_resi( ir );
		Size jumpfrom = ( jnum && !tvb.is_single( ft.upstream_jump_residue( jnum ) ) ) ? ft.  upstream_jump_residue( jnum ) : 0;
		Size jumpto   = ( jnum && !tvb.is_single( ft.downstream_jump_residue( jnum ) ) ) ? ft.downstream_jump_residue( jnum ) : 0;
		auto markjump = ( char )NULL;
		for ( Size jr = lb; jr <= ub; ++jr ) if ( mark_jump_to_res.count( jr ) ) markjump = mark_jump_to_res.find( jr )->second;
		Size follows = jump_follows.find( jnum ) != jump_follows.end() ? jump_follows.find( jnum )->second : ( Size ) 0;
		nodemap[ res_to_nodename[ ir ] ] = NodeOP( new Node( res_to_nodename[ ir ], jnum, jumpfrom, jumpto, markjump, follows ) );
	}
	// set tree topology
	for ( Size i = 1; i <= ft.num_jump(); ++i ) {
		std::string up = res_to_nodename[   ft.upstream_jump_residue( i ) ];
		std::string dn = res_to_nodename[ ft.downstream_jump_residue( i ) ];
		if ( up == dn ) {
			std::cerr << "bad jump " << i << " from " << ft.upstream_jump_residue( i ) << " to " << ft.downstream_jump_residue( i ) << std::endl;
			utility_exit_with_message( "BAD FOLD TREE NODES!!!!!!!" );
		}
		nodemap[ dn ]->setparent( nodemap[ up ] );
	}
	// sanity check: make sure tree is connected. This only works if the fold tree has no chemical edges
	if ( ft.get_chemical_edges().size() == 0 ) {
		for ( std::map< std::string,NodeOP >::const_iterator i = nodemap.begin(); i != nodemap.end(); ++i ) {
			std::string rootname0 = nodemap.begin()->second->root().name;
			if ( rootname0 != i->second->root().name ) utility_exit_with_message( "Nodes not connected!!!" );
			// std::cerr << "=========================== " << i->second->name << " ===========================" << std::endl;
			// std::cerr << i->second->str() << std::endl;
			// std::cerr << "========================================================================" << std::endl;
		}
	}
	return nodemap.begin()->second->root().str();
}

std::string
visualize_fold_tree( FoldTree const & ft ) {
	std::map< Size, std::string > empty_labels;
	std::map< Size, char > empty_marks;
	std::map< Size, Size > empty_follows;
	return visualize_fold_tree( ft, empty_labels, empty_marks, empty_follows );
}

std::string
visualize_fold_tree( FoldTree const & ft, std::map< Size, std::string > const & node_labels_partial ) {
	std::map< Size, char > empty_marks;
	std::map< Size, Size > empty_follows;
	return visualize_fold_tree( ft, node_labels_partial, empty_marks, empty_follows );
}

std::string
visualize_fold_tree( FoldTree const & ft, std::map<Size,char> const & mark_jump_to_res ) {
	std::map< Size, std::string > empty_labels;
	std::map< Size, Size > empty_follows;
	return visualize_fold_tree( ft, empty_labels, mark_jump_to_res, empty_follows );
}

//////////////////////////// END sheffler visualize fold tree

/// @brief remodel a fold tree to account for a large insertion by adding the size of the insert to upstream positions
/// @author Steven Lewis smlewi@gmail.com as a favor for Jared
core::kinematics::FoldTree
remodel_fold_tree_to_account_for_insertion(
	core::kinematics::FoldTree const & input_tree, //return a remodeled version of this tree
	core::Size insert_after, //add insert_size to points after this in primary sequence in the tree
	core::Size insert_size){
	if ( input_tree.is_jump_point(insert_after) ) {
		throw CREATE_EXCEPTION(utility::excn::Exception, "FoldTree utility remodel_fold_tree_to_account_for_insertion: I do not know how to handle insertion points that are also jump points - does the jump stay where it was or move to the end of the insert?");
	}
	core::kinematics::FoldTree return_tree;
	for ( auto const & it : input_tree ) {
		//get a copy of the old Edge's start/stop, and update them as necessary
		core::Size start(it.start()), stop(it.stop());
		if ( start>insert_after ) start = start+insert_size;
		if ( stop>insert_after )  stop  = stop+insert_size;
		//copy old edge to new edge
		core::kinematics::Edge const new_edge(
			start,
			stop,
			it.label(),
			it.start_atom(),
			it.stop_atom(),
			it.keep_stub_in_residue());
		//put in new fold tree
		return_tree.add_edge(new_edge);
	}
	//return_tree.reorder(input_tree.root()); //I am not convinced this is necessary or useful here but welcome input on the topic
	return return_tree;
}

/// @brief Get a vector of residues matching the id from a movemap.
utility::vector1< core::Size >
get_residues_from_movemap_with_id( id::TorsionType query_torsion, MoveMap const & movemap){
	utility::vector1< core::Size > residues;

	using TorsionType = id::TorsionType;
	typedef std::pair< Size, TorsionType > MoveMapTorsionID;

	for ( auto it = movemap.movemap_torsion_id_begin(), it_end = movemap.movemap_torsion_id_end();
			it != it_end; ++it ) {
		MoveMapTorsionID mmtorsionID = it->first;
		Size res = mmtorsionID.first;
		TorsionType torsiontype = mmtorsionID.second;

		// Jumps are handled under a separate heading.
		if ( torsiontype == query_torsion && it->second == true ) {
			residues.push_back( res );
			continue;
		}
	}
	return residues;
}

utility::vector1< core::Size >
get_residues_from_movemap_bb_any_torsion(MoveMap const & movemap, Size total_resnum){
	utility::vector1< bool > bb_on(total_resnum, false);
	utility::vector1< core::Size > final_vec;

	for ( core::Size i = 1; i <= total_resnum; ++i ) {
		if ( movemap.get_bb( i ) ) {
			bb_on[ i ] = true;
			continue;
		}

		for ( core::Size x = 1; x <= 4; ++x ) {
			if ( movemap.get_bb( i, x ) ) {
				bb_on[ i ] = true;
				break;
			}
		}
	}

	//Probably a better way for this, but I don't know it.
	for ( core::Size resnum = 1; resnum <= total_resnum; ++resnum ) {
		if ( bb_on[ resnum ] ) {
			final_vec.push_back( resnum );
		}
	}
	return final_vec;

}

utility::vector1< core::Size >
get_residues_from_movemap_bb_or_chi(MoveMap const & movemap, Size total_resnum){

	utility::vector1< bool > on(total_resnum, false);
	utility::vector1< core::Size > final_vec;

	for ( core::Size i = 1; i <= total_resnum; ++i ) {
		if ( movemap.get_bb( i )  || movemap.get_chi( i ) ) {
			on[ i ] = true;
			continue;
		}

		//Torsion specific
		for ( core::Size x = 1; x <= 4; ++x ) {
			if ( movemap.get_bb( i, x ) ) {
				on[ i ] = true;
				break;
			}
		}
	}

	for ( core::Size resnum = 1; resnum <= total_resnum; ++resnum ) {
		if ( on[ resnum ] ) {
			final_vec.push_back( resnum );
		}
	}
	return final_vec;
}

} // namespace kinematics
} // namespace core
