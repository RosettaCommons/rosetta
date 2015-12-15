// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*- // vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/FoldGraph.cc
/// @brief FoldGraph - a fold-tree which takes two-way dependencies into account
/// @detailed
/// @author Tom Linsky

//Unit Headers
#include <protocols/denovo_design/components/FoldGraph.hh>

// Project Headers
#include <protocols/denovo_design/components/Segment.hh>

// Protocol Headers
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Core Headers
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>

//Basic/Utility/Numeric Headers
#include <basic/Tracer.hh>

//C++ Headers

static basic::Tracer TR("protocols.denovo_design.components.FoldGraph");

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace components {

///////////////////////////////////////////////////////////////////////////////
/// FOLDGRAPH CLASS METHODS
///////////////////////////////////////////////////////////////////////////////

FoldGraph::FoldGraph( StructureData const & perm, core::pose::PoseCOP pose ) :
	utility::pointer::ReferenceCount()
{
	TR.Debug << "Constructing foldgraph for " << perm.id() << " : " << perm  << std::endl;
	TR.Debug << "Chains:" << std::endl;
	for ( core::Size i=1; i<=pose->total_residue(); ++i ) {
		TR.Debug << pose->residue( i ).name() << " " << i << " : " << pose->chain( i ) << std::endl;
	}
	g_.drop_all_edges();
	gpeptide_.drop_all_edges();
	seg2node_.clear();
	node2seg_.clear();

	//must have at least one segment
	assert( std::distance( perm.segments_begin(), perm.segments_end() ) );

	// populate seg2node
	for ( StringList::const_iterator r = perm.segments_begin(), end = perm.segments_end(); r != end; ++r ) {
		if ( seg2node_.find( *r ) == seg2node_.end() ) {
			core::Size const nodenum = seg2node_.size() + 1;
			seg2node_[*r] = nodenum;
			TR.Debug << "Added " << *r << std::endl;
		}
	}

	// populate node2seg
	for ( std::map< std::string, core::Size >::const_iterator c2n = seg2node_.begin(); c2n != seg2node_.end(); ++c2n ) {
		node2seg_[c2n->second] = c2n->first;
	}
	assert( seg2node_.size() );
	assert( node2seg_.size() == seg2node_.size() );
	g_.set_num_nodes( seg2node_.size() );
	gpeptide_.set_num_nodes( seg2node_.size() );

	// add edges between movable groups
	std::set< core::Size > const movable_groups = perm.movable_groups();
	TR.Debug << "Num groups = " << movable_groups.size() << std::endl;
	assert( movable_groups.size() );
	for ( std::set< core::Size >::const_iterator grp = movable_groups.begin(); grp != movable_groups.end(); ++grp ) {
		utility::vector1< std::string > segments = perm.segments_in_movable_group( *grp );
		for ( core::Size i=1; i<=segments.size(); ++i ) {
			TR.Debug << "Group " << *grp << " Segment " << segments[i] << std::endl;
			for ( core::Size j=i+1; j<=segments.size(); ++j ) {
				if ( segments[i] != segments[j] ) {
					add_edge( segments[i], segments[j] );
				}
			}
		}
	}

	// look for non-canonical bonds
	for ( utility::vector1< BondInfo >::const_iterator bi = perm.covalent_bonds_begin(); bi != perm.covalent_bonds_end(); ++bi ) {
		TR << "Found non-canonical connection from " << bi->seg1 << ":" << bi->res1 << ":" << bi->atom1
			<<" to " << bi->seg2 << ":" << bi->res2 << ":" << bi->atom2 << std::endl;
		add_edge( bi->seg1, bi->seg2 );
	}

	// add peptide edges
	for ( StringList::const_iterator n = perm.segments_begin(); n != perm.segments_end(); ++n ) {
		Segment const & resis = perm.segment(*n);
		if ( !resis.has_free_lower_terminus() ) {
			TR.Debug << "Segment1 " << *n << " res " << resis.safe() << " chain " << pose->chain( resis.safe() ) << std::endl;
			TR.Debug << "Segment2 " << resis.lower_segment() << " res " << perm.segment( resis.lower_segment() ).safe() << " chain " << pose->chain( perm.segment( resis.lower_segment() ).safe() ) << std::endl;
			if ( pose->chain( resis.safe() ) == pose->chain( perm.segment( resis.lower_segment() ).safe() ) ) {
				add_peptide_edge( *n, resis.lower_segment() );
			}
		}
		if ( !resis.has_free_upper_terminus() ) {
			TR.Debug << "Segment1 " << *n << " res " << resis.safe() << " chain " << pose->chain( resis.safe() ) << std::endl;
			TR.Debug << "Segment2 " << resis.upper_segment() << " res " << perm.segment( resis.upper_segment() ).safe() << " chain " << pose->chain( perm.segment( resis.upper_segment() ).safe() ) << std::endl;
			if ( pose->chain( resis.safe() ) == pose->chain( perm.segment( resis.upper_segment() ).safe() ) ) {
				add_peptide_edge( *n, resis.upper_segment() );
			}
		}
	}
}

/// @brief given a segment name, returns the associated graph node index
core::Size
FoldGraph::nodeidx( std::string const & segment ) const
{
	std::map< std::string, core::Size >::const_iterator it = seg2node_.find(segment);
	if ( it == seg2node_.end() ) {
		TR.Error << "Component " << segment << " not found in map: " << seg2node_ << std::endl;
		assert( it != seg2node_.end() );
	}
	return it->second;
}

/// @brief given a node index, returns the associated segment
std::string const &
FoldGraph::segment( core::Size const nodeidx ) const
{
	std::map< core::Size, std::string >::const_iterator it = node2seg_.find(nodeidx);
	if ( it == node2seg_.end() ) {
		TR.Error << "Node idx " << nodeidx << " not found in map: " << seg2node_ << std::endl;
		assert( it != node2seg_.end() );
	}
	return it->second;
}

/// @brief adds an non-peptide edge to the foldgraph to indicate a non-covalent interaction
void
FoldGraph::add_edge( std::string const & comp1, std::string const & comp2 )
{
	if ( ! g_.get_edge_exists( nodeidx(comp1), nodeidx(comp2) ) ) {
		TR.Debug << "Adding edge from " << comp1 << " to " << comp2 << std::endl;
		g_.add_edge( nodeidx(comp1), nodeidx(comp2) );
	}
}

/// @brief adds a peptide edge to the foldgraph, indicating a direct covalent interaction
void
FoldGraph::add_peptide_edge( std::string const & comp1, std::string const & comp2 )
{
	if ( ! gpeptide_.get_edge_exists( nodeidx(comp1), nodeidx(comp2) ) ) {
		TR.Debug << "Adding peptide edge from " << comp1 << " to " << comp2 << std::endl;
		gpeptide_.add_edge( nodeidx(comp1), nodeidx(comp2) );
	}
	add_edge( comp1, comp2 );
}

/// @brief returns true if there is a peptide edge between the two given segments
bool
FoldGraph::has_edge( std::string const & seg1, std::string const & seg2 ) const
{
	return has_edge( nodeidx(seg1), nodeidx(seg2) );
}

/// @brief returns true if there is a peptide edge between the two given nodes
bool
FoldGraph::has_edge( core::Size const n1, core::Size const n2 ) const
{
	return g_.get_edge_exists( n1, n2 ) && ( ! gpeptide_.get_edge_exists( n1, n2 ) );
}

/// @brief returns true if there is a peptide edge between the two given segments
bool
FoldGraph::has_peptide_edge( std::string const & seg1, std::string const & seg2 ) const
{
	return has_peptide_edge( nodeidx(seg1), nodeidx(seg2) );
}

/// @brief returns true if there is a peptide edge between the two given nodes
bool
FoldGraph::has_peptide_edge( core::Size const n1, core::Size const n2 ) const
{
	return gpeptide_.get_edge_exists( n1, n2 );
}

/// @brief gives a fold tree based on the segments in the given permutation
core::kinematics::FoldTree
FoldGraph::fold_tree( StructureData const & perm, utility::vector1< std::string > const & root_segments ) const
{
	debug_assert( root_segments.size() >= 1 );
	core::kinematics::FoldTree ft;
	NodeSet visited;
	std::stack< core::Size > to_search;

	// add to fold tree based on peptide graph
	fold_tree_rec( ft, visited, to_search, perm, root_segments[1], 0, true );
	Segment root_resis = perm.segment(root_segments[1]);
	for ( core::Size i=2; i<=root_segments.size(); ++i ) {
		Segment const & resis = perm.segment(root_segments[i]);
		if ( visited.find( nodeidx(root_segments[i]) ) == visited.end() ) {
			core::kinematics::Edge e( root_resis.safe(), resis.safe(), ft.num_jump()+1 );
			TR.Debug << " adding unconnected segment from peptide graph rooted at " << root_segments[i] << " " << e << std::endl;
			ft.add_edge( e );
			fold_tree_rec( ft, visited, to_search, perm, root_segments[i], 0, true );
		}
	}

	// add to fold tree based on interaction graph -- look for reachable nodes that are not yet visited.
	while ( to_search.size() > 0 ) {
		core::Size const node = to_search.top();
		to_search.pop();
		core::graph::Node const * n = g_.get_node(node);
		assert( n );
		for ( core::graph::Graph::EdgeListConstIter e=n->const_edge_list_begin(); e != n->const_edge_list_end(); ++e ) {
			core::Size const othernode = (*e)->get_other_ind(node);
			if ( visited.find(othernode) != visited.end() ) {
				TR.Debug << "ALready visited " << segment(node) << "__" << segment(othernode) << std::endl;
				continue;
			}
			core::kinematics::Edge edge( perm.segment(segment(node)).safe(), perm.segment(segment(othernode)).safe(), ft.num_jump()+1 );
			utility::vector1< BondInfo >::const_iterator bi = perm.non_peptidic_bond( segment(node), segment(othernode) );
			if ( bi == perm.covalent_bonds_end() ) {
				TR.Debug << " adding unconnected segment from interaction graph rooted at " << segment(othernode) << " " << edge << std::endl;
			} else {
				std::string this_seg = "", other_seg = "", this_atom = "", other_atom = "";
				core::Size this_res = 0, other_res = 0;
				if ( bi->seg1 == segment(node) ) {
					debug_assert( bi->seg2 == segment(othernode) );
					this_seg = bi->seg1;
					other_seg = bi->seg2;
					this_res = bi->res1;
					other_res = bi->res2;
					this_atom = bi->atom1;
					other_atom = bi->atom2;
				} else {
					debug_assert( bi->seg2 == segment(node) );
					this_seg = bi->seg2;
					other_seg = bi->seg1;
					this_res = bi->res2;
					other_res = bi->res1;
					this_atom = bi->atom2;
					other_atom = bi->atom1;
				}
				TR.Debug << " adding non-peptide-bonded segment from " << this_seg << "#" << this_res << "," << this_atom
					<< " --> " << other_seg << "#" << other_res << "," << other_atom << std::endl;
				edge.start() = perm.pose_residue( this_seg, this_res );
				edge.stop() = perm.pose_residue( other_seg, other_res );
				edge.start_atom() = this_atom;
				edge.stop_atom() = other_atom;
				core::kinematics::FoldTree const & ft_const = ft;
				utility::vector1< core::kinematics::Edge > newedges;
				utility::vector1< core::kinematics::Edge > deledges;
				for ( core::kinematics::FoldTree::const_iterator e=ft_const.begin(); e != ft_const.end(); ++e ) {
					if ( e->label() > 0 ) continue;
					if ( ( edge.start() >= e->stop() ) && ( edge.start() <= e->start() ) ) {
						TR << "Removing edge: " << *e << " to insert " << edge << std::endl;
						newedges.push_back( core::kinematics::Edge( e->stop(), edge.start(), core::kinematics::Edge::PEPTIDE ) );
						newedges.push_back( core::kinematics::Edge( edge.start(), e->start(), core::kinematics::Edge::PEPTIDE ) );
						deledges.push_back( *e );
					} else if ( ( edge.start() >= e->start() ) && ( edge.start() <= e->stop() ) ) {
						TR << "2Adding edge: " << *e << " to insert " << edge << std::endl;
						newedges.push_back( core::kinematics::Edge( e->start(), edge.start(), core::kinematics::Edge::PEPTIDE ) );
						newedges.push_back( core::kinematics::Edge( edge.start(), e->stop(), core::kinematics::Edge::PEPTIDE ) );
						deledges.push_back( *e );
					}
					if ( ( edge.stop() >= e->stop() ) && ( edge.stop() <= e->start() ) ) {
						TR << "Removing edge: " << *e << " to insert " << edge << std::endl;
						newedges.push_back( core::kinematics::Edge( e->stop(), edge.stop(), core::kinematics::Edge::PEPTIDE ) );
						newedges.push_back( core::kinematics::Edge( edge.stop(), e->start(), core::kinematics::Edge::PEPTIDE ) );
						deledges.push_back( *e );
					} else if ( ( edge.stop() >= e->start() ) && ( edge.stop() <= e->stop() ) ) {
						TR << "2Adding edge: " << *e << " to insert " << edge << std::endl;
						newedges.push_back( core::kinematics::Edge( e->start(), edge.stop(), core::kinematics::Edge::PEPTIDE ) );
						newedges.push_back( core::kinematics::Edge( edge.stop(), e->stop(), core::kinematics::Edge::PEPTIDE ) );
						deledges.push_back( *e );
					}
				}
				for ( utility::vector1< core::kinematics::Edge >::const_iterator e=newedges.begin(); e != newedges.end(); ++e ) {
					ft.add_edge( *e );
				}
				for ( utility::vector1< core::kinematics::Edge >::const_iterator e=deledges.begin(); e != deledges.end(); ++e ) {
					ft.delete_edge( *e );
				}
			}
			ft.add_edge( edge );
			fold_tree_rec( ft, visited, to_search, perm, segment(othernode), 0, true );
		}
	}

	// add unconnected foldtree elements directly to the start node
	for ( core::Size i=1; i<=seg2node_.size(); ++i ) {
		if ( visited.find(i) == visited.end() ) {
			std::string const & unconnected_seg = segment(i);
			Segment const & resis = perm.segment(unconnected_seg);
			core::kinematics::Edge e( root_resis.safe(), resis.safe(), ft.num_jump()+1 );
			TR.Debug << " unconnected_comp = " << unconnected_seg << " " << e << std::endl;
			ft.add_edge( e );
			fold_tree_rec( ft, visited, to_search, perm, unconnected_seg, 0, true );
		}
	}

	TR.Debug << "Before deleting vertices: " << ft << std::endl;
	ft.delete_extra_vertices();
	TR.Debug << "After deleting vertices: " << ft << std::endl;
	debug_assert( ft.check_fold_tree() );
	return ft;
}

/// @brief recursive function to traverse graphs and build fold tree
/// @param[parent_direction] -1 if the previous edge is a peptide edge going backward, 1 if the previous edge is going forward, and 0 if the previous edge is a jump
void
FoldGraph::fold_tree_rec(
	core::kinematics::FoldTree & ft,
	NodeSet & visited,
	std::stack< core::Size > & node_stack,
	StructureData const & perm,
	std::string const & segment_name,
	int const parent_direction,
	bool const polymer_only ) const
{
	core::Size const nodenum = nodeidx( segment_name );
	if ( visited.find(nodenum) != visited.end() ) {
		return;
	}

	TR.Debug << "Starting traversal at " << segment_name <<  " node " << nodenum << " visited " << visited << std::endl;

	Segment const & resis = perm.segment(segment_name);
	if ( resis.safe() != resis.cterm_resi() ) {
		if ( parent_direction > -1 ) {
			TR.Debug << "Adding " << resis.safe() << " " << resis.cterm_resi() << " -1 for " << segment_name << std::endl;
			ft.add_edge( resis.safe(), resis.cterm_resi(), core::kinematics::Edge::PEPTIDE );
		} else {
			TR.Debug << "Adding " << resis.cterm_resi() << " " << resis.safe() << " -1 for " << segment_name << std::endl;
			ft.add_edge( resis.cterm_resi(), resis.safe(), core::kinematics::Edge::PEPTIDE );
		}
	}
	if ( resis.safe() != resis.nterm_resi() ) {
		if ( parent_direction < 1 ) {
			TR.Debug << "Adding " << resis.safe() << " " << resis.nterm_resi() << " -1 for " << segment_name << std::endl;
			ft.add_edge( resis.safe(), resis.nterm_resi(), core::kinematics::Edge::PEPTIDE );
		} else {
			TR.Debug << "Adding " << resis.nterm_resi() << " " << resis.safe() << " -1 for " << segment_name << std::endl;
			ft.add_edge( resis.nterm_resi(), resis.safe(), core::kinematics::Edge::PEPTIDE );
		}
	}

	visited.insert( nodenum );
	node_stack.push( nodenum );

	// traverse connected edges
	core::graph::Node const * n;
	if ( polymer_only ) {
		n = gpeptide_.get_node(nodenum);
	} else {
		n = g_.get_node(nodenum);
	}
	assert( n );
	for ( core::graph::Graph::EdgeListConstIter e=n->const_edge_list_begin(); e != n->const_edge_list_end(); ++e ) {
		core::Size const othernode = (*e)->get_other_ind(nodenum);
		if ( visited.find(othernode) != visited.end() ) {
			TR.Debug << "Skipping " << othernode << " as it has been visited." << std::endl;
			continue;
		}

		// add connecting edge
		std::string const & otherseg = segment(othernode);
		Segment const & other_resis = perm.segment( otherseg );
		if ( resis.cterm_resi() <= other_resis.nterm_resi() ) {
			int direction = 1;
			core::kinematics::Edge e( resis.cterm_resi(), other_resis.nterm_resi(), core::kinematics::Edge::PEPTIDE );
			if ( !polymer_only ) {
				e = core::kinematics::Edge( resis.safe(), other_resis.safe(), ft.num_jump()+1 );
				direction = 0;
			}
			TR.Debug << "Adding " << e << " for " << segment_name << "-->" << otherseg << std::endl;
			ft.add_edge( e );
			fold_tree_rec( ft, visited, node_stack, perm, otherseg, direction, true );
		} else if ( resis.nterm_resi() >= other_resis.cterm_resi() ) {
			int direction = -1;
			core::kinematics::Edge e( resis.nterm_resi(), other_resis.cterm_resi(), core::kinematics::Edge::PEPTIDE );
			if ( !polymer_only ) {
				e = core::kinematics::Edge( resis.safe(), other_resis.safe(), ft.num_jump()+1 );
				direction = 0;
			}
			TR.Debug << "Adding " << e << " for " << segment_name << "-->" << otherseg << std::endl;
			ft.add_edge( e );
			fold_tree_rec( ft, visited, node_stack, perm, otherseg, direction, true );
		} else {
			TR << "Something is wrong. THere are probably two segments that overlap in residue numbering. Nterm=" << resis.nterm_resi() << " Cterm=" << resis.cterm_resi() << " safe=" << resis.safe() << " otherNterm=" << other_resis.nterm_resi() << " otherCterm=" << other_resis.cterm_resi() << " safe=" << other_resis.safe() << std::endl;
			assert( false );
		}
	}
}

bool
is_non_polymeric_connection(
	StructureData const & perm,
	std::string const & seg1,
	std::string const & seg2 )
{
	for ( utility::vector1< BondInfo >::const_iterator bi = perm.covalent_bonds_begin(); bi != perm.covalent_bonds_end(); ++bi ) {
		if ( ( ( seg1 == bi->seg1 ) && ( seg2 == bi->seg2 ) )
				|| ( ( seg2 == bi->seg1 ) && ( seg1 == bi->seg2 ) ) ) {
			TR << "is_non_polymeric_connection " << seg1 << " " << seg2 << " 1" << std::endl;
			return true;
		}
	}
	TR << "is_non_polymeric_connection " << seg1 << " " << seg2 << " 0" << std::endl;
	return false;
}

template< class T >
void
insert_mg(
	std::map< core::Size, std::list< T > > & fixed_mgs,
	core::Size const mg,
	T const & segname )
{
	typename std::map< core::Size, std::list< T > >::iterator fixed_mg = fixed_mgs.find( mg );
	if ( fixed_mg == fixed_mgs.end() ) {
		fixed_mg = fixed_mgs.insert( std::make_pair( mg, std::list< T >() ) ).first;
	}
	debug_assert( fixed_mg != fixed_mgs.end() );
	fixed_mg->second.push_back( segname );
}

/// @brief checks a solution to ensure that covalently bound segments are not found both inside and outside of loops
bool
FoldGraph::check_solution(
	StructureData const & perm,
	Solution const & solution ) const
{
	std::set< std::string > loop_segments;
	for ( core::Size i=1; i<=solution.size(); ++i ) {
		for ( NodeSet::const_iterator lseg = solution[i].begin(); lseg != solution[i].end(); ++lseg ) {
			std::string const & segmentname = segment( *lseg );
			std::pair< std::set< std::string >::iterator, bool > ins_result = loop_segments.insert( segmentname );
			if ( !ins_result.second ) {
				TR << "Segment " << segmentname << " is present in more than one loop... this should not be possible. Skipping." << std::endl;
				return false;
			}
		}
	}

	// matches a fixed mg with a list of segments
	std::map< core::Size, StringList > fixed_mgs;

	// ensure that the same MG isn't present both inside and outside of loops
	for ( StringList::const_iterator c = perm.segments_begin(); c != perm.segments_end(); ++c ) {
		core::Size const mg = perm.segment(*c).movable_group;
		if ( loop_segments.find(*c) == loop_segments.end() ) {
			insert_mg( fixed_mgs, mg, *c );
		}
	}

	for ( Solution::const_iterator nodeset = solution.begin(); nodeset != solution.end(); ++nodeset ) {
		for ( NodeSet::const_iterator n = nodeset->begin(); n != nodeset->end(); ++n ) {
			std::string const & segname = segment( *n );
			core::Size const mg = perm.segment( segname ).movable_group;
			std::map< core::Size, StringList >::const_iterator fixed_mg = fixed_mgs.find( mg );
			if ( fixed_mg != fixed_mgs.end() ) {
				StringList const connected = perm.connected_segments( segname );
				for ( StringList::const_iterator c=fixed_mg->second.begin(); c!=fixed_mg->second.end(); ++c ) {
					// check to see if fixed segment is connected to anything connected to *n
					if ( std::count( connected.begin(), connected.end(), *c ) > 0 ) {
						TR.Debug << "Movable group for node " << *n << " --> " << segment( *n ) << " : " << mg << " is present in both loops and non-loop regions. Skipping." << std::endl;
						return false;
					} else {
						TR.Debug << "Movable group for node " << *n << " --> " << segment( *n ) << " : " << mg << " is present in both loops and non-loop regions. Allowing because there is no covalent connection." << std::endl;
					}
				}
			}
		}
	}

	// check for the same MG present in different loops
	NodeSet seen;
	for ( Solution::const_iterator nodeset = solution.begin(); nodeset != solution.end(); ++nodeset ) {
		for ( NodeSet::const_iterator n = nodeset->begin(); n != nodeset->end(); ++n ) {
			std::string const & segname = segment( *n );
			core::Size const mg = perm.segment(segname).movable_group;
			if ( seen.find( mg ) != seen.end() ) {
				TR.Debug << "Movable group for node " << *n << " --> " << segname << " : " << mg << " is present in two different loops. Skipping." << std::endl;
				return false;
			}
		}
		for ( NodeSet::const_iterator n = nodeset->begin(); n != nodeset->end(); ++n ) {
			seen.insert( perm.segment( segment( *n ) ).movable_group );
		}
		TR.Debug << "Seen is now: " << seen << std::endl;
	}

	for ( Solution::const_iterator nodeset = solution.begin(); nodeset != solution.end(); ++nodeset ) {
		std::set< std::string > segments_in_set;
		for ( NodeSet::const_iterator n = nodeset->begin(); n != nodeset->end(); ++n ) {
			std::string const & segname = segment( *n );
			segments_in_set.insert( segname );
			core::Size const mg = perm.segment(segname).movable_group;
			// check to see if this mg has already been "used" i.e. referred to
			if ( fixed_mgs.find(mg) != fixed_mgs.end() ) {
				// check to see if the violating segments are connected by a bond -- this makes it OK
				std::string const & seg1 = segment( *n );
				for ( std::set< std::string >::const_iterator s2 = segments_in_set.begin(); s2 != segments_in_set.end(); ++s2 ) {
					if ( perm.segment( *s2 ).movable_group != mg ) continue;
					if ( seg1 == *s2 ) continue;
					if ( ! is_non_polymeric_connection( perm, seg1, *s2 ) ) {
						TR.Debug << "Movable group for node " << *n << " --> " << segname << " : " << mg << " has been used twice in the same loop. Skipping." << std::endl;
						return false;
					}
				}
			}
			insert_mg( fixed_mgs, mg, segname );
		}
	}

	// each set in the vector represents one loop object
	for ( core::Size i=1; i<=solution.size(); ++i ) {
		// create a set of movable groups to be included in this loop
		std::set< core::Size > mgs;
		NodeSet visited;
		std::stack< core::Size > idxs;
		for ( std::set< core::Size >::const_iterator lseg = solution[i].begin(); lseg != solution[i].end(); ++lseg ) {
			idxs.push( *lseg );
			std::string const segmentname = segment( *lseg );
			std::pair< std::set< core::Size >::iterator, bool > mg_insert_result = mgs.insert( perm.segment(segmentname).movable_group );
			if ( !mg_insert_result.second ) {
				TR.Debug << "Movable group for " << segmentname << " : " << perm.segment(segmentname).movable_group << " already existed. Skipping this solution." << std::endl;
				//return false;
			}
		}
		TR.Debug << "Movable groups in loop are: " << mgs << std::endl;
		while ( idxs.size() ) {
			core::Size const cur = idxs.top();
			idxs.pop();
			if ( visited.find( cur ) != visited.end() ) {
				continue;
			}
			visited.insert( cur );
			bool cur_in_solution = ( solution[i].find(cur) != solution[i].end() );
			if ( !cur_in_solution && ( mgs.find( perm.segment( segment(cur) ).movable_group ) != mgs.end() ) ) {
				// in this case, there must be a cutpoint between the two objects of same movable group
				TR.Debug << segment(cur) << " with mg " << perm.segment(segment(cur)).movable_group << " is covalently bound to something in the loop with the same movable group without a cutpoint in between: " << mgs << std::endl;
				return false;
			}

			// look for cutpoints within this segment
			bool has_cutpoint = false;
			if ( perm.segment(segment(cur)).cutpoint() ) {
				has_cutpoint = true;
			}
			// only continue traversing past this segment if it doesn't contain a cutpoint
			if ( has_cutpoint ) {
				continue;
			}

			//look for this segment's movable group in the solution
			core::graph::Node const * node = gpeptide_.get_node(cur);
			assert( node );
			for ( core::graph::Graph::EdgeListConstIter e=node->const_edge_list_begin(); e != node->const_edge_list_end(); ++e ) {
				core::Size const other = (*e)->get_other_ind(cur);
				idxs.push( other );
			}
		}
	}
	return true;
}

/// @brief checks to see whether sets within solutions can be combined.
void
FoldGraph::add_combined_solutions( utility::vector1< Solution > & solutions ) const
{
	utility::vector1< Solution > new_solutions;

	for ( core::Size s=1; s<=solutions.size(); ++s ) {
		utility::vector1< std::set< core::Size > > const & solution = solutions[s];
		for ( core::Size i=1; i<solution.size(); ++i ) {
			for ( core::Size j=i+1; j<=solution.size(); ++j ) {
				NodeSet merged = solution[i];
				core::Size overlap_count = 0;
				for ( NodeSet::const_iterator n = solution[j].begin(), end=solution[j].end(); n != end; ++n ) {
					core::graph::Node const * nptr = gpeptide_.get_node(*n);
					assert( nptr );
					for ( core::graph::Graph::EdgeListConstIter e=nptr->const_edge_list_begin(); e != nptr->const_edge_list_end(); ++e ) {
						core::Size const othernode = (*e)->get_other_ind(*n);
						if ( solution[i].find( othernode ) != solution[i].end() ) {
							++overlap_count;
						}
						merged.insert( *n );
					}
				}
				if ( overlap_count == 1 ) {
					Solution new_solution;
					new_solution.push_back( merged );
					for ( core::Size ii=1; ii<=solution.size(); ++ii ) {
						if ( ii == i || ii == j ) {
							continue;
						}
						new_solution.push_back( solution[ii] );
					}
					assert( new_solution.size() + 1 == solution.size() );
					new_solutions.push_back( new_solution );
				}
			}
		}
	}
	TR.Debug << " Added " << new_solutions.size() << " combined solutions." << std::endl;
	for ( core::Size i=1; i<=new_solutions.size(); ++i ) {
		solutions.push_back( new_solutions[i] );
	}
}

std::ostream & operator<<( std::ostream & os, NamedSolution const & s )
{
	for ( NamedSolution::const_iterator nset=s.begin(); nset!=s.end(); ++nset ) {
		os << *nset;
	}
	return os;
}

std::map< core::Size, std::list< core::Size > >
map_mgs( StructureData const & sd, FoldGraph const & fg )
{
	std::map< core::Size, std::list< core::Size > > mg_map;
	for ( StringList::const_iterator c=sd.segments_begin(); c!=sd.segments_end(); ++c ) {
		insert_mg( mg_map, sd.segment( *c ).movable_group, fg.nodeidx( *c ) );
	}
	return mg_map;
}

std::map< core::Size, core::Size >
map_nodes( StructureData const & sd, FoldGraph const & fg )
{
	std::map< core::Size, core::Size > node_to_mg;
	for ( StringList::const_iterator c=sd.segments_begin(); c!=sd.segments_end(); ++c ) {
		debug_assert( node_to_mg.find( fg.nodeidx( *c ) ) == node_to_mg.end() );
		node_to_mg[ fg.nodeidx( *c ) ] = sd.segment( *c ).movable_group;
	}
	return node_to_mg;
}

class SolutionSorter {
public:
	SolutionSorter( StructureData const & sd, FoldGraph const & fg ):
		mg_map_( map_mgs( sd, fg ) ), node_to_mg_( map_nodes( sd, fg ) )
	{}

	// we want the most nodes with shared MGS in the fixed set
	bool operator()( Solution const & s1, Solution const & s2 ) const
	{
		TR << "Comparing " << s1 << " to " << s2 << std::endl;
		std::set< core::Size > const nodes_in_a_loop1 = nodes_in_solution( s1 );
		std::set< core::Size > const nodes_in_a_loop2 = nodes_in_solution( s2 );
		std::set< core::Size > const fixed1 = fixed_nodes( nodes_in_a_loop1 );
		std::set< core::Size > const fixed2 = fixed_nodes( nodes_in_a_loop2 );
		return shared_mgs( fixed1 ) > shared_mgs( fixed2 );
	}

private:
	std::set< core::Size > nodes_in_solution( Solution const & s ) const
	{
		std::set< core::Size > in_a_loop;
		for ( Solution::const_iterator ns=s.begin(); ns!=s.end(); ++ns ) {
			for ( NodeSet::const_iterator n=ns->begin(); n!=ns->end(); ++n ) {
				in_a_loop.insert( *n );
			}
		}
		return in_a_loop;
	}

	std::set< core::Size > fixed_nodes( std::set< core::Size > const & nodes_in_solution ) const
	{
		std::set< core::Size > fixed;
		for ( std::map< core::Size, core::Size >::const_iterator n_mg=node_to_mg_.begin(); n_mg!=node_to_mg_.end(); ++n_mg ) {
			if ( nodes_in_solution.find( n_mg->first ) == nodes_in_solution.end() ) {
				fixed.insert( n_mg->first );
			}
		}
		return fixed;
	}

	core::Size shared_mgs( std::set< core::Size > const & fixed ) const
	{
		std::set< core::Size > seen;
		core::Size shared_count = 0;
		for ( std::set< core::Size >::const_iterator fnode=fixed.begin(); fnode!=fixed.end(); ++fnode ) {
			core::Size const mg = node_to_mg_.find( *fnode )->second;
			if ( seen.find( mg ) != seen.end() ) {
				++shared_count;
			} else {
				seen.insert( mg );
			}
		}
		TR << "Found " << shared_count << " shared mgs in fixed regions" << std::endl;
		return shared_count;
	}

	std::map< core::Size, std::list< core::Size > > mg_map_;
	std::map< core::Size, core::Size > node_to_mg_;
};

/// @brief sorts solutions such that the ones maximizing segments w/ same MG are fixed.
void
FoldGraph::sort_solutions(
	utility::vector1< Solution > & solutions,
	StructureData const & sd ) const
{
	SolutionSorter sorter( sd, *this );
	std::sort( solutions.begin(), solutions.end(), sorter );
}

/// @brief convert a solution to a named solution
NamedSolution
FoldGraph::named_solution( Solution const & solution ) const
{
	NamedSolution nsolution;
	for ( Solution::const_iterator s=solution.begin(), ends=solution.end(); s!=ends; ++s ) {
		NamedNodeSet nset;
		for ( NodeSet::const_iterator n=s->begin(), endn=s->end(); n!=endn; ++n ) {
			nset.insert( segment( *n ) );
		}
		nsolution.push_back( nset );
		debug_assert( nset.size() == s->size() );
	}
	debug_assert( nsolution.size() == solution.size() );
	return nsolution;
}

void
FoldGraph::add_non_polymeric_connections(
	utility::vector1< Solution > & solutions,
	StructureData const & perm ) const
{
	for ( utility::vector1< Solution >::iterator s = solutions.begin(); s != solutions.end(); ++s ) {
		for ( Solution::iterator sol = s->begin(); sol != s->end(); ++sol ) {
			NodeSet connected_nodes;
			for ( NodeSet::const_iterator n = sol->begin(); n != sol->end(); ++n ) {
				for ( utility::vector1< BondInfo >::const_iterator bi = perm.covalent_bonds_begin();
						bi != perm.covalent_bonds_end(); ++bi ) {
					core::Size othernode = 0;
					if ( segment( *n ) == bi->seg1 ) {
						othernode = nodeidx( bi->seg2 );
					} else if ( segment( *n ) == bi->seg2 ) {
						othernode = nodeidx( bi->seg1 );
					}
					if ( othernode ) connected_nodes.insert( othernode );
				}
			}
			if ( !connected_nodes.empty() ) {
				for ( NodeSet::const_iterator n = connected_nodes.begin(); n != connected_nodes.end(); ++n ) {
					sol->insert( *n );
				}
			} // if !connected_nodes.empty()
		} // for nodeset in s
	} // for s in solution
}

/// @brief generates a solution of movable segments to be used in remodel, based on the foldgraph
Solution
FoldGraph::compute_best_solution(
	StructureData const & perm,
	utility::vector1< std::string > const & staple_loops,
	utility::vector1< std::string > const & cut_loops ) const
{
	TR.Debug << "staple loops= " << staple_loops << " cut_loops= " << cut_loops << std::endl;

	std::set< core::Size > cut_loop_nodes;
	for ( core::Size i=1; i<=cut_loops.size(); ++i ) {
		assert( cut_loop_nodes.find(nodeidx(cut_loops[i])) == cut_loop_nodes.end() );
		cut_loop_nodes.insert( nodeidx( cut_loops[i] ) );
	}

	utility::vector1< Solution > solutions;
	for ( core::Size i=1; i<=staple_loops.size(); ++i ) {
		TR.Debug << "Running for " << staple_loops[i] << std::endl;
		core::Size nruns = solutions.size();
		bool first_run = false;
		if ( i == 1 ) {
			nruns = 1;
			first_run = true;
		}
		for ( core::Size sol=1; sol<=nruns; ++sol ) {
			// do a dfs of the peptide graph
			NodeSet visited;
			if ( !first_run ) {
				for ( core::Size sidx=1; sidx<=solutions[sol].size(); ++sidx ) {
					for ( NodeSet::const_iterator n=solutions[sol][sidx].begin(); n != solutions[sol][sidx].end(); ++n ) {
						visited.insert( *n );
					}
				}
			}
			Solution dfs_solutions;
			NodeSet const visited_copy = visited;
			TR.Debug << "Input visited " << visited << std::endl;
			create_loops_dfs( dfs_solutions, visited, nodeidx(staple_loops[i]), cut_loop_nodes, perm );
			TR.Debug << "SUCCESS: " << dfs_solutions.size() << std::endl;
			for ( core::Size j=1; j<=dfs_solutions.size(); ++j ) {
				NodeSet new_visited;
				for ( NodeSet::const_iterator v=dfs_solutions[j].begin(); v != dfs_solutions[j].end(); ++v ) {
					if ( visited_copy.find(*v) == visited_copy.end() ) {
						new_visited.insert( *v );
					}
				}
				TR.Debug << "New visited: " << new_visited << std::endl;
				if ( first_run ) {
					Solution tmpset;
					tmpset.push_back( new_visited );
					solutions.push_back( tmpset );
				} else {
					if ( j == 1 ) {
						solutions[sol].push_back( new_visited );
					} else {
						Solution tmpset;
						for ( core::Size t=1; t<solutions[sol].size(); ++t ) {
							tmpset.push_back( solutions[sol][t] );
						}
						tmpset.push_back( new_visited );
						solutions.push_back( tmpset );
					}
				}
			}
		}
	}
	// size of solutions must be >= 1, so add an empty list if solutions.size() == 0
	if ( !staple_loops.size() ) {
		solutions.push_back( Solution() );
	}
	for ( core::Size i=1; i<=solutions.size(); ++i ) {
		NodeSet visited;
		for ( core::Size j=1, endj=solutions[i].size(); j<=endj; ++j ) {
			TR.Debug << solutions[i][j] << " ";
			for ( NodeSet::const_iterator g=solutions[i][j].begin(); g != solutions[i][j].end(); ++g ) {
				visited.insert( *g );
			}
		}
		TR.Debug << std::endl;
		TR.Debug << "Adding unvisited loop nodes.  Visited nodes are " << visited << std::endl;
		for ( NodeSet::const_iterator c=cut_loop_nodes.begin(); c != cut_loop_nodes.end(); ++c ) {
			if ( visited.find(*c) == visited.end() ) {
				NodeSet tmpset;
				tmpset.insert(*c);
				solutions[i].push_back( tmpset );
			}
		}
	}

	// combine solution sets which are connected by one and only one peptide edge
	add_combined_solutions( solutions );
	sort_solutions( solutions, perm );

	// add non-polymeric connections to solutions
	utility::vector1< Solution > const solutions_return = solutions;

	add_non_polymeric_connections( solutions, perm );

	core::Size const bestidx = select_best_solution( perm, solutions );
	if ( bestidx == 0 ) {
		TR << "No valid solutions found..." << std::endl;
		return Solution();
	}
	TR << "Chose the following loops from solution " << bestidx << ":"
		<< named_solution( solutions[bestidx] ) << std::endl;

	TR.Debug << "perm=" << perm << std::endl;
	return solutions_return[ bestidx ];
}

core::Size
FoldGraph::select_best_solution(
	StructureData const & perm,
	utility::vector1< Solution > const & solutions ) const
{
	// print candidates
	for ( core::Size i=1; i<=solutions.size(); ++i ) {
		TR.Debug << "S" << i << ": " << named_solution( solutions[i] ) << std::endl;
	}

	// choose best solution == fewest loops
	core::Size bestidx = 0;
	for ( core::Size i=1; i<=solutions.size(); ++i ) {
		if ( ! check_solution( perm, solutions[i] ) ) {
			TR.Debug << "Skipping solution " << i << " because it has a movable groups conflict." << std::endl;
			continue;
		}
		if ( !bestidx || ( solutions[i].size() < solutions[bestidx].size() ) ) {
			bestidx = i;
		}
	}
	return bestidx;
}

protocols::loops::LoopsOP
FoldGraph::create_loops(
	StructureData const & perm,
	utility::vector1< std::string > const & staple_loops,
	utility::vector1< std::string > const & bridge_loops ) const
{
	Solution solution = compute_best_solution( perm, staple_loops, bridge_loops );
	return create_loops( perm, solution, bridge_loops );
}

protocols::loops::LoopsOP
FoldGraph::create_loops(
	StructureData const & perm,
	Solution const & solution,
	utility::vector1< std::string > const & cut_loops ) const
{
	if ( !solution.size() ) {
		return protocols::loops::LoopsOP();
	}

	// now make a loops object
	protocols::loops::LoopsOP loops = protocols::loops::LoopsOP( new protocols::loops::Loops() );
	for ( core::Size j=1; j<=solution.size(); ++j ) {
		int min_res = -1;
		int max_res = -1;
		for ( NodeSet::const_iterator n=solution[j].begin(); n != solution[j].end(); ++n ) {
			std::string const & seg = segment(*n);
			int const start = perm.segment(seg).nterm_resi();
			int const stop = perm.segment(seg).cterm_resi();
			if ( min_res < 0 ) {
				min_res = start;
			} else {
				if ( start < min_res ) {
					min_res = start;
				}
			}
			if ( stop < min_res ) {
				min_res = stop;
			}
			if ( max_res < 0 ) {
				max_res = start;
			} else {
				if ( start > max_res ) {
					max_res = start;
				}
			}
			if ( stop > max_res ) {
				max_res = stop;
			}
		}
		assert( min_res > 0 );
		assert( max_res > 0 );
		assert( min_res <= max_res );
		// find cutpoint
		int cut = 0;
		for ( core::Size cuti=1, endcut=cut_loops.size(); cuti<=endcut; ++cuti ) {
			int const intra_loop_cut = perm.get_data_int( cut_loops[cuti], "cut_resi" );
			int const cut_to_check = perm.segment(cut_loops[cuti]).start() + intra_loop_cut - 1;
			TR.Debug << "Checking whether " << cut_to_check << " is within " << min_res << " and " << max_res << std::endl;
			if ( ( min_res <= cut_to_check ) && ( cut_to_check <= max_res ) ) {
				// this assertion makes sure there aren't multiple cutpoints in a loop -- very important
				assert( cut == 0 );
				cut = cut_to_check;
				assert( cut == static_cast< int >( perm.segment(cut_loops[cuti]).cutpoint() ) );
			}
		}
		// if the permutation has a cut residue of 0, we set it here
		protocols::loops::Loop l( min_res, max_res, cut );
		loops->add_loop( l );
		TR << "Created " << l << std::endl;
	}
	return loops;
}

/// @brief recursive inner function that traverses the foldgraph to create loops objects
void
FoldGraph::create_loops_dfs(
	Solution & solutions,
	NodeSet & visited,
	core::Size const current_node,
	NodeSet const & cut_loop_nodes,
	StructureData const & perm ) const
{
	if ( visited.find(current_node) != visited.end() ) {
		return;
	}
	TR.Debug << "Visiting " << segment(current_node) << " " << current_node << " " << visited << std::endl;
	visited.insert(current_node);

	// if we're at a cut loop, detect it and return true
	for ( core::Size i=1; i<=cut_loop_nodes.size(); ++i ) {
		if ( cut_loop_nodes.find(current_node) != cut_loop_nodes.end() ) {
			solutions.push_back( visited );
			return;
		}
	}

	core::graph::Node const * n = gpeptide_.get_node(current_node);
	assert( n );
	bool terminal_node = true;
	for ( core::graph::Graph::EdgeListConstIter e=n->const_edge_list_begin(); e != n->const_edge_list_end(); ++e ) {
		core::Size const othernode = (*e)->get_other_ind(current_node);
		if ( visited.find( othernode ) == visited.end() ) {
			terminal_node = false;
		}
		NodeSet visited_new = visited;
		create_loops_dfs( solutions, visited_new, othernode, cut_loop_nodes, perm );
	}
	// add non-polymeric bonded things
	if ( terminal_node ) {
		solutions.push_back( visited );
	}
}

} // namespace components
} // namespace denovo_design
} // namespace protocols
