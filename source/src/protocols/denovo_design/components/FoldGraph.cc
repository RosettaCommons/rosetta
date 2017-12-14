// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/components/FoldGraph.cc
/// @brief FoldGraph - a fold-tree which takes two-way dependencies into account
/// @details
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
#include <core/pose/Pose.hh>

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

FoldGraph::FoldGraph( StructureData const & perm ):
	utility::pointer::ReferenceCount(),
	sd_( StructureDataOP( new StructureData( perm ) ) )
{
	TR.Debug << "Constructing foldgraph for " << sd().id() << " : " << sd()  << std::endl;
	g_.drop_all_edges();
	gpeptide_.drop_all_edges();
	seg2node_.clear();
	node2seg_.clear();

	build_seg2node();
	build_node2seg( seg2node_ );
	find_cutpoints();

	g_.set_num_nodes( seg2node_.size() );
	gpeptide_.set_num_nodes( seg2node_.size() );

	//must have at least one segment to move forward
	if ( sd().segments_begin() == sd().segments_end() ) {
		return;
	}

	// add edges between movable groups
	MovableGroups const movable_groups = sd().movable_groups();
	TR.Debug << "Num groups = " << movable_groups.size() << std::endl;
	debug_assert( movable_groups.size() );
	for ( unsigned long movable_group : movable_groups ) {
		SegmentNames const segments = sd().segments_in_movable_group( movable_group );
		for ( auto s1=segments.begin(); s1!=segments.end(); ++s1 ) {
			TR.Debug << "Group " << movable_group << " Segment " << *s1 << std::endl;
			auto next = s1;
			++next;
			for ( auto s2=next; s2!=segments.end(); ++s2 ) {
				debug_assert( *s1 != *s2 );
				add_edge( *s1, *s2 );
			}
		}
	}

	// look for non-canonical bonds
	for ( auto bi=sd().covalent_bonds_begin(); bi!=sd().covalent_bonds_end(); ++bi ) {
		TR << "Found non-canonical connection from " << bi->seg1 << ":" << bi->res1 << ":" << bi->atom1
			<<" to " << bi->seg2 << ":" << bi->res2 << ":" << bi->atom2 << std::endl;
		if ( bi->seg1 != bi->seg2 ) {
			add_edge( bi->seg1, bi->seg2 );
		}
	}

	// add peptide edges
	for ( auto n=sd().segments_begin(); n!=sd().segments_end(); ++n ) {
		Segment const & resis = sd().segment(*n);
		if ( !resis.has_free_lower_terminus() ) {
			TR.Debug << "Segment1 " << *n << " res " << resis.safe() << std::endl;
			TR.Debug << "Segment2 " << resis.lower_segment() << " res " << sd().segment( resis.lower_segment() ).safe() << std::endl;
			add_peptide_edge( *n, resis.lower_segment() );
		}
		if ( !resis.has_free_upper_terminus() ) {
			TR.Debug << "Segment1 " << *n << " res " << resis.safe() << std::endl;
			TR.Debug << "Segment2 " << resis.upper_segment() << " res " << sd().segment( resis.upper_segment() ).safe() << std::endl;
			add_peptide_edge( *n, resis.upper_segment() );
		}
	}
}

FoldGraph::~FoldGraph() = default;

void
FoldGraph::find_cutpoints()
{
	cutpoints_.clear();
	for ( auto s=sd().segments_begin(); s!=sd().segments_end(); ++s ) {
		core::Size cut = sd().segment( *s ).cutpoint();
		if ( cut != 0 ) cutpoints_.insert( nodeidx( *s ) );
	}
}

void
FoldGraph::build_seg2node()
{
	// populate seg2node
	for ( auto r=sd().segments_begin(); r!=sd().segments_end(); ++r ) {
		core::Size const nodenum = seg2node_.size() + 1;
		auto s2n = seg2node_.find( *r );
		if ( s2n == seg2node_.end() ) {
			s2n = seg2node_.insert( std::make_pair( *r, nodenum ) ).first;
			TR.Debug << "Added " << *r << std::endl;
		} else {
			std::stringstream msg;
			msg << "FoldGraph::build_node2seg(): Attempted to add a segment (" << *r
				<< ") and node number (" << nodenum << "), but the segment name is already present in the seg2node_ map. Map = "
				<< seg2node_ << std::endl;
			msg << "SD = " << sd() << std::endl;
			utility_exit_with_message( msg.str() );
		}
		debug_assert( s2n != seg2node_.end() );
	}
	debug_assert( seg2node_.size() == static_cast< core::Size >( std::distance( sd().segments_begin(), sd().segments_end() ) ) );
}

void
FoldGraph::build_node2seg( SegmentToNodeMap const & seg2node )
{
	// populate node2seg
	for ( auto const & s2n : seg2node ) {
		auto n2s = node2seg_.find( s2n.second );
		if ( n2s == node2seg_.end() ) {
			n2s = node2seg_.insert( std::make_pair( s2n.second, s2n.first ) ).first;
		} else {
			std::stringstream msg;
			msg << "FoldGraph::build_node2seg(): Attempted to add a segment (" << s2n.first
				<< ") and node number (" << s2n.second << "), but the node number is already present in the node2seg_ map. SegmentToNodeMap = "
				<< seg2node_ << std::endl;
			utility_exit_with_message( msg.str() );
		}
		debug_assert( n2s != node2seg_.end() );
	}
	debug_assert( seg2node.size() == node2seg_.size() );
}

/// @brief given a segment name, returns the associated graph node index
core::Size
FoldGraph::nodeidx( std::string const & segment ) const
{
	auto it = seg2node_.find(segment);
	if ( it == seg2node_.end() ) {
		TR.Error << "Component " << segment << " not found in map: " << seg2node_ << std::endl;
		debug_assert( it != seg2node_.end() );
	}
	return it->second;
}

/// @brief given a node index, returns the associated segment
std::string const &
FoldGraph::segment( core::Size const nodeidx ) const
{
	auto it = node2seg_.find(nodeidx);
	if ( it == node2seg_.end() ) {
		TR.Error << "Node idx " << nodeidx << " not found in map: " << seg2node_ << std::endl;
		debug_assert( it != node2seg_.end() );
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

StructureData const &
FoldGraph::sd() const
{
	return *sd_;
}

/// @brief gives a fold tree based on the segments in the SD
core::kinematics::FoldTree
FoldGraph::fold_tree( SegmentNames const & root_segments ) const
{
	debug_assert( root_segments.size() >= 1 );
	core::kinematics::FoldTree ft;
	NodeSet visited;
	std::stack< core::Size > to_search;

	// add to fold tree based on peptide graph
	Segment const & root_resis = sd().segment( *root_segments.begin() );
	fold_tree_rec( ft, visited, to_search, *root_segments.begin(), root_resis.safe(), 0, true );
	for ( auto s=++root_segments.begin(); s!=root_segments.end(); ++s ) {
		if ( visited.find( nodeidx( *s ) ) != visited.end() ) continue;
		Segment const & resis = sd().segment( *s );
		core::kinematics::Edge e( root_resis.safe(), resis.safe(), ft.num_jump()+1 );
		TR.Debug << " adding unconnected segment from peptide graph rooted at " << *s << " " << e << std::endl;
		ft.add_edge( e );
		fold_tree_rec( ft, visited, to_search, *s, e.stop(), 0, true );
	}

	// add to fold tree based on interaction graph -- look for reachable nodes that are not yet visited.
	while ( to_search.size() > 0 ) {
		core::Size const node = to_search.top();
		to_search.pop();
		utility::graph::Node const * n = g_.get_node(node);
		debug_assert( n );
		for ( utility::graph::Graph::EdgeListConstIter e=n->const_edge_list_begin(); e != n->const_edge_list_end(); ++e ) {
			core::Size const othernode = (*e)->get_other_ind(node);
			if ( visited.find(othernode) != visited.end() ) {
				TR.Debug << "ALready visited " << segment(node) << "__" << segment(othernode) << std::endl;
				continue;
			}
			core::kinematics::Edge edge( sd().segment(segment(node)).safe(), sd().segment(segment(othernode)).safe(), ft.num_jump()+1 );
			auto bi = sd().find_non_polymer_bond( segment(node), segment(othernode) );
			if ( bi == sd().covalent_bonds_end() ) {
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
				edge.start() = sd().pose_residue( this_seg, this_res );
				edge.stop() = sd().pose_residue( other_seg, other_res );
				edge.start_atom() = this_atom;
				edge.stop_atom() = other_atom;
				core::kinematics::FoldTree const & ft_const = ft;
				utility::vector1< core::kinematics::Edge > newedges;
				utility::vector1< core::kinematics::Edge > deledges;
				for ( auto const & e : ft_const ) {
					if ( e.label() > 0 ) continue;
					if ( ( edge.start() >= e.stop() ) && ( edge.start() <= e.start() ) ) {
						TR << "Removing edge: " << e << " to insert " << edge << std::endl;
						newedges.push_back( core::kinematics::Edge( e.stop(), edge.start(), core::kinematics::Edge::PEPTIDE ) );
						newedges.push_back( core::kinematics::Edge( edge.start(), e.start(), core::kinematics::Edge::PEPTIDE ) );
						deledges.push_back( e );
					} else if ( ( edge.start() >= e.start() ) && ( edge.start() <= e.stop() ) ) {
						TR << "2Adding edge: " << e << " to insert " << edge << std::endl;
						newedges.push_back( core::kinematics::Edge( e.start(), edge.start(), core::kinematics::Edge::PEPTIDE ) );
						newedges.push_back( core::kinematics::Edge( edge.start(), e.stop(), core::kinematics::Edge::PEPTIDE ) );
						deledges.push_back( e );
					}
					if ( ( edge.stop() >= e.stop() ) && ( edge.stop() <= e.start() ) ) {
						TR << "Removing edge: " << e << " to insert " << edge << std::endl;
						newedges.push_back( core::kinematics::Edge( e.stop(), edge.stop(), core::kinematics::Edge::PEPTIDE ) );
						newedges.push_back( core::kinematics::Edge( edge.stop(), e.start(), core::kinematics::Edge::PEPTIDE ) );
						deledges.push_back( e );
					} else if ( ( edge.stop() >= e.start() ) && ( edge.stop() <= e.stop() ) ) {
						TR << "2Adding edge: " << e << " to insert " << edge << std::endl;
						newedges.push_back( core::kinematics::Edge( e.start(), edge.stop(), core::kinematics::Edge::PEPTIDE ) );
						newedges.push_back( core::kinematics::Edge( edge.stop(), e.stop(), core::kinematics::Edge::PEPTIDE ) );
						deledges.push_back( e );
					}
				}
				for ( utility::vector1< core::kinematics::Edge >::const_iterator e=newedges.begin(); e!=newedges.end(); ++e ) {
					ft.add_edge( *e );
				}
				for ( utility::vector1< core::kinematics::Edge >::const_iterator e=deledges.begin(); e!=deledges.end(); ++e ) {
					ft.delete_edge( *e );
				}
			}
			ft.add_edge( edge );
			fold_tree_rec( ft, visited, to_search, segment(othernode), edge.stop(), 0, true );
		}
	}

	// add unconnected foldtree elements directly to the start node
	for ( auto const & s2n : seg2node_ ) {
		if ( visited.find( s2n.second ) == visited.end() ) {
			std::string const & unconnected_seg = segment(s2n.second);
			Segment const & resis = sd().segment(unconnected_seg);
			core::kinematics::Edge e( root_resis.safe(), resis.safe(), ft.num_jump()+1 );
			TR.Debug << " unconnected_comp = " << unconnected_seg << " " << e << std::endl;
			ft.add_edge( e );
			fold_tree_rec( ft, visited, to_search, unconnected_seg, e.stop(), 0, true );
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
	std::string const & segment_name,
	core::Size const parent_residue,
	int const parent_direction,
	bool const polymer_only ) const
{
	core::Size const nodenum = nodeidx( segment_name );
	if ( visited.find(nodenum) != visited.end() ) {
		return;
	}

	TR.Debug << "Starting traversal at " << segment_name <<  " node " << nodenum << " visited " << visited << " parent=" << parent_residue << " direction=" << parent_direction << std::endl;

	Segment const & resis = sd().segment(segment_name);
	debug_assert( resis.lower() <= parent_residue );
	debug_assert( parent_residue <= resis.upper() );

	if ( parent_direction == 0 ) {
		//previous residue is jump
		if ( parent_residue != resis.safe() ) {
			TR.Debug << "Adding " << parent_residue << " " << resis.safe() << " -1 for " << segment_name << std::endl;
			ft.add_edge( parent_residue, resis.safe(), core::kinematics::Edge::PEPTIDE );
		}
		if ( parent_residue > resis.safe() ) {
			TR.Debug << "Adding " << parent_residue << " " << resis.upper() << " -1 for " << segment_name << std::endl;
			ft.add_edge( parent_residue, resis.upper(), core::kinematics::Edge::PEPTIDE );
			TR.Debug << "Adding " << resis.safe() << " " << resis.lower() << " -1 for " << segment_name << std::endl;
			ft.add_edge( resis.safe(), resis.lower(), core::kinematics::Edge::PEPTIDE );
		} else {
			TR.Debug << "Adding " << parent_residue << " " << resis.lower() << " -1 for " << segment_name << std::endl;
			ft.add_edge( parent_residue, resis.lower(), core::kinematics::Edge::PEPTIDE );
			TR.Debug << "Adding " << resis.safe() << " " << resis.upper() << " -1 for " << segment_name << std::endl;
			ft.add_edge( resis.safe(), resis.upper(), core::kinematics::Edge::PEPTIDE );
		}
	} else if ( parent_direction > 0 ) {
		// previous is a peptide edge going up
		if ( resis.safe() != resis.lower() ) {
			TR.Debug << "Adding " << resis.lower() << " " << resis.safe() << " -1 for " << segment_name << std::endl;
			ft.add_edge( resis.lower(), resis.safe(), core::kinematics::Edge::PEPTIDE );
		}
		if ( resis.safe() != resis.upper() ) {
			TR.Debug << "Adding " << resis.safe() << " " << resis.upper() << " -1 for " << segment_name << std::endl;
			ft.add_edge( resis.safe(), resis.upper(), core::kinematics::Edge::PEPTIDE );
		}
	} else {
		// previous is a peptide edge going down
		if ( resis.safe() != resis.upper() ) {
			TR.Debug << "Adding " << resis.upper() << " " << resis.safe() << " -1 for " << segment_name << std::endl;
			ft.add_edge( resis.upper(), resis.safe(), core::kinematics::Edge::PEPTIDE );
		}
		if ( resis.safe() != resis.lower() ) {
			TR.Debug << "Adding " << resis.safe() << " " << resis.lower() << " -1 for " << segment_name << std::endl;
			ft.add_edge( resis.safe(), resis.lower(), core::kinematics::Edge::PEPTIDE );
		}
	} // if direction

	visited.insert( nodenum );
	node_stack.push( nodenum );

	// traverse connected edges
	utility::graph::Node const * n;
	if ( polymer_only ) {
		n = gpeptide_.get_node(nodenum);
	} else {
		n = g_.get_node(nodenum);
	}
	debug_assert( n );
	for ( utility::graph::Graph::EdgeListConstIter e=n->const_edge_list_begin(); e!=n->const_edge_list_end(); ++e ) {
		core::Size const othernode = (*e)->get_other_ind(nodenum);
		if ( visited.find(othernode) != visited.end() ) {
			TR.Debug << "Skipping " << othernode << " as it has been visited." << std::endl;
			continue;
		}

		// add connecting edge
		std::string const & otherseg = segment(othernode);
		Segment const & other_resis = sd().segment( otherseg );
		if ( resis.upper() <= other_resis.lower() ) {
			int direction = 1;
			core::kinematics::Edge e( resis.upper(), other_resis.lower(), core::kinematics::Edge::PEPTIDE );
			if ( !polymer_only ) {
				e = core::kinematics::Edge( resis.safe(), other_resis.safe(), ft.num_jump()+1 );
				direction = 0;
			}
			TR.Debug << "Adding " << e << " for " << segment_name << "-->" << otherseg << std::endl;
			ft.add_edge( e );
			fold_tree_rec( ft, visited, node_stack, otherseg, e.stop(), direction, true );
		} else if ( resis.lower() >= other_resis.upper() ) {
			int direction = -1;
			core::kinematics::Edge e( resis.lower(), other_resis.upper(), core::kinematics::Edge::PEPTIDE );
			if ( !polymer_only ) {
				e = core::kinematics::Edge( resis.safe(), other_resis.safe(), ft.num_jump()+1 );
				direction = 0;
			}
			TR.Debug << "Adding " << e << " for " << segment_name << "-->" << otherseg << std::endl;
			ft.add_edge( e );
			fold_tree_rec( ft, visited, node_stack, otherseg, e.stop(), direction, true );
		} else {
			TR.Fatal << "Something is wrong. THere are probably two segments that overlap in residue numbering. Nterm=" << resis.lower() << " Cterm=" << resis.upper() << " safe=" << resis.safe() << " otherNterm=" << other_resis.lower() << " otherCterm=" << other_resis.upper() << " safe=" << other_resis.safe() << std::endl;
			utility_exit_with_message("Likely residue overlap in FoldGraph.");
		}
	}
}

template< class T >
void
insert_mg(
	std::map< core::Size, std::list< T > > & fixed_mgs,
	core::Size const mg,
	T const & segname )
{
	auto fixed_mg = fixed_mgs.find( mg );
	if ( fixed_mg == fixed_mgs.end() ) {
		fixed_mg = fixed_mgs.insert( std::make_pair( mg, std::list< T >() ) ).first;
	}
	debug_assert( fixed_mg != fixed_mgs.end() );
	fixed_mg->second.push_back( segname );
}

/// @brief checks a solution to ensure that covalently bound segments are not found both inside and outside of loops
bool
FoldGraph::check_solution( Solution const & solution ) const
{
	SegmentNameSet loop_segments;
	for ( auto const & s : solution ) {
		for ( unsigned long lseg : s ) {
			SegmentName const & segmentname = segment( lseg );
			std::pair< SegmentNameSet::iterator, bool > const ins_result = loop_segments.insert( segmentname );
			if ( !ins_result.second ) {
				TR << "Segment " << segmentname << " is present in more than one loop... this should not be possible. Skipping." << std::endl;
				return false;
			}
		}
	}

	// matches a fixed mg with a list of segments
	MovableGroupToSegmentNameListMap fixed_mgs;

	// ensure that the same MG isn't present both inside and outside of loops
	for ( auto c=sd().segments_begin(); c!=sd().segments_end(); ++c ) {
		MovableGroup const mg = sd().segment(*c).movable_group();
		if ( loop_segments.find(*c) == loop_segments.end() ) {
			insert_mg( fixed_mgs, mg, *c );
		}
	}

	for ( auto const & nodeset : solution ) {
		for ( unsigned long n : nodeset ) {
			std::string const & segname = segment( n );
			core::Size const mg = sd().segment( segname ).movable_group();
			MovableGroupToSegmentNameListMap::const_iterator fixed_mg = fixed_mgs.find( mg );
			if ( fixed_mg != fixed_mgs.end() ) {
				SegmentNameList const connected = sd().connected_segments( segname, false );
				for ( auto const & c : fixed_mg->second ) {
					// check to see if fixed segment is connected to anything connected to *n
					if ( std::count( connected.begin(), connected.end(), c ) > 0 ) {
						TR.Debug << "Movable group for node " << n << " --> " << segment( n ) << " : " << mg << " is present in both loops and non-loop regions. Skipping." << std::endl;
						return false;
					} else {
						TR.Debug << "Movable group for node " << n << " --> " << segment( n ) << " : " << mg << " is present in both loops and non-loop regions. Allowing because there is no covalent connection." << std::endl;
					}
				}
			}
		}
	}

	// check for the same MG present in different loops
	NodeSet seen;
	for ( auto const & nodeset : solution ) {
		for ( unsigned long n : nodeset ) {
			SegmentName const & segname = segment( n );
			MovableGroup const mg = sd().segment(segname).movable_group();
			if ( seen.find( mg ) != seen.end() ) {
				TR.Debug << "Movable group for node " << n << " --> " << segname << " : " << mg << " is present in two different loops. Skipping." << std::endl;
				return false;
			}
		}
		for ( unsigned long n : nodeset ) {
			seen.insert( sd().segment( segment( n ) ).movable_group() );
		}
		TR.Debug << "Seen is now: " << seen << std::endl;
	}

	for ( auto const & nodeset : solution ) {
		SegmentNames segments_in_set;
		for ( unsigned long n : nodeset ) {
			segments_in_set.push_back( segment( n ) );
		}

		for ( unsigned long n : nodeset ) {
			SegmentName const & segname = segment( n );
			MovableGroup const mg = sd().segment(segname).movable_group();
			// check to see if this mg has already been "used" i.e. referred to
			if ( fixed_mgs.find(mg) == fixed_mgs.end() ) {
				insert_mg( fixed_mgs, mg, segname );
				continue;
			}

			// after here, we have a possible movable group violation
			// check to see if the violating segments are connected by a bond -- this makes it OK
			SegmentNameList const connected_segment_list = sd().connected_segments( segname, false );
			SegmentNameSet const connected( connected_segment_list.begin(), connected_segment_list.end() );

			for ( SegmentNames::const_iterator s2=segments_in_set.begin(); s2!=segments_in_set.end(); ++s2 ) {
				// different MG - skip
				if ( sd().segment( *s2 ).movable_group() != mg ) continue;
				// same segment - skip
				if ( segname == *s2 ) continue;
				// connected to current segment - skip
				if ( connected.find( *s2 ) != connected.end() ) continue;
				// connected to current segment via non-polymer bond - skip
				if ( sd().non_polymer_bond_exists( segname, *s2 ) ) continue;

				TR.Debug << "Movable group for node " << n << " --> " << segname << " : " << mg << " has been used twice in the same loop. Skipping." << std::endl;
				return false;
			}
			insert_mg( fixed_mgs, mg, segname );
		}
	}

	// each set in the vector represents one loop object
	for ( auto const & s : solution ) {
		// create a set of movable groups to be included in this loop
		std::set< core::Size > mgs;
		NodeSet visited;
		std::stack< core::Size > idxs;
		for ( unsigned long lseg : s ) {
			idxs.push( lseg );
			std::string const & segmentname = segment( lseg );
			std::pair< NodeSet::iterator, bool > mg_insert_result = mgs.insert( sd().segment(segmentname).movable_group() );
			if ( !mg_insert_result.second ) {
				TR.Debug << "Movable group for " << segmentname << " : " << sd().segment(segmentname).movable_group() << " already existed. Skipping this solution." << std::endl;
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
			bool cur_in_solution = ( s.find(cur) != s.end() );
			if ( !cur_in_solution && ( mgs.find( sd().segment( segment(cur) ).movable_group() ) != mgs.end() ) ) {
				// in this case, there must be a cutpoint between the two objects of same movable group
				TR.Debug << segment(cur) << " with mg " << sd().segment(segment(cur)).movable_group() << " is covalently bound to something in the loop with the same movable group without a cutpoint in between: " << mgs << std::endl;
				return false;
			}

			// only continue traversing past this segment if it doesn't contain a cutpoint
			if ( sd().segment( segment(cur) ).cutpoint() != 0 ) {
				continue;
			}

			//look for this segment's movable group in the solution
			utility::graph::Node const * node = gpeptide_.get_node(cur);
			debug_assert( node );
			for ( utility::graph::Graph::EdgeListConstIter e=node->const_edge_list_begin(); e!=node->const_edge_list_end(); ++e ) {
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
				for ( auto n = solution[j].begin(), end=solution[j].end(); n != end; ++n ) {
					utility::graph::Node const * nptr = gpeptide_.get_node(*n);
					debug_assert( nptr );
					for ( utility::graph::Graph::EdgeListConstIter e=nptr->const_edge_list_begin(); e != nptr->const_edge_list_end(); ++e ) {
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
					debug_assert( new_solution.size() + 1 == solution.size() );
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

std::map< core::Size, std::list< core::Size > >
map_mgs( StructureData const & sd, FoldGraph const & fg )
{
	std::map< core::Size, std::list< core::Size > > mg_map;
	for ( auto c=sd.segments_begin(); c!=sd.segments_end(); ++c ) {
		insert_mg( mg_map, sd.segment( *c ).movable_group(), fg.nodeidx( *c ) );
	}
	return mg_map;
}

std::map< core::Size, core::Size >
map_nodes( StructureData const & sd, FoldGraph const & fg )
{
	std::map< core::Size, core::Size > node_to_mg;
	for ( auto c=sd.segments_begin(); c!=sd.segments_end(); ++c ) {
		debug_assert( node_to_mg.find( fg.nodeidx( *c ) ) == node_to_mg.end() );
		node_to_mg[ fg.nodeidx( *c ) ] = sd.segment( *c ).movable_group();
	}
	return node_to_mg;
}

class SolutionSorter {
public:
	SolutionSorter( StructureData const & sd, FoldGraph const & fg ):
		mg_map_( map_mgs( sd, fg ) ), node_to_mg_( map_nodes( sd, fg ) )
	{}

	// we want the most nodes with shared MGS in the fixed set
	bool
	operator()( FoldGraph::Solution const & s1, FoldGraph::Solution const & s2 ) const
	{
		TR << "Comparing " << s1 << " to " << s2 << std::endl;
		FoldGraph::NodeSet const nodes_in_a_loop1 = nodes_in_solution( s1 );
		FoldGraph::NodeSet const nodes_in_a_loop2 = nodes_in_solution( s2 );
		FoldGraph::NodeSet const fixed1 = fixed_nodes( nodes_in_a_loop1 );
		FoldGraph::NodeSet const fixed2 = fixed_nodes( nodes_in_a_loop2 );
		return shared_mgs( fixed1 ) > shared_mgs( fixed2 );
	}

private:
	FoldGraph::NodeSet
	nodes_in_solution( FoldGraph::Solution const & s ) const
	{
		FoldGraph::NodeSet in_a_loop;
		for ( auto const & ns : s ) {
			for ( unsigned long n : ns ) {
				in_a_loop.insert( n );
			}
		}
		return in_a_loop;
	}

	FoldGraph::NodeSet
	fixed_nodes( FoldGraph::NodeSet const & nodes_in_solution ) const
	{
		FoldGraph::NodeSet fixed;
		for ( auto n_mg : node_to_mg_ ) {
			if ( nodes_in_solution.find( n_mg.first ) == nodes_in_solution.end() ) {
				fixed.insert( n_mg.first );
			}
		}
		return fixed;
	}

	core::Size
	shared_mgs( FoldGraph::NodeSet const & fixed ) const
	{
		FoldGraph::NodeSet seen;
		core::Size shared_count = 0;
		for ( unsigned long fnode : fixed ) {
			core::Size const mg = node_to_mg_.find( fnode )->second;
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
FoldGraph::sort_solutions( Solutions & solutions ) const
{
	SolutionSorter sorter( sd(), *this );
	std::sort( solutions.begin(), solutions.end(), sorter );
}

void
FoldGraph::add_non_polymeric_connections( Solutions & solutions ) const
{
	for ( auto & solution : solutions ) {
		for ( auto & sol : solution ) {
			NodeSet connected_nodes;
			for ( unsigned long n : sol ) {
				for ( auto bi = sd().covalent_bonds_begin();
						bi != sd().covalent_bonds_end(); ++bi ) {
					core::Size othernode = 0;
					if ( segment( n ) == bi->seg1 ) {
						othernode = nodeidx( bi->seg2 );
					} else if ( segment( n ) == bi->seg2 ) {
						othernode = nodeidx( bi->seg1 );
					}
					if ( othernode ) connected_nodes.insert( othernode );
				}
			}
			if ( !connected_nodes.empty() ) {
				for ( unsigned long connected_node : connected_nodes ) {
					sol.insert( connected_node );
				}
			} // if !connected_nodes.empty()
		} // for nodeset in s
	} // for s in solution
}

/// @brief generates a solution of movable segments to be used in remodel, based on the foldgraph
FoldGraph::Solution
FoldGraph::compute_best_solution( SegmentNames const & staple_loops ) const
{
	TR.Debug << "staple loops= " << staple_loops << " cutpoints= " << cutpoints_ << std::endl;

	Solutions solutions;
	for ( auto const & staple_loop : staple_loops ) {
		TR.Debug << "Running for " << staple_loop << std::endl;
		bool const first_run = ( solutions.size() == 0 );

		core::Size nruns;
		if ( first_run ) {
			nruns = 1;
		} else {
			nruns = solutions.size();
		}

		for ( core::Size sol=1; sol<=nruns; ++sol ) {
			// do a dfs of the peptide graph
			NodeSet visited;
			if ( !first_run ) {
				for ( Solution::const_iterator s=solutions[sol].begin(); s!=solutions[sol].end(); ++s ) {
					for ( unsigned long n : *s ) {
						visited.insert( n );
					}
				}
			}
			Solution dfs_solution;
			NodeSet const visited_copy = visited;
			TR.Debug << "Input visited " << visited << std::endl;
			create_loops_dfs( dfs_solution, visited, nodeidx(staple_loop) );
			TR.Debug << "SUCCESS: " << dfs_solution.size() << std::endl;
			for ( Solution::const_iterator dfs_sol=dfs_solution.begin(); dfs_sol!=dfs_solution.end(); ++dfs_sol ) {
				NodeSet new_visited;
				for ( unsigned long v : *dfs_sol ) {
					if ( visited_copy.find(v) == visited_copy.end() ) {
						new_visited.insert( v );
					}
				}
				TR.Debug << "New visited: " << new_visited << std::endl;
				if ( first_run ) {
					Solution tmpset;
					tmpset.push_back( new_visited );
					solutions.push_back( tmpset );
				} else {
					if ( dfs_sol == dfs_solution.begin() ) {
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
	for ( auto & solution : solutions ) {
		NodeSet visited;
		for ( Solution::const_iterator nset=solution.begin(); nset!=solution.end(); ++nset ) {
			TR.Debug << *nset << " ";
			for ( unsigned long g : *nset ) {
				visited.insert( g );
			}
		}
		TR.Debug << std::endl;
		TR.Debug << "Adding unvisited loop nodes.  Visited nodes are " << visited << std::endl;
		for ( unsigned long cutpoint : cutpoints_ ) {
			if ( visited.find(cutpoint) == visited.end() ) {
				NodeSet tmpset;
				tmpset.insert(cutpoint);
				solution.push_back( tmpset );
			}
		}
	}

	// combine solution sets which are connected by one and only one peptide edge
	add_combined_solutions( solutions );
	sort_solutions( solutions );

	// add non-polymeric connections to solutions
	Solutions const solutions_return = solutions;

	add_non_polymeric_connections( solutions );

	Solution const & best_solution = select_best_solution( solutions );

	TR << "Chose the following loops from best solution " << NamedSolution( *this, best_solution ) << std::endl;

	TR.Debug << "perm=" << sd() << std::endl;
	return best_solution;
}

FoldGraph::Solution const &
FoldGraph::select_best_solution( Solutions const & solutions ) const
{
	// print candidates
	TR.Debug << NamedSolutions( *this, solutions ) << std::endl;

	// choose best solution == fewest loops
	auto best = solutions.end();
	core::Size i = 1;
	for ( auto s=solutions.begin(); s!=solutions.end(); ++s, ++i ) {
		if ( ! check_solution( *s ) ) {
			TR.Debug << "Skipping solution " << i << " because it has a movable groups conflict." << std::endl;
			continue;
		}
		if ( ( best == solutions.end() ) || ( s->size() < best->size() ) ) {
			best = s;
		}
	}

	if ( best == solutions.end() ) {
		std::stringstream msg;
		msg << "FoldGraph::select_best_solution(): No valid solutions could be selected." << std::endl;
		msg << "Solutions found: " << NamedSolutions( *this, solutions ) << std::endl;
		msg << "SD: " << sd() << std::endl;
		utility_exit_with_message( msg.str() );
	}

	return *best;
}

protocols::loops::LoopsOP
FoldGraph::create_loops( SegmentNames const & staple_loops ) const
{
	Solution solution = compute_best_solution( staple_loops );
	return create_loops( solution );
}

protocols::loops::LoopsOP
FoldGraph::create_loops( Solution const & solution ) const
{
	if ( !solution.size() ) {
		return protocols::loops::LoopsOP();
	}

	// now make a loops object
	protocols::loops::LoopsOP loops = protocols::loops::LoopsOP( new protocols::loops::Loops() );
	for ( auto const & nset : solution ) {
		int min_res = -1;
		int max_res = -1;
		for ( unsigned long n : nset ) {
			std::string const & seg = segment(n);
			int const start = sd().segment(seg).lower();
			int const stop = sd().segment(seg).upper();
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
		debug_assert( min_res > 0 );
		debug_assert( max_res > 0 );
		debug_assert( min_res <= max_res );

		// find cutpoint
		int cut = find_cutpoint_in_range( min_res, max_res );

		protocols::loops::Loop l( min_res, max_res, cut );
		loops->add_loop( l );
		TR << "Created " << l << std::endl;
	}
	return loops;
}

int
FoldGraph::find_cutpoint_in_range( int const min_res, int const max_res ) const
{
	int cut = 0;
	for ( unsigned long cutpoint : cutpoints_ ) {
		std::string const & segmentname = segment( cutpoint );
		int const cut_to_check = sd().segment( segmentname ).cutpoint();
		TR.Debug << "Checking whether " << cut_to_check << " is within " << min_res << " and " << max_res << std::endl;
		if ( ( min_res <= cut_to_check ) && ( cut_to_check <= max_res ) ) {
			if ( cut != 0 ) {
				std::stringstream msg;
				msg << "FoldGraph::create_loops(): More than one cutpoint found inside loop from "
					<< min_res << " to " << max_res << ". The first was " << cut << " and the current is "
					<< cut_to_check << "." << std::endl;
				msg << "SD: " << sd() << std::endl;
				utility_exit_with_message( msg.str() );
			}
			cut = cut_to_check;
			debug_assert( cut == static_cast< int >( sd().segment( segmentname ).cutpoint() ) );
		}
	}
	return cut;
}

/// @brief recursive inner function that traverses the foldgraph to create loops objects
void
FoldGraph::create_loops_dfs(
	Solution & solutions,
	NodeSet & visited,
	core::Size const current_node ) const
{
	if ( visited.find(current_node) != visited.end() ) {
		return;
	}
	TR.Debug << "Visiting " << segment(current_node) << " " << current_node << " " << visited << std::endl;
	visited.insert(current_node);

	// stop if we have reached a cutpoint
	if ( cutpoints_.find( current_node ) != cutpoints_.end() ) {
		solutions.push_back( visited );
		return;
	}

	utility::graph::Node const * n = gpeptide_.get_node(current_node);
	debug_assert( n );
	bool terminal_node = true;
	for ( utility::graph::Graph::EdgeListConstIter e=n->const_edge_list_begin(); e!=n->const_edge_list_end(); ++e ) {
		core::Size const othernode = (*e)->get_other_ind(current_node);
		if ( visited.find( othernode ) == visited.end() ) {
			terminal_node = false;
		}
		NodeSet visited_new = visited;
		create_loops_dfs( solutions, visited_new, othernode );
	}
	// add non-polymeric bonded things
	if ( terminal_node ) {
		solutions.push_back( visited );
	}
}

NamedSolutions::NamedSolutions( FoldGraph const & fg, FoldGraph::Solutions const & solutions ):
	utility::vector1< NamedSolution >()
{
	for ( auto const & solution : solutions ) {
		push_back( NamedSolution( fg, solution ) );
	}
}

std::ostream &
operator<<( std::ostream & os, NamedSolutions const & solutions )
{
	core::Size i = 1;
	for ( auto s=solutions.begin(); s!=solutions.end(); ++s, ++i ) {
		os << "S" << i << ": " << *s << std::endl;
	}
	return os;
}

NamedSolution::NamedSolution( FoldGraph const & fg, FoldGraph::Solution const & solution ):
	utility::vector1< NamedNodeSet >()
{
	for ( auto s=solution.begin(); s!=solution.end(); ++s ) {
		NamedNodeSet nset;
		for ( unsigned long n : *s ) {
			nset.insert( fg.segment( n ) );
		}
		push_back( nset );
		debug_assert( nset.size() == s->size() );
	}
}

std::ostream &
operator<<( std::ostream & os, NamedSolution const & s )
{
	for ( auto const & nset : s ) {
		os << nset;
	}
	return os;
}


} // namespace components
} // namespace denovo_design
} // namespace protocols
