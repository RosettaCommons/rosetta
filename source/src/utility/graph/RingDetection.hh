// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
/// @begin RingDetection
///
/// @brief determines all possible rings in a graph.
///
/// Hanser, Th., Jauffret, Ph., Kaufmann, G., "A New Algorithm for Exhaustive Ring Perception in
/// a Molecular Graph", Laboratoire de Modeles Informatiques Appliques a la Synthese, URA 405 du CNRS, Universite
/// Louis Pasteur, 67000 Strasbourg, France, Received May 15, 1996
/// @date 02/13/2014
///
/// This detects rings in residues!
///
/// @author
/// Steven Combs (steven.combs1@gmail.com), Ralf Mueller, Jeff Menden
///
/////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifndef INCLUDED_utility_graph_RingDetection_hh
#define INCLUDED_utility_graph_RingDetection_hh


// Unit headers
#include <platform/types.hh>
#include <boost/graph/adjacency_list.hpp>
#include <vector>
#include <utility/graph/BFS_prune.hh>


namespace utility {
namespace graph {


/// @brief basic chemical Bond
///
/// @details name, element, certain properties and parameters from .params file
///
template< class Graph >
class RingDetection: public utility::graph::null_bfs_prune_visitor {
	typedef typename boost::graph_traits<Graph>::edge_descriptor Edge;
	typedef typename boost::graph_traits<Graph>::vertex_descriptor VD;
	typedef typename boost::graph_traits<Graph>::vertex_iterator VIter;
	typedef typename boost::graph_traits<Graph>::edge_iterator EIter;
	typedef typename boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;

public:

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	//! @brief default constructor
	RingDetection() :
	paths_(),
	rings_(),
	index_to_vd_(),
	vd_to_index_()
	{}

	//! @brief constructor from a graph (either GraphWithData or ConstGraph)
	//! @param graph graph for exhaustive ring detection
	//! @note prefer using ConstGraph here, it is often many times faster than GraphWithData
	RingDetection( const Graph &graph) :
		graph_size_( boost::num_vertices(graph) )
	{
		// initialize paths_ with all edges
		Initialize( graph);

		// remove all vertices, update paths, collect rings
		while( !paths_.empty())
		{
			// in the paper it's suggested to remove vertices with smaller amount of edges first
			// such a heuristic should go here instead of taking the 'next' vertex
			Remove( paths_.front().front());
		}
	}

	//! clone the object
	RingDetection *Clone() const
	{
		return new RingDetection( *this);
	}



	/////////////////
	// data access //
	/////////////////

	//! @brief Get the paths between vertices in the graph
	//! @return list of paths
	const std::list< std::vector< size_t> > &GetPaths() const{
		return paths_;
	}

	//! @brief Get all rings in the graph
	//! @return list of rings (all closed paths)
	const utility::vector1<utility::vector1<VD> > GetRings() const{
		utility::vector1<utility::vector1<VD> > rings;
		utility::vector1<VD> vertices;
		for(std::list<std::vector<size_t> >::const_iterator it = rings_.begin(); it != rings_.end(); ++it){
			std::vector<size_t> list_ring = *it;
			for(size_t j=0; j < list_ring.size(); ++j){
				size_t number = list_ring[j];
				VD lwrg_vd = index_to_vd_.find(number)->second;
				vertices.push_back( lwrg_vd );
			}
			rings.push_back(vertices);
		}
		//return rings_;
		return rings;
	}

	/////////////
	// methods //
	/////////////

private:

	//////////////////////
	// helper functions //
	//////////////////////

	//! @brief initialize paths_ with all edges from graph (which is either a GraphWithData or a ConstGraph)

	void Initialize( const Graph &graph) {
		//setup the vertex and edge iterators. This is because the following code uses index based
			//system for everything. Set them up once, dont have to worry anymore
			VIter v_start, v_end;
			boost::tie(v_start, v_end) = boost::vertices(graph);
			size_t index(0);
			for(VIter it = v_start; it != v_end; ++it){
				index_to_vd_[index] = *it;
				vd_to_index_[*it] = index;
				++index;
			}
			// determine the vertex girths; e.g. the shortest cycle containing a particular vertex
			std::vector< std::size_t > shortest_cycles(graph_size_ );
			for(size_t vertex_number=0; vertex_number < graph_size_; ++vertex_number){
				shortest_cycles[vertex_number] = LengthOfSmallestCycleWithVertex(graph, vertex_number);
			}

			// our graphs are not directed // const bool is_directed( false ); // determine whether the graph is undirected
			// initialize paths_
			for(
					size_t source_vertex_number = 0, graph_size = graph_size_;
					source_vertex_number != graph_size;
					++source_vertex_number
			)
			{
				if( shortest_cycles[source_vertex_number] > graph_size) // this vertex is not part of any cycle
				{
					continue;
				}

				//we need to get the indices for all the vertex neighbors. We do this by first getting the
				//outer edge (in an undirected graph, out edge is the same as in edge) which will tell you
				//which is bonded to
				OutEdgeIter oe_start, oe_end;
				boost::tie(oe_start, oe_end) = boost::out_edges(index_to_vd_[source_vertex_number], graph);
				std::vector<size_t> source_vertex_number_neighbors;
				for(OutEdgeIter itr_edge_target = oe_start; itr_edge_target != oe_end; ++itr_edge_target){
					VD target = boost::target(*itr_edge_target, graph); //get the vertex that is the target of edge (boost::source would return original vertex)
					source_vertex_number_neighbors.push_back(vd_to_index_[target]);
				}
				for(utility::vector1<size_t>::iterator itr_edge_target = source_vertex_number_neighbors.begin(); itr_edge_target != source_vertex_number_neighbors.end(); ++itr_edge_target){
					if( shortest_cycles[ *itr_edge_target ] > graph_size) // the edge target vertex is not part of any cycle
					{
						continue;
					}
					if( source_vertex_number <= *itr_edge_target)
					{
						std::vector<size_t> initial_path(2);
						initial_path[0] = source_vertex_number;
						initial_path[1] = *itr_edge_target;
						paths_.push_back(initial_path);
					}
				}
			}
	}

	//! @brief Remove vertices, store new paths and rings
	//! @param vertex vertex to be removed
	void Remove( const std::size_t vertex){
		// BCL_MessageDbg
			//(
			//"Now " + util::Format()( paths_.GetSize()) + " paths and " + util::Format()( rings_.GetSize()) + " rings"
			//);

			// storage for the updated paths
			std::list<std::vector<size_t> > paths_with_vertex;

			// transfer any paths that end in vertex into paths_with_vertex
			for(std::list<std::vector<size_t> >::iterator itr_paths(paths_.begin()),
					itr_paths_end(paths_.end());
					itr_paths != itr_paths_end;
			){
				if( vertex == itr_paths->front() || vertex == itr_paths->back() ){
					std::list<std::vector<size_t> >::iterator old_itr_paths(itr_paths);
					++itr_paths;
					paths_with_vertex.splice(paths_with_vertex.end(), paths_, old_itr_paths);
				}
				else{
					++itr_paths;
				}


			}


			// make an incidence vector for path_a
			// vertex_is_in_path[ x] -> tells whether vertex #x is in the path

			//update paths
			std::vector< bool> vertex_is_in_path( graph_size_);

			//update paths
			for(std::list<std::vector<size_t> >::iterator itr_paths(paths_with_vertex.begin()),
					itr_paths_end(paths_with_vertex.end());
					itr_paths != itr_paths_end;
					++itr_paths)
			{
				const std::size_t path_size( itr_paths->size() - 2);
				SetupAdjacencyVector( vertex_is_in_path, *itr_paths);

				std::list< std::vector<size_t> >::iterator itr_paths_match(itr_paths);
				++itr_paths_match;

				while( itr_paths_match != itr_paths_end)
				{
					if
					(
							path_size + itr_paths_match->size() <= graph_size_ // eliminate obvious overlaps
							&& !Overlap( vertex_is_in_path, *itr_paths_match)
					)
					{
						// combine the paths and add them to paths_
						paths_.insert( paths_.begin(), CombinePaths( vertex, *itr_paths, *itr_paths_match));

						std::list< std::vector< size_t> >::iterator itr_begin( paths_.begin());

						// check whether the new ring was actually a path
						if( itr_begin->front() == itr_begin->back())
						{
							// yep, so splice it into rings instead (this removes it from paths_)
							rings_.splice( rings_.end(), paths_, itr_begin);
						}
					}

					++itr_paths_match;
				}
			}
	}

	//! @brief Combine two paths at their common vertex
	//! @param COMMON_vertex vertex to join paths at
	//! @param PATH_A first path to be combined
	//! @param PATH_B second path to be combined
	//! @return combined path
	std::vector< size_t> CombinePaths
	(
			const std::size_t COMMON_vertex,
			const std::vector< size_t> &path_a,
			const std::vector< size_t> &path_b
	) const {
		std::vector< size_t> new_path( path_b.size() + path_a.size() - 1);

		std::vector< size_t>::iterator new_path_itr( new_path.begin());

		// copy path A such that the common vertex is last
		if( COMMON_vertex == path_a.back())
		{ // copy forward
			new_path_itr = std::copy( path_a.begin(), --path_a.end(), new_path_itr);
		}
		else // COMMON_vertex == path_a.FirstElement()
		{ // copy in reverse
			new_path_itr = std::copy( path_a.rbegin(), --path_a.rend(), new_path_itr);
		}

		// copy path B such that the common vertex is first
		if( COMMON_vertex == path_b.front()) // follow with path B, omitting the common element
		{ //copy forward
			std::copy( path_b.begin(), path_b.end(), new_path_itr);
		}
		else
		{ // copy in reverse
			std::copy( path_b.rbegin(), path_b.rend(), new_path_itr);
		}

		return new_path;
	}

	//! @brief Check if two paths overlap
	//! @param vertex_IS_IN_PATH_A adjacency vector for path a
	//! @param PATH_B second path to be combined
	//! @return if paths overlap
	bool Overlap
	(
			const std::vector< bool> &vertex_is_in_path_a,
			const std::vector< size_t > &path_b
	) const {
		// check all internal elements in path_a for occurrence in path_b
		for
		(
				std::vector< size_t>::const_iterator
				itr_path( path_b.begin() ),
				itr_path_end( path_b.end());
				itr_path != itr_path_end;
				++itr_path
		)
		{
			if( vertex_is_in_path_a[ *itr_path])
			{
				return true;
			}
		}

		return false;
	}

	//! @brief set the adjacency vector with a given path
	//! @param vertex_IS_IN_PATH the vector to setup such that vertex_IS_IN_PATH[ x] == true iff x is in PATH
	//! @param PATH the path to use in setting up the adjacency vector
	void SetupAdjacencyVector
	(
			std::vector< bool> &vertex_is_in_path,
			const std::vector< size_t> &path
	) const {
		vertex_is_in_path.assign( graph_size_, false);
		for
		(
				std::vector< size_t>::const_iterator itr_path( ++path.begin()), itr_path_end( --path.end());
				itr_path != itr_path_end;
				++itr_path
		)
		{
			vertex_is_in_path[ *itr_path] = true;
		}
	}


	//! @brief LengthOfSmallestCycleWithVertex finds the length of the shortest cycle beginning at vertex
	//! @param graph the graph in which the vertex resides
	//! @param vertex index of the vertex to find the girth from
	//! @param CAN_VISIT 0 if the associated vertex cannot be visited (normally vertices that are already known to be non-cyclical), non-zero if it can
	//! @return length of the shortest cycle beginning at vertex (may be undefined)
	//! @note any vertices that are unreachable from the vertex at vertex are skipped
	//! @note this works appropriately on both directed and undirected graphs
	size_t LengthOfSmallestCycleWithVertex
	(
			const Graph &graph,
			const size_t &vertex,
			const std::vector< size_t> CAN_VISIT = std::vector< size_t>()
	) {
		size_t shortest_branch( std::numeric_limits< size_t>::max());

		// if there are only 0 or 1 vertices, then there can be no cycle, so return immediately
		if( boost::in_degree(vertex, graph) <= 1)
		{
			return shortest_branch;
		}

		const size_t unseen_flag( std::numeric_limits< size_t>::max() );
		const size_t size( graph_size_);

		// distances will hold the distance of each vertex to the current vertex
		std::vector< size_t> distances( size, unseen_flag);

		// initialize seen_vertices with the first vertex
		std::vector< size_t> seen_vertices_queue( 1, vertex);

		// store the first branch # (index of the vertex connected to vertex which first reached this particular vertex
		// If multiple branches reached a particular vertex simultaneously, then they are pushed back into the list
		std::vector< size_t> branch_number( size, unseen_flag);

		// allocate enough memory for the queue to contain all vertices in the graph
		seen_vertices_queue.reserve(size);

		size_t vertex_queue_position( 0); // index of the active vertex in the breadth-first-search

		distances[vertex] = 0; // vertex is where the search starts, so its distance is 0
		branch_number[vertex] = vertex;

		bool cycle_found( false);

		size_t distance( 1);    // initial distance will be 0

		std::vector< size_t> can_visit;
		if( CAN_VISIT.size() != size)
		{
			can_visit = std::vector< size_t>( size, size_t( 1));
		}
		const std::vector< size_t> &can_visit_ref( CAN_VISIT.size() != size ? can_visit : CAN_VISIT);

		{

			//we need to get the indices for all the vertex neighbors. We do this by first getting the
			//outer edge (in an undirected graph, out edge is the same as in edge) which will tell you
			//which is bonded to
			OutEdgeIter oe_start, oe_end;
			boost::tie(oe_start, oe_end) = boost::out_edges(index_to_vd_[vertex], graph);
			std::vector<size_t> row;
			for(OutEdgeIter itr_edge_target = oe_start; itr_edge_target != oe_end; ++itr_edge_target)
			{
				VD target = boost::target(*itr_edge_target, graph);
				row.push_back(vd_to_index_[target]);
			}

			//I think I did this wrong. I should be looking at the vertex not the edges. row[i] returns
			//an edge. I think that it should return a vertex
			for(size_t i(0), number_seen(boost::in_degree(vertex, graph)); i < number_seen; i++){
				// assign the 1st-degree nieghbors of a branch number = their index
				if( can_visit_ref[  row[i] ] ){
					distances[ row[i]  ] = distance;
					branch_number[ row[i]  ] = row[i];
					seen_vertices_queue.push_back( row[i] );
				}
				++vertex_queue_position;
				++distance;
			}

		}

		// if there are only 0 or 1 visitable vertices, then there can be no cycle, so return immediately
		if(  seen_vertices_queue.size() <= 1)
		{
			return shortest_branch;
		}

		// So long as there are vertices left in the queue whose connections haven't been examined, this loop will continue
		// unless all vertices are put into the graph.
		while( vertex_queue_position < seen_vertices_queue.size() && !cycle_found)
		{
			// loop over all vertices left in the queue that are at the current distance.  If they connect to any vertices
			// not already in the queue, then add them to the queue and record their distance.
			// Stop if all vertices in the graph are in the queue or we reach the last vertex at this distance
			for
			(
					const std::size_t last_vertex_at_distance( seen_vertices_queue.size());
					vertex_queue_position < last_vertex_at_distance;
					vertex_queue_position++
			)
			{
				const std::size_t current_vertex( seen_vertices_queue[ vertex_queue_position]);
				// target row is a reference to the edges reachable from the current vertex
				OutEdgeIter oe_start, oe_end;
				boost::tie(oe_start, oe_end) = boost::out_edges(index_to_vd_[vertex], graph);
				std::vector<size_t> target_row;
				for(OutEdgeIter itr_edge_target = oe_start; itr_edge_target != oe_end; ++itr_edge_target)
				{
					VD target = boost::target(*itr_edge_target, graph);
					target_row.push_back(vd_to_index_[target]);
				}

				for( std::size_t i( 0), number_seen( target_row.size()); i < number_seen; i++)
				{
					const std::size_t new_vertex( target_row[ i]);
					if( !can_visit_ref[ new_vertex])
					{
						// do nothing; can't visit the associated vertex
					}
					else if( distances[ new_vertex] == unseen_flag) // found a vertex in the target list of vertex seen_vertices_queue(vertex_queue_position)
					{
						seen_vertices_queue.push_back(new_vertex);
						distances[ new_vertex] = distance;
						branch_number[ new_vertex] = branch_number[ current_vertex];
					}
					else if // check for cycles
					(
							new_vertex != vertex
							&& branch_number[current_vertex] != branch_number[ new_vertex]
					)
					{
						cycle_found = true;
						shortest_branch = std::min( shortest_branch, distance + distances[ new_vertex]);
					}
				}
			}

			distance++;
		}

		return shortest_branch;
	}


private:

	//////////
	// data //
	//////////

	//! paths between vertices in the graph
	std::list< std::vector< size_t> > paths_;

	//! rings in the graph
	std::list< std::vector< size_t> > rings_;

	//! size of the graph
	size_t graph_size_;

	boost::unordered_map<size_t, VD> index_to_vd_;

	boost::unordered_map<VD, size_t> vd_to_index_;

};


} // graph
} // utility



#endif
