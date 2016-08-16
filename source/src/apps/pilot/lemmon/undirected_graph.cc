// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/** @page readAndWrite
	This simple file describes reading and writing a PDB with Rosetta
	Run it like this:
	"readAndWritePDB.cc -in::file::s <PDB file name> -in::path::database <database root dir>"

	If you provide a PDB with a ligand that has multiple residues, this code will "glue" them together
	Of course you would need to provide information about these additional residues:
		in::file::extra_res_fa <list of extra params files, one per residue type>
 */

/// @file   apps/pilot/lemmon/readAndWrite.cc
///
/// @brief This is to illustrate reading a PDB
/// @detail Run this script with the following arguments:
/// 1) in::file::s <list of one PDB to score>
/// 2) in::path::database <list of one database root directory>
/// 3) in::file::extra_res_fa <list of extra params files, one per residue type>
/// @author Gordon Lemmon (glemmon@gmail.com)

#include <boost/graph/undirected_graph.hpp>
#include <boost/graph/property_iter_range.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/named_function_params.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
//#include <boost/graph/graph_utility.hpp> // The print_graph function here is lame so I write my own
#include <iostream>
#include <utility/vector1.hh>

struct Atom{
	int i;
	char a;
	bool b;

	Atom():i(1),a('a'),b(true){};
	Atom(int i, char a, bool b):i(i),a(a),b(b){};
	Atom(Atom const & atom):i(atom.i),a(atom.a),b(atom.b){}
};

struct Bond{
	int i;
	char a;
	bool b;

	Bond():i(4),a('d'),b(false){};
	Bond(int i, char a, bool b):i(i),a(a),b(b){};
	Bond(Bond const & bond):i(bond.i),a(bond.a),b(bond.b){}
};

std::ostream & operator<<(std::ostream& os, Atom a){os << a.i << " " << a.a << " " << a.b; return os; }
std::ostream & operator<<(std::ostream& os, Bond b){os << b.i << " " << b.a << " " << b.b; return os; }

///////////////////////////////////////////////////////////////////

typedef boost::undirected_graph<
		Atom, // struct with properties of a node
		Bond // struct with properties of an edge
		/*,ResidueType*/
> Graph;

typedef Graph::vertex_descriptor VD;
typedef utility::vector1<VD> VDs;
typedef Graph::edge_descriptor ED;
typedef std::pair<ED, bool> EdgeBoolPair;
typedef boost::graph_traits<Graph>::vertex_iterator VIter;
typedef boost::graph_traits<Graph>::edge_iterator EIter;
typedef std::pair<VIter, VIter> VIterPair;
typedef std::pair<EIter, EIter> EIterPair;
typedef boost::graph_traits<Graph>::out_edge_iterator OutEdgeIter;
typedef boost::graph_traits<Graph>::in_edge_iterator InEdgeIter;
typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIter;
typedef std::pair<OutEdgeIter, OutEdgeIter> OutEdgeIterPair;
typedef std::pair<InEdgeIter, InEdgeIter> InEdgeIterPair;
typedef std::pair<AdjacencyIter, AdjacencyIter> AdjacencyIterPair;

// This is used as a predicate to make a filtered graph
class AtomNumFilter{
public:
	AtomNumFilter(){}
	AtomNumFilter(Graph& graph):graph_(&graph){}
	bool operator()(VD const vd) const{
		return (*graph_)[vd].i < 3;
	};
private:
	Graph * graph_; // Cannot use a reference because 0-arg constructor needed by boost::iterators
};

// This "visitor" is used with breadth first search
class PathFinder:public boost::default_bfs_visitor {
public:
	PathFinder(VD vd_to_find, VDs & path):
			found_(false), vd_to_find_(vd_to_find), path_(path) {
		assert(path_.empty());// We will fill it
	}
	void discover_vertex(VD vd, Graph const & /*g*/)
	{
		if( found_ ) return; // short circuit the BFS
		if( vd == vd_to_find_){
			found_ = true;
			calculate_path();
		}
	}
	void tree_edge(ED ed, Graph const & g){ // invoked on edges as they become part of the search tree
		parents_[ boost::target(ed, g)] = boost::source(ed,g);
	}
private:
	void calculate_path(){
		VD current = vd_to_find_;
		path_.push_back(current);
		while( parents_.find( current ) != parents_.end()){
			path_.push_back( parents_[current] );
			current = parents_[current];
		}
	}
	bool found_;
	VD vd_to_find_;
	std::map<VD, VD> parents_; // map of child to parent
	VDs & path_;
};

/// A better print function that works with any graph type
template< typename graph_t >
void print_graph(graph_t const & g){
	typedef typename graph_t::vertex_descriptor vd_t;
	typedef typename boost::graph_traits<graph_t>::vertex_iterator v_iter_t;

	std::cout << "Iterate over all vertices" << std::endl;
	for( std::pair<v_iter_t, v_iter_t> vp = boost::vertices(g); vp.first != vp.second; ++vp.first){
		v_iter_t v_iter= vp.first;
		vd_t vd = *v_iter;
		Atom a = g[vd];
		std::cout << vd << " " << a << std::endl;
	}

	std::cout << std::endl;

	typedef typename graph_t::edge_descriptor ed_t;
	typedef typename boost::graph_traits<graph_t>::edge_iterator e_iter_t;
	std::cout << "Iterate over all edges" << std::endl;
	for( std::pair<e_iter_t, e_iter_t> ep = boost::edges(g); ep.first != ep.second; ++ep.first){
		e_iter_t e_iter = ep.first;
		ed_t ed = *e_iter;
		assert( boost::source(ed, g) ); // source returns a vertex_descriptor
		assert( boost::target(ed, g) ); // target returns a vertex_descriptor
		Bond b = g[ed];
		std::cout << ed << " " << b << std::endl;
	}
	std::cout << std::endl;
}

void manual_copy(Graph const & g_old, Graph & g_new){
	// Manual Copy
	{
		std::map<VD, VD> old_to_new;
		for( VIterPair vp = boost::vertices(g_old); vp.first != vp.second; ++vp.first){
			VIter v_iter= vp.first;
			VD v_old = *v_iter;
			Atom a = g_old[v_old];
			VD v_new = g_new.add_vertex(a);
			old_to_new[v_old] = v_new;
		}
		for( EIterPair ep = boost::edges(g_old); ep.first != ep.second; ++ep.first){
			EIter e_iter = ep.first;
			ED ed = *e_iter;
			Bond b = g_old[ed];
			VD source = old_to_new[ boost::source(ed, g_old) ]; /// Todo replace with safe find function
			VD target = old_to_new[ boost::target(ed, g_old) ]; /// Todo replace with safe find function
			assert(source);
			assert(target);
			g_new.add_edge( source, target, b);
		}
		std::cout << "g_new" << std::endl;
		print_graph(g_new);
	}
}

//////////////////////////////////////////////////////////////////////////////
int main()
{
	Graph g;

	std::cout << "Setup a graph with 3 vertices and two edges" << std::endl;
	VD v1 = g.add_vertex( Atom() );
	VD v2 = g.add_vertex( Atom(2,'b',false) );
	VD v3 = g.add_vertex( Atom(3,'c',true) );
	assert(v1 && v2 && v3);
	Graph::vertex_index_type v1_index = boost::get_vertex_index(v1, g);
	Graph::vertex_index_type v2_index = boost::get_vertex_index(v2, g);
	Graph::vertex_index_type v3_index = boost::get_vertex_index(v3, g);

	std::cout << "num_vertices: " << g.num_vertices() << ", num_edges: " << g.num_edges() << std::endl;

	EdgeBoolPair ebp1= g.add_edge(v1, v2, Bond() );
	EdgeBoolPair ebp2= g.add_edge(v2, v3, Bond(5,'e',true) );
	assert(ebp1.second);
	assert(ebp2.second);
	ED e1 = ebp1.first;
	ED e2 = ebp2.first;

	Graph::edge_index_type e1_index = boost::get_edge_index(e1, g);
	Graph::edge_index_type e2_index = boost::get_edge_index(e2, g);

	std::cout << "Get the Atoms and Bonds from this graph using descriptors" << std::endl;
	Atom a1 = g[v1];
	Atom a2 = g[v2];
	Atom a3 = g[v3];
	Bond b1 = g[e1];
	Bond b2 = g[e2];

	std::cout << v1 << " " << v1_index << " " << a1 << std::endl;
	std::cout << v2 << " " << v2_index << " " << a2 << std::endl;
	std::cout << v3 << " " << v3_index << " " << a3 << std::endl;
	std::cout << e1 << " " << e1_index << " " << b1 << std::endl;
	std::cout << e2 << " " << e2_index << " " << b2 << std::endl;
	std::cout << std::endl;

	std::cout << "Access the vertices connected to an edge"  << std::endl;
	VD source = boost::source(e1, g);
	VD target = boost::target(e1, g);
	std::cout << "source: "<< g[source] << "; target: " << g[target] << std::endl;

	std::cout << std::endl;
	std::cout << "g" << std::endl;
	print_graph(g);
	std::cout << std::endl;

	// Copying the graph...
	{
		Graph g2 = g;
		std::cout << "g2" << std::endl;
		print_graph(g2);
		Graph g3(g);
		std::cout << "g3" << std::endl;
		print_graph(g3);
		Graph g4;
		manual_copy(g, g4);
		std::cout << "g4" << std::endl;
		print_graph(g4);
	}


	/// In the case of an undirected graph, out and in edges are identical
	std::cout << "Iterate over out edges" << std::endl;
	for(OutEdgeIterPair ep = boost::out_edges(v1, g); ep.first != ep.second; ++ep.first){
		OutEdgeIter e_iter= ep.first;
		ED ed = *e_iter;
		Bond b = g[ed];
		std::cout << ed << " " << b << std::endl;
	}
	std::cout << std::endl;

	std::cout << "Iterate over in edges" << std::endl;
	for( InEdgeIterPair ep = boost::in_edges(v1, g); ep.first != ep.second; ++ep.first){
		InEdgeIter e_iter = ep.first;
		ED ed = *e_iter;
		Bond b = g[ed];
		std::cout << ed << " " << b << std::endl;
	}
	std::cout << std::endl;

	std::cout << "Iterate over adjacent vertices" << std::endl;
	for( AdjacencyIterPair ep = boost::adjacent_vertices(v2, g); ep.first != ep.second; ++ep.first){
		AdjacencyIter adj_iter = ep.first;
		VD vd = *adj_iter;
		Atom a = g[vd];
		std::cout << vd << " " << a << std::endl;
	}
	std::cout << std::endl;

	std::cout << "Iterate over an atom property with property_map & graph access" << std::endl;
	typedef boost::property_map<Graph, int Atom::*>::type AtomIntMap; // This works
	AtomIntMap atom_int_map = boost::get(&Atom::i, g);
	for( VIterPair vp = boost::vertices(g); vp.first != vp.second; ++vp.first){
		VIter v_iter= vp.first;
		VD vd = *v_iter;
		int atom_int = atom_int_map[vd];
		std::cout << vd << " " << atom_int << std::endl;
	}
	std::cout << std::endl;

	///Figure out how to iterate over a map without a graph
//	std::cout << "Iterate over a vector of atom properties made from a property map (Doesn't work this way)" << std::endl;
//	std::vector<int> atom_ints(boost::num_vertices(g));
//	boost::make_property_map_iterator(atom_int_map, atom_ints.begin()); // This doesn't do anything!
//	for(
//			std::vector<int>::iterator begin = atom_ints.begin(), end = atom_ints.end();
//			begin != end;
//			++begin
//	){
//		std::cout << "Int: " << *begin << std::endl;
//	}
//	std::cout << std::endl;

	////////////////////////////////////////////////////////////////////////////////////
	std::cout << "Find shortest unweighted path. Use BFS and a custom visitor" << std::endl;
	VDs path; // In order from v3 to v1 (child to parent)
	PathFinder path_finder(v3, path);
	boost::breadth_first_search(g, v1, visitor(path_finder));
	for(
			VDs::iterator iter = path.begin(),
			end= path.end();
			iter != end;
			++iter
	){
		VD vd = *iter;
		std::cout << vd << " ";
		// Mark these as "mainchain" in ResidueType
	}
	std::cout << std::endl << std::endl;;

	//////////////////////////////////////////////////////////////////////////////
	std::cout << "Create a filtered graph" << std::endl;
	typedef boost::filtered_graph<Graph, boost::keep_all, AtomNumFilter> AtomNumGraph;
	AtomNumFilter atom_num_filter(g);
	AtomNumGraph filtered_graph(g, boost::keep_all(), atom_num_filter); // holds a reference to g so it always returns what you need.
	print_graph(filtered_graph);

	/// Removing graph vertices MUST be preceeded by first removing their edges with clear_vertex
	/// Otherwise serious problems may occur
	std::cout << "Delete a vertex and print the removed vertex, demonstrating it still exists!" << std::endl;
	g.clear_vertex(v1);
	g.remove_vertex(v1);
	std::cout << g[v1] << std::endl;
	print_graph(g);
	std::cout << std::endl;

	std::cout << "The filtered graph stays up to date!" << std::endl;
	print_graph(filtered_graph);
	std::cout << std::endl;

	std::cout << "Restore the deleted vertex and edge!" << std::endl;
	VD new_v1 = g.add_vertex(g[v1]);
	g.add_edge(new_v1, v2, Bond());
	print_graph(g);
	std::cout << "The filtered_graph stays up to date!" << std::endl;
	print_graph(filtered_graph);

	/////// Find the shortest weighted path. Use Dijkstras...
//	std::vector<VD> parents(boost::num_vertices(g));
//	std::vector<int> distances(boost::num_vertices(g));
//	boost::predecessor_map(&parents[0]).distance_map(&distances[0]);
//	boost::dijkstra_shortest_paths(g, v1, boost::predecessor_map(&parents[0]).distance_map(&distances[0]));
//
//	std::cout << "distances and parents:" << std::endl;
//	boost::graph_traits < Graph >::vertex_iterator vertexIterator, vend;
//	std::cout << "distance(" << v3 << ") = " << distances[ v3 ] << ", ";
//	std::cout << "parent(" << v3 << ") = " << parents[ v3 ] << std::endl;

}
