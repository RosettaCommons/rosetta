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
#include <core/chemical/Bond.fwd.hh>
#include <core/types.hh>
#include <boost/graph/adjacency_list.hpp>
#include <core/chemical/ResidueGraphTypes.hh>
#include <vector>


namespace utility {
namespace graph {


/// @brief basic chemical Bond
///
/// @details name, element, certain properties and parameters from .params file
///
class RingDetection {

public:

	//////////////////////////////////
	// construction and destruction //
	//////////////////////////////////

	//! @brief default constructor
	RingDetection();

	//! @brief constructor from a graph (either GraphWithData or ConstGraph)
	//! @param graph graph for exhaustive ring detection
	//! @note prefer using ConstGraph here, it is often many times faster than GraphWithData
	RingDetection( const core::chemical::LightWeightResidueGraph &graph);

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
	const utility::vector1<utility::vector1<core::chemical::lwrg_VD> > GetRings() const{
		utility::vector1<utility::vector1<core::chemical::lwrg_VD> > rings;
		utility::vector1<core::chemical::lwrg_VD> vertices;
		for(std::list<std::vector<size_t> >::const_iterator it = rings_.begin(); it != rings_.end(); ++it){
			std::vector<size_t> list_ring = *it;
			for(size_t j=0; j < list_ring.size(); ++j){
				size_t number = list_ring[j];
				core::chemical::lwrg_VD lwrg_vd = index_to_vd_.find(number)->second;
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

	void Initialize( const core::chemical::LightWeightResidueGraph &graph);

	//! @brief Remove vertices, store new paths and rings
	//! @param vertex vertex to be removed
	void Remove( const size_t vertex);

	//! @brief Combine two paths at their common vertex
	//! @param COMMON_vertex vertex to join paths at
	//! @param PATH_A first path to be combined
	//! @param PATH_B second path to be combined
	//! @return combined path
	std::vector< size_t> CombinePaths
	(
			const size_t COMMON_vertex,
			const std::vector< size_t> &PATH_A,
			const std::vector< size_t> &PATH_B
	) const;

	//! @brief Check if two paths overlap
	//! @param vertex_IS_IN_PATH_A adjacency vector for path a
	//! @param PATH_B second path to be combined
	//! @return if paths overlap
	bool Overlap
	(
			const std::vector< bool> &vertex_IS_IN_PATH_A,
			const std::vector< size_t > &PATH_B
	) const;

	//! @brief set the adjacency vector with a given path
	//! @param vertex_IS_IN_PATH the vector to setup such that vertex_IS_IN_PATH[ x] == true iff x is in PATH
	//! @param PATH the path to use in setting up the adjacency vector
	void SetupAdjacencyVector
	(
			std::vector< bool> &vertex_IS_IN_PATH,
			const std::vector< size_t> &PATH
	) const;


	//! @brief LengthOfSmallestCycleWithVertex finds the length of the shortest cycle beginning at vertex
	//! @param graph the graph in which the vertex resides
	//! @param vertex index of the vertex to find the girth from
	//! @param CAN_VISIT 0 if the associated vertex cannot be visited (normally vertices that are already known to be non-cyclical), non-zero if it can
	//! @return length of the shortest cycle beginning at vertex (may be undefined)
	//! @note any vertices that are unreachable from the vertex at vertex are skipped
	//! @note this works appropriately on both directed and undirected graphs
	size_t LengthOfSmallestCycleWithVertex
	(
			const core::chemical::LightWeightResidueGraph &graph,
			const size_t &vertex,
			const std::vector< size_t> CAN_VISIT = std::vector< size_t>()
	);

	//////////////////////
	// input and output //
	//////////////////////

protected:

	//! @brief read from std::istream
	//! @param ISTREAM input stream
	//! @return istream which was read from
	// virtual std::istream &Read( std::istream &ISTREAM);

	//! @brief write to std::ostream
	//! @param OSTREAM output stream to write to
	//! @param INDENT number of indentations
	//! @return output stream which was written to
	//virtual std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

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

	boost::unordered_map<size_t, core::chemical::lwrg_VD> index_to_vd_;

	boost::unordered_map<core::chemical::lwrg_VD, size_t> vd_to_index_;

};


} // graph
} // utility



#endif // INCLUDED_core_chemical_Bond_HH
